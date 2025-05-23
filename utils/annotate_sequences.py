import os
import subprocess
import gzip
import shutil
import concurrent.futures
from pathlib import Path
from Bio import SeqIO

# ==== GLOBAL SETTINGS ====
DIAMOND_BIN = os.environ.get("DIAMOND_BIN", "diamond")
TRNASCAN_BIN = os.environ.get("TRNASCAN_BIN", "/usr/local/tRNAscan-SE-2.0.12/bin/tRNAscan-SE")
MINIPROT_BIN = os.environ.get("MINIPROT_BIN", "miniprot")
AHRD_BIN = os.environ.get("AHRD_BIN", "/home/diriano/sugarcaneMapperIDs/AHRD/dist/ahrd.jar")
DBCAN_ENV_NAME = os.environ.get("DBCAN_ENV", "CAZyme_annotation_main")
### Set the path to the databases
DIAMOND_DATABASES = ["DBs/uniprot_sprot.fasta","DBs/uniprot_trembl.fasta"]
DBCAN_DATABASES_PATH = os.environ.get("DBCAN_DATABASES_PATH", "/home/diriano/sugarcaneMapperIDs/DBs/DB_CAN/")
# ==== END GLOBAL SETTINGS ====


def extract_genotype_from_filename(filename):
    filename_prefix = 'sugarcanePanTranscriptome_06052025_'
    filename_suffix = '_transcript_4120596.fix.fasta.gz'
    part = filename.replace(filename_prefix, "").replace(filename_suffix, "")
    return part

def check_exec_path(name, path):
    if shutil.which(path) is None:
        raise FileNotFoundError(f"Executable for {name} not found at '{path}'. Check PATH or environment variable.")
    
def check_conda_env_exists(env_name):
    try:
        result = subprocess.run(["conda", "env", "list"], capture_output=True, text=True, check=True)
        if env_name not in result.stdout:
            raise EnvironmentError(f"Conda environment '{env_name}' not found.")
        print(f" Conda environment '{env_name}' found.")
    except FileNotFoundError:
        raise EnvironmentError("Conda is not installed or not in PATH.")
    except subprocess.CalledProcessError as e:
        raise EnvironmentError(f"Failed to list conda environments: {e}")

def modify_transcript_headers(transcripts_fasta_gz, genotype, output_fasta):
    organism_str = f"[moltype=transcribed_RNA] [tech=TSA] [organism=Saccharum hybrid cultivar {genotype}]"
    with gzip.open(transcripts_fasta_gz, "rt") as input_handle, open(output_fasta, "w") as output_handle:
        for record in SeqIO.parse(input_handle, "fasta"):
            record.description = f"{record.id} {organism_str}"
            SeqIO.write(record, output_handle, "fasta")
    return output_fasta

def ensure_diamond_db(fasta_path):
    db_path = Path(fasta_path).with_suffix(".dmnd")
    print(f'  Checking DIAMOND database: {db_path}')
    if db_path.exists():
        print(f"  DIAMOND database already exists: {db_path}")
        return db_path
    print(f"  Creating DIAMOND database for {fasta_path}...")
    subprocess.run([DIAMOND_BIN, "makedb", "--in", str(fasta_path), "--db", Path(fasta_path).with_suffix(".dmnd")], check=True)
    return db_path

def run_diamond(query_fasta, db_paths, output_prefix, threads=20):
    for db in db_paths:
        db_name = Path(db).stem
        output_file = Path(f"{output_prefix}.blast.{db_name}.txt")
        if output_file.exists():
            print(f"  Skipping DIAMOND for {db_name} — output exists: {output_file}")
            continue
        print(f"  writing in {output_file}")
        ensure_diamond_db(db)
        print(f"  Running DIAMOND: {query_fasta} vs {db_name} with {threads} threads")
        subprocess.run([
            DIAMOND_BIN, "blastp",
            "--query", query_fasta,
            "--db", str(db).replace(".fasta", ".dmnd"),
            "--out", str(output_file),
            "--outfmt", "6",
            "--evalue", "1e-5",
            "--threads", str(threads),
            "--quiet",
            "--sensitive"
        ], check=True)

def split_fasta(input_fasta, output_prefix, chunk_size=500):
    chunk_files = []
    records = list(SeqIO.parse(input_fasta, "fasta"))
    for i in range(0, len(records), chunk_size):
        chunk_records = records[i:i + chunk_size]
        chunk_file = f"{output_prefix}.chunk{i // chunk_size + 1}.fasta"
        with open(chunk_file, "w") as f:
            SeqIO.write(chunk_records, f, "fasta")
        chunk_files.append(chunk_file)
    return chunk_files

def trnascan_chunk_worker(chunk_file, output_file, stats_file):
    if Path(output_file).exists() and Path(stats_file).exists():
        print(f"  Skipping tRNAscan chunk (already exists): {output_file}")
        return output_file, stats_file
    print(f"  Running tRNAscan-SE on {chunk_file}")
    cmd = [
        TRNASCAN_BIN,
        "-E", "--thread", "10",
        "-o", output_file,
        "-m", stats_file,
        chunk_file
    ]
    subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    return output_file, stats_file

def run_trnascan_chunked_parallel(transcript_fasta, output_prefix, chunk_size=500, max_workers=4):
    merged_txt = f"{output_prefix}.trnascan.txt"
    merged_stats = f"{output_prefix}.trnascan.stats.txt"
    if Path(merged_txt).exists() and Path(merged_stats).exists():
        print(f"  Skipping tRNAscan-SE (merged outputs already exist): {merged_txt}")
        return
    chunk_files = split_fasta(transcript_fasta, output_prefix, chunk_size)
    jobs = []
    trna_txt_chunks = []
    stats_chunks = []

    for i, chunk in enumerate(chunk_files):
        chunk_id = i + 1
        out_txt = f"{output_prefix}.chunk{chunk_id}.trnascan.txt"
        stats_txt = f"{output_prefix}.chunk{chunk_id}.trnascan.stats.txt"
        trna_txt_chunks.append(out_txt)
        stats_chunks.append(stats_txt)
        jobs.append((chunk, out_txt, stats_txt))

    print(f"  Submitting {len(jobs)} tRNAscan chunk jobs for {Path(transcript_fasta).name}...")
    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(trnascan_chunk_worker, *job) for job in jobs]
        for future in concurrent.futures.as_completed(futures):
            try:
                future.result()
            except Exception as e:
                print(f"    ❌ Error in tRNAscan chunk: {e}")

    print(f"  Merging tRNAscan chunks into: {merged_txt} and {merged_stats}")
    with open(merged_txt, "w") as out:
        for f in trna_txt_chunks:
            with open(f) as part:
                shutil.copyfileobj(part, out)
            Path(f).unlink()
    with open(merged_stats, "w") as out:
        for f in stats_chunks:
            with open(f) as part:
                shutil.copyfileobj(part, out)
            Path(f).unlink()
    for f in chunk_files:
        Path(f).unlink()

def run_miniprot(proteins_fasta, transcripts_fasta, output_file):
    if Path(output_file).exists():
        print(f"  Skipping Miniprot (already exists): {output_file}")
        return
    print(f"  Running Miniprot...")
    with open(output_file, "w") as out:
        subprocess.run([MINIPROT_BIN, "-t", "15", transcripts_fasta, proteins_fasta], stdout=out, check=True)

def run_ahrd(protein_fasta, output_dir, output_prefix):
    config_path = Path(output_dir) / f"{output_prefix}.ahrd.yml"
    output_csv = Path(output_dir) / f"{output_prefix}.ahrd.csv"
    swiss_blast = Path(output_dir) / f"{output_prefix}.blast.uniprot_sprot.txt"
    trembl_blast = Path(output_dir) / f"{output_prefix}.blast.uniprot_trembl.txt"

    if Path(output_csv).exists():
        print(f"  Skipping AHRD (already exists)")
        return
    
    if not Path(swiss_blast).exists() or not Path(trembl_blast).exists():
        print(f" AHRD will not run because similarity search files are not present. Check that Diamond results are present!")
        return
    
    # Caminhos dos arquivos fixos
    go_annotation = "DBs/goa_uniprot_all.gaf"
    blacklist = "DBs/ahrd/blacklist_descline.txt"
    filter_sprot = "DBs/ahrd/filter_descline_sprot.txt"
    filter_trembl = "DBs/ahrd/filter_descline_trembl.txt"
    token_blacklist = "DBs/ahrd/blacklist_token.txt"

    config_yaml = f"""\
proteins_fasta: {protein_fasta}
gene_ontology_result: {go_annotation}
reference_go_regex: ^UniProtKB\\s+(?<shortAccession>\\S+)\\s+\\S+\\s+(?<goTerm>GO:\\d{7})
prefer_reference_with_go_annos: false
token_score_bit_score_weight: 0.468
token_score_database_score_weight: 0.2098
token_score_overlap_score_weight: 0.3221
output: {output_csv}
blast_dbs:
  swissprot:
    weight: 653
    description_score_bit_score_weight: 2.717061
    file: {swiss_blast}
    database: DBs/uniprot_sprot.fasta
    blacklist: {blacklist}
    filter: {filter_sprot}
    token_blacklist: {token_blacklist}
  trembl:
    weight: 904
    description_score_bit_score_weight: 2.590211
    file: {trembl_blast}
    database: DBs/uniprot_trembl.fasta
    blacklist: {blacklist}
    filter: {filter_trembl}
    token_blacklist: {token_blacklist}
"""

    with open(config_path, "w") as f:
        f.write(config_yaml)

    print(f"  Running AHRD for {protein_fasta}")
    subprocess.run(["java", "-jar", '-Xmx10g','-XX:ActiveProcessorCount=10', AHRD_BIN, str(config_path)], check=True)

def run_dbcan(input_fasta, output_dir, mode="protein"):
    check_conda_env_exists(DBCAN_ENV_NAME)
    output_dir = Path(output_dir)
    if output_dir.exists():
        print(f"  Skipping run_dbcan: output directory already exists: {output_dir}")
        return

    print(f"  Running run_dbcan on {input_fasta}")

    # Construir o comando com ativação do conda
    cmd = f"""
    source $(conda info --base)/etc/profile.d/conda.sh && \
    conda activate {DBCAN_ENV_NAME} && \
    run_dbcan CAZyme_annotation --input_raw_data {input_fasta} --threads 60 --mode {mode} --output_dir {output_dir} --db_dir {DBCAN_DATABASES_PATH}
    """

    subprocess.run(cmd, shell=True, executable="/bin/bash", check=True)

def parse_miniprot_results(miniprot_file):
    res={}
    with open(miniprot_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if fields[0] == fields[5] and int(fields[2]) == 0 and int(fields[3]) == int(fields[1]) and int(fields[9]) // 3 == int(fields[1]):
                if fields[0] not in res:
                    res[fields[0]] = []
                res[fields[0]].append((int(fields[8]), int(fields[7]), fields[4]))

                # print(f"  Miniprot match: {line.strip()}")
    return res

def parse_ahrd_results(ahrd_file):
    res = {}
    with open(ahrd_file) as f:
        for line in f:
            if line.startswith("#") or line.startswith("Protein-Accession"):
                continue
            fields = line.strip().split("\t")
            if len(fields) <= 3:
                continue
            if fields[2]=='***':
                res[fields[0]]=fields[3]
    return res

def parse_trnascan_results(trnascan_file):
    res = {}
    with open(trnascan_file) as f:
        for line in f:
            if line.startswith("#") or line.startswith("Sequence    ") or line.startswith("Name    ") or line.startswith("------"):
                continue
            fields = line.strip().split()
            type='normal'
            strand='+'
            if len(fields) < 8:
                continue
            if len(fields) > 9:
                if fields[9] == 'pseudo':
                    type='pseudogene'
            if fields[2]>fields[3]:
                strand='-'
            if fields[0] not in res:
                res[fields[0]] = []
            res[fields[0]].append((type, strand, int(fields[2]), int(fields[3]), fields[4], fields[5]))
    return res

def parse_results_and_generate_tbl(genotype, transcript_file, ahrd_file, paf_file, trnascan_file, output_tbl, new_transcript_file, id_map_file):
    if Path(output_tbl).exists():
        print(f"  Skipping TBL generation (already exists): {output_tbl}")
        return

    protein2pos = parse_miniprot_results(paf_file)
    ahrd_desc = parse_ahrd_results(ahrd_file)
    trnascan_genes= parse_trnascan_results(trnascan_file)
    new_records = []
    id_map = []

    with open(output_tbl, "w") as out:
        for idx, record in enumerate(SeqIO.parse(transcript_file, "fasta"), start=1):
            old_id = record.id
            new_id = f"{genotype}_{idx:08d}"
            seq_len = len(record.seq)
            modified_seq = record.seq  # default

            out.write(f">Feature {new_id}\n")
            out.write(f"1\t{seq_len}\tgene\n")

            if old_id in protein2pos:
                # Detect if any match is in reverse strand
                for start, end, strand in protein2pos[old_id]:
                    if strand == "-":
                        # print(f"  Reverse strand detected for {old_id}, computing reverse complement.")
                        modified_seq = record.seq.reverse_complement()
                        # Adjust coordinates
                        adj_start = seq_len - end + 1
                        adj_end = seq_len - start + 1
                        start, end = sorted((adj_start, adj_end))
                    else:
                        start, end = sorted((start, end))
                    
                    out.write(f"{start}\t{end}\tCDS\n")
                    if old_id in ahrd_desc:
                        out.write(f"\t\t\tproduct\t{ahrd_desc[old_id]}\n")
                    else:
                        out.write(f"\t\t\tnote\thypothetical protein\n")
            elif old_id in trnascan_genes:
                for type, strand, start, end, gene_type, codon in trnascan_genes[old_id]:
                    out.write(f"{start}\t{end}\ttRNA\n")
                    out.write(f"\t\t\tproduct\ttRNA-{gene_type}\n")
                    if type == 'pseudogene':
                        out.write(f"\t\t\tpseudogene\tunknown\n")
                # No match found
                # out.write(f"{new_id}\t0\t0\n")

            # Update record
            record.id = new_id
            record.name = new_id
            record.description = ""
            record.seq = modified_seq
            new_records.append(record)

            id_map.append((old_id, new_id))

    # Save the new transcript fasta with modified IDs and sequences
    with open(new_transcript_file, "w") as f_out:
        SeqIO.write(new_records, f_out, "fasta")

    print(f"  Saved modified transcripts to {new_transcript_file}")

    # Save ID mapping
    with open(id_map_file, "w") as map_out:
        for original_id, new_id in id_map:
            map_out.write(f"{original_id}\t{new_id}\n")

    print(f"  Saved ID map to {id_map_file}")


    # if Path(output_tbl).exists():
    #     print(f"  Skipping TBL generation (already exists)")
    #     return
    # protein2pos=parse_miniprot_results(paf_file)
    # with open(output_tbl, "w") as out:
    #     for record in SeqIO.parse(transcript_file, "fasta"):
    #         out.write(f">Feature {record.id}\n")
    #         seq_len = len(record.seq)
    #         out.write(f"1\t{seq_len}\tgene\n")
    #         if record.id in protein2pos:
    #             for start, end, strand in protein2pos[record.id]:
    #                 out.write(f"{start}\t{end}\tCDS\n")
    #         else:
    #             # If no match found, write a hypothetical sequence
    #             out.write(f"{record.id}\t0\t0\n")
        

    #     # TODO: replace with real annotation logic

def main(transcript_list_file, output_dir):

    # Validate paths for required executables
    check_exec_path("DIAMOND", DIAMOND_BIN)
    check_exec_path("tRNAscan-SE", TRNASCAN_BIN)
    check_exec_path("Miniprot", MINIPROT_BIN)

    os.makedirs(output_dir, exist_ok=True)

    with open(transcript_list_file) as f:
        for line in f:
            transcript_path = Path(line.strip())
            if not transcript_path.name.endswith(".fix.fasta.gz"):
                continue

            base_prefix = transcript_path.name.replace("_transcript_", "_protein_").replace(".fix.fasta.gz", ".fasta.gz")
            protein_path = transcript_path.parent / base_prefix
            genotype = extract_genotype_from_filename(transcript_path.name)

            if protein_path.suffix == ".gz":
                uncompressed_protein_path = Path(output_dir) / protein_path.with_suffix('').name
                if not uncompressed_protein_path.exists():
                    print(f"  Decompressing protein file: {protein_path.name}")
                    with gzip.open(protein_path, 'rt') as f_in, open(uncompressed_protein_path, 'w') as f_out:
                        f_out.writelines(f_in)
            else:
                uncompressed_protein_path = protein_path

            print(f"\nProcessing genotype: {genotype}")
            print(f"Transcript file: {transcript_path}")
            print(f"Protein file: {protein_path}")

            modified_transcript_filename = transcript_path.name.replace("fix", "fix2").removesuffix(".gz")
            modified_transcript_path = Path(output_dir) / modified_transcript_filename

            if modified_transcript_path.exists():
                print("  Skipping transcript header modification (already exists)")
            else:
                modify_transcript_headers(transcript_path, genotype, modified_transcript_path)

            base_output_name = protein_path.stem
            output_prefix = Path(output_dir) / base_output_name

            run_trnascan_chunked_parallel(
                str(modified_transcript_path),
                str(modified_transcript_path.with_suffix('')),
                chunk_size=5000,
                max_workers=12
            )

            paf_out = Path(output_dir) / f"{base_output_name}.miniprot.paf"
            run_miniprot(str(protein_path), str(modified_transcript_path), str(paf_out))

            run_diamond(
                str(uncompressed_protein_path),
                DIAMOND_DATABASES,
                str(output_prefix),
                threads=60
            )

            run_ahrd(
                str(uncompressed_protein_path),
                str(output_dir),
                base_output_name
            )

            dbcan_output_dir = Path(output_dir) / f"{base_output_name}.dbcan"
            run_dbcan(str(uncompressed_protein_path), str(dbcan_output_dir))
            


            if protein_path.suffix == ".gz" and uncompressed_protein_path.exists():
                print(f"  Removing temporary uncompressed protein file: {uncompressed_protein_path}")
                uncompressed_protein_path.unlink()

            output_tbl= Path(output_dir) / transcript_path.name.replace(".fix.fasta.gz", ".2ncbi.tbl")
            new_transcript_file = Path(output_dir) / transcript_path.name.replace(".fix.fasta.gz", ".2ncbi.fna")
            id_map_file = Path(output_dir) / transcript_path.name.replace(".fix.fasta.gz", ".2ncbi.id_map.txt")
            ahrd_out=Path(output_dir) / f"{base_output_name}.ahrd.csv"
            trnascan_file=f"{str(modified_transcript_path.with_suffix(''))}.trnascan.txt"
            #Processes the results
            parse_results_and_generate_tbl(genotype,str(modified_transcript_path), ahrd_out, paf_out, trnascan_file, output_tbl, new_transcript_file, id_map_file)
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Annotate protein sequences using DIAMOND.")
    parser.add_argument("--transcript_file_list", required=True, help="File of filename with a list of transcript .fix.fasta.gz files")
    parser.add_argument("--output_dir", required=True, help="Directory to store outputs")
    args = parser.parse_args()

    main(args.transcript_file_list, args.output_dir)
