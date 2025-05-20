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
DIAMOND_DATABASES = ["DBs/uniprot_sprot.fasta","DBs/uniprot_trembl.fasta"]
# ==== END GLOBAL SETTINGS ====

def extract_genotype_from_filename(filename):
    filename_prefix = 'sugarcanePanTranscriptome_06052025_'
    filename_suffix = '_transcript_4120596.fix.fasta.gz'
    part = filename.replace(filename_prefix, "").replace(filename_suffix, "")
    return part

def check_exec_path(name, path):
    if shutil.which(path) is None:
        raise FileNotFoundError(f"Executable for {name} not found at '{path}'. Check PATH or environment variable.")
    
def modify_transcript_headers(transcripts_fasta_gz, genotype, output_fasta):
    organism_str = f"[moltype=mRNA] [organism=Saccharum hybrid cultivar {genotype}]"
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
        print(f"  Skipping Miniprot (already exists)")
        return
    print(f"  Running Miniprot...")
    with open(output_file, "w") as out:
        subprocess.run([MINIPROT_BIN, "-t", "15", transcripts_fasta, proteins_fasta], stdout=out, check=True)

def parse_miniprot_results(miniprot_file):
    with open(miniprot_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if fields[0] == fields[5] and fields[1]*3 == fields[6]:
                continue

def parse_results_and_generate_tbl(blast_file, paf_file, output_tbl):
    if Path(output_tbl).exists():
        print(f"  Skipping TBL generation (already exists)")
        return
    with open(output_tbl, "w") as out:
        out.write(">Feature hypothetical_sequence\n")
        # TODO: replace with real annotation logic

def main(transcript_list_file, swissprot_db, output_dir):

    # Validate paths for required executables
    check_exec_path("DIAMOND", DIAMOND_BIN)
    check_exec_path("tRNAscan-SE", TRNASCAN_BIN)
    check_exec_path("Miniprot", MINIPROT_BIN)

    os.makedirs(output_dir, exist_ok=True)
    additional_db = "DBs/uniprot_trembl.fasta"

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

            if protein_path.suffix == ".gz" and uncompressed_protein_path.exists():
                print(f"  Removing temporary uncompressed protein file: {uncompressed_protein_path}")
                uncompressed_protein_path.unlink()

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Annotate protein sequences using DIAMOND.")
    parser.add_argument("--transcript_file_list", required=True, help="List of transcript .fix.fasta.gz files")
    parser.add_argument("--output_dir", required=True, help="Directory to store outputs")
    args = parser.parse_args()

    main(args.transcript_file_list, args.output_dir)
