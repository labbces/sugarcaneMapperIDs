import os
import subprocess
import gzip
from pathlib import Path
from Bio import SeqIO

def extract_genotype_from_filename(filename):
    """Extract genotype name (e.g., US851008) from transcript filename."""
    filename_prefix = 'sugarcanePanTranscriptome_06052025_'
    filename_suffix = '_transcript_4120596.fix.fasta.gz'
    part = filename.replace(filename_prefix, "").replace(filename_suffix, "")
    return part

def modify_transcript_headers(transcripts_fasta_gz, genotype, output_fasta):
    """Add [moltype=mRNA] and [organism=...] to FASTA headers."""
    organism_str = f"[moltype=mRNA] [organism=Saccharum hybrid cultivar {genotype}]"
    with gzip.open(transcripts_fasta_gz, "rt") as input_handle, open(output_fasta, "w") as output_handle:
        for record in SeqIO.parse(input_handle, "fasta"):
            record.description = f"{record.id} {organism_str}"
            SeqIO.write(record, output_handle, "fasta")
    return output_fasta

def ensure_blast_db(fasta_path):
    """Create BLAST database if it does not exist (supports both split and non-split DBs)."""
    db_path = Path(fasta_path)
    parent = db_path.parent
    basename = db_path.name

    # Match both split (e.g., .00.pin) and non-split (.pin) databases
    pin_files = list(parent.glob(f"{basename}*.pin"))
    psq_files = list(parent.glob(f"{basename}*.psq"))
    phr_files = list(parent.glob(f"{basename}*.phr"))

    if pin_files and psq_files and phr_files:
        print(f"  BLAST database for {fasta_path} already exists.")
        return

    print(f"  Creating BLAST database for {fasta_path}...")
    subprocess.run([
        "makeblastdb",
        "-in", str(fasta_path),
        "-dbtype", "prot"
    ], check=True)

def run_blastp_multi(query_fasta, db_paths, output_prefix):
    """Run BLASTP against multiple databases."""
    print(f'{query_fasta} -> {output_prefix}')
    for db in db_paths:
        ensure_blast_db(db)
        db_name = Path(db).stem
        output_file = f"{output_prefix}.blast.{db_name}.txt"
        if Path(output_file).exists():
            print(f"  Skipping BLASTP against {db_name} (already exists)")
            continue
        print(f"  Running BLASTP against {db_name}...")
        subprocess.run([
            "blastp",
            "-query", query_fasta,
            "-db", db,
            "-out", output_file,
            "-outfmt", "6",
            "-evalue", "1e-5",
            "-num_threads", "10"
        ], check=True)

def run_trnascan(input_fasta, output_file, stats_file):
    """
    Run tRNAscan-SE 2.0 with output and statistics file paths.
    """
    if Path(output_file).exists() and Path(stats_file).exists():
        print(f"  Skipping tRNAscan-SE (outputs already exist)")
        return

    print(f"  Running tRNAscan-SE on {input_fasta}")
    cmd = [
        "/usr/local/tRNAscan-SE-2.0.12/bin/tRNAscan-SE",
        "-E",
        "--thread", "10",
        "-o", output_file,
        "-m", stats_file,
        input_fasta
    ]

    subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)

def run_miniprot(proteins_fasta, transcripts_fasta, output_file):
    if Path(output_file).exists():
        print(f"  Skipping Miniprot (already exists)")
        return
    print(f"  Running Miniprot...")
    with open(output_file, "w") as out:
        subprocess.run([
            "miniprot", "-t", "15",
            transcripts_fasta, proteins_fasta
        ], stdout=out, check=True)

def parse_miniprot_results(miniprot_file):
    """Parse Miniprot results and return a list of tuples."""
    with open(miniprot_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if fields[0] == fields[5] and fields[1]*3 == fields[6]:#Revisar
                continue

def parse_results_and_generate_tbl(blast_file, paf_file, output_tbl):
    if Path(output_tbl).exists():
        print(f"  Skipping TBL generation (already exists)")
        return
    with open(output_tbl, "w") as out:
        out.write(">Feature hypothetical_sequence\n")
        # TODO: replace with real annotation logic

def main(transcript_list_file, swissprot_db, output_dir):
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
                uncompressed_protein_path = Path(output_dir) / protein_path.with_suffix('').name  # remove .gz
                if not uncompressed_protein_path.exists():
                    print(f"  Decompressing protein file: {protein_path.name}")
                    with gzip.open(protein_path, 'rt') as f_in, open(uncompressed_protein_path, 'w') as f_out:
                        f_out.writelines(f_in)
            else:
                uncompressed_protein_path = protein_path
            print(f"\nProcessing genotype: {genotype}")
            print(f"Transcript file: {transcript_path}")
            print(f"Protein file: {protein_path}")

            # Generate modified transcript filename with "fix2"
            modified_transcript_filename = transcript_path.name.replace("fix", "fix2").removesuffix(".gz")
            modified_transcript_path = Path(output_dir) / modified_transcript_filename

            if modified_transcript_path.exists():
                print("  Skipping transcript header modification (already exists)")
            else:
                modify_transcript_headers(transcript_path, genotype, modified_transcript_path)

            base_output_name = protein_path.stem
            output_prefix = Path(output_dir) / base_output_name
            paf_out = Path(output_dir) / f"{base_output_name}.miniprot.paf"
            tbl_out = Path(output_dir) / f"{base_output_name}.tbl"
            trnascan_out = Path(output_dir) / f"{base_output_name}.trnascan.txt"
            trnascan_stats = Path(output_dir) / f"{base_output_name}.trnascan.stats.txt"

            run_trnascan(str(modified_transcript_path), str(trnascan_out), str(trnascan_stats))


            run_blastp_multi(str(uncompressed_protein_path), [swissprot_db, additional_db], str(output_prefix))

            if protein_path.suffix == ".gz" and uncompressed_protein_path.exists():
                print(f"  Removing temporary uncompressed protein file: {uncompressed_protein_path}")
                uncompressed_protein_path.unlink()
            run_miniprot(str(protein_path), str(modified_transcript_path), str(paf_out))


            parse_results_and_generate_tbl(str(output_prefix) + f".blast.{Path(swissprot_db).stem}.txt", str(paf_out), str(tbl_out))

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Annotate protein sequences using BLAST and Miniprot.")
    parser.add_argument("--transcript_file_list", required=True, help="List of transcript .fix.fasta.gz files")
    parser.add_argument("--swissprot_db", required=True, help="Path to uniprot_sprot.fasta (will auto-index if needed)")
    parser.add_argument("--output_dir", required=True, help="Directory to store outputs")
    args = parser.parse_args()

    main(args.transcript_file_list, args.swissprot_db, args.output_dir)
