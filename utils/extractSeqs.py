import sys
import argparse
#import re
import os
import logging
#from Bio import SeqIO
from peewee import fn
import gzip

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from db.model import Sequence, panTranscriptomeGroup, SequenceSet, Sequence2Set

# Get the current process ID
current_pid = os.getpid()

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Initialize the parser
parser = argparse.ArgumentParser(description='Extract sequences from the database')
parser.add_argument('--cds', action='store_true', help='Boolean - Extract CDS sequences')
parser.add_argument('--proteins', action='store_true', help='Boolean - Extract Protein sequences')
parser.add_argument('--transcripts', action='store_true', help='Boolean - Extract Transcript sequences')
parser.add_argument('--genes', action='store_true', help='Boolean - Extract Gene sequences')
parser.add_argument('--representatives', action='store_true', help='Boolean - Extract only Representative sequences, at least one of --cds, --proteins, or --transcripts must be specified')
parser.add_argument('--prefix', type=str, required=True, help='Prefix to be used to create output files')
parser.add_argument('--genotypes', action='store_true', help='Separate sequence by genotype. Boolean')
parser.add_argument('--gzip', action='store_true', help='Save output in gzip-compressed FASTA format')
args = parser.parse_args()

# Argument validation
if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)

if args.representatives and not (args.cds or args.proteins or args.transcripts):
    logging.error("At least one of --cds, --proteins, or --transcripts must be specified when using --representatives")
    sys.exit(1)

if not (args.cds or args.proteins or args.transcripts or args.genes):
    logging.error("At least one of --cds, --proteins, or --transcripts, or --genes must be specified")
    sys.exit(1)

# Create a list of sequence classes to process
sequence_classes = []
if args.cds:
    sequence_classes.append('CDS')
if args.proteins:
    sequence_classes.append('protein')
if args.transcripts:
    sequence_classes.append('transcript')
if args.genes:
    sequence_classes.append('gene')

# Function to extract sequences and write to a FASTA file
def extract_sequences(sequence_class, representative, output_file, genotype=None):
    page_size = 100000
    page_number = 1
    has_more_results = True
    logging.info(f"Extracting {sequence_class} sequences" + (f" for genotype: {genotype.nameSet}" if genotype else ""))
    while has_more_results:
        logging.info(f"Processing page {page_number} for {sequence_class} sequences"+ (f" for genotype: {genotype.nameSet}" if genotype else ""))

        query = Sequence.select(Sequence.sequenceIdentifier, Sequence.sequence)

        if genotype:
            query = query.join(Sequence2Set, on=(Sequence.ID == Sequence2Set.sequenceID)).where(
                (Sequence.sequenceClass == sequence_class) &
                (Sequence2Set.seID == genotype.seID)
            )
        else:
            query = query.where(Sequence.sequenceClass == sequence_class)

        if representative and sequence_class != 'protein':
            logging.info(f"Filtering by representative sequences for {sequence_class}")
            subquery = (Sequence
                        .select(Sequence.ID)
                        .join(panTranscriptomeGroup, on=(Sequence.ID == panTranscriptomeGroup.sequenceID))
                        .where(panTranscriptomeGroup.representative == True))
            query = query.where(Sequence.ID.in_(subquery))
        elif representative and sequence_class == 'protein':
            logging.info(f"Filtering by representative protein sequences")
            query = query.join(panTranscriptomeGroup, JOIN.LEFT_OUTER,
                               on=(Sequence.ID == panTranscriptomeGroup.sequenceID)) \
                         .where(panTranscriptomeGroup.representative == True)

        sequences = query.paginate(page_number, page_size).execute()

        if not sequences:
            has_more_results = False
        else:
            with (gzip.open(output_file, 'at') if args.gzip else open(output_file, 'a')) as fasta_file:
                for sequence in sequences:
                    fasta_file.write(f'>{sequence.sequenceIdentifier}\n{sequence.sequence}\n')
                    fasta_file.flush()
            page_number += 1

if args.genotypes:
    genotypes = SequenceSet.select()
    for genotype in genotypes:
        for sequence_class in sequence_classes:
            outfile = f"{args.prefix}_{genotype.nameSet}_{sequence_class}_{current_pid}"
            if args.representatives:
                outfile += '_representatives'
            outfile += '.fasta.gz' if args.gzip else '.fasta'
            extract_sequences(sequence_class, args.representatives, outfile, genotype=genotype)

else:
# Loop over sequence classes and extract the corresponding sequences
    for sequence_class in sequence_classes:
        outfile = f"{args.prefix}_{sequence_class}_{current_pid}"
        if args.representatives:
            outfile += '_representatives'
        
        # Add .fasta or .fasta.gz extension based on gzip flag
        outfile += '.fasta.gz' if args.gzip else '.fasta'

        extract_sequences(sequence_class, args.representatives, outfile)

logging.info("Sequence extraction completed.")
