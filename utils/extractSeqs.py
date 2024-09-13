import sys
import argparse
#import re
import os
import logging
#from Bio import SeqIO
from peewee import fn
import gzip

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from db.model import Sequence, panTranscriptomeGroup

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
def extract_sequences(sequence_class, representative, output_file):
    page_size = 100000
    page_number = 1
    has_more_results = True

    while has_more_results:
        logging.info(f'M1. Extracting {sequence_class} sequences')

        query = Sequence.select(Sequence.sequenceIdentifier, Sequence.sequence).where(Sequence.sequenceClass == sequence_class)

        # Apply representative filter if necessary
        if representative and sequence_class != 'protein':
            logging.info(f'M2. Filtering by representative sequences for {sequence_class}')
            subquery = (Sequence
                        .select(Sequence.sequenceIdentifier)
                        .join(panTranscriptomeGroup, on=(Sequence.ID == panTranscriptomeGroup.sequenceID))
                        .where(panTranscriptomeGroup.representative == True))
            query = query.where(Sequence.sequenceIdentifier.in_(subquery)).paginate(page_number, page_size)
        elif representative and sequence_class == 'protein':
            logging.info(f'M3. Filtering by representative sequences for {sequence_class}')
            query = query.join(panTranscriptomeGroup, on=(Sequence.ID == panTranscriptomeGroup.sequenceID)).where(panTranscriptomeGroup.representative == True).paginate(page_number, page_size)
        else:
            query = query.paginate(page_number, page_size)

        sequences = query.execute()

        if not sequences:
            has_more_results = False
        else:
            with (gzip.open(output_file, 'at') if args.gzip else open(output_file, 'a')) as fasta_file:
                for sequence in sequences:
                    fasta_file.write(f'>{sequence.sequenceIdentifier}\n{sequence.sequence}\n')
                    fasta_file.flush()
            page_number += 1

# Loop over sequence classes and extract the corresponding sequences
for sequence_class in sequence_classes:
    outfile = f"{args.prefix}_{sequence_class}_{current_pid}"
    if args.representatives:
        outfile += '_representatives'
    
    # Add .fasta or .fasta.gz extension based on gzip flag
    outfile += '.fasta.gz' if args.gzip else '.fasta'

    extract_sequences(sequence_class, args.representatives, outfile)

logging.info("Sequence extraction completed.")