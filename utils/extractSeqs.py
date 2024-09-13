# extractSeqs.py

import sys
import argparse
import re
import os
from Bio import SeqIO

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

current_pid = os.getpid()

from db.model import Sequence, Sequence2Set, SequenceSet, panTranscriptomeGroup

#initialize the parser 
parser = argparse.ArgumentParser(description='Extract sequences from the database')
parser.add_argument('--cds', action='store_true', help='Boolean - Extract CDD sequences')
parser.add_argument('--proteins', action='store_true', help='Boolean - Extract Protein sequences')
parser.add_argument('--transcripts', action='store_true', help='Boolean - Extract Transcript sequences')
parser.add_argument('--genes', action='store_true', help='Boolean - Extract Gene sequences')
parser.add_argument('--representatives', action='store_true', help='Boolean - Extract only Representative sequences, at least one fof --cds, --proteins, or --transcripts must be specified')
parser.add_argument('--prefix', type=str, required=True, help='Prefix to be used to create output files')
args = parser.parse_args()

if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)

if args.representatives and not (args.cds or args.proteins or args.transcripts):
    print("At least one of --cds, --proteins, or --transcripts must be specified when using --representatives")
    sys.exit(1)

if not (args.cds or args.proteins or args.transcripts):
    print("At least one of --cds, --proteins, or --transcripts must be specified")
    sys.exit(1)

sequenceClasses=[]

if args.cds:
    sequenceClasses.append('CDS')
if args.proteins:
    sequenceClasses.append('protein')
if args.transcripts:
    sequenceClasses.append('transcript')
if args.genes:
    sequenceClasses.append('gene')

for sequenceClass in sequenceClasses:
    page_size = 100000
    page_number = 1
    has_more_results = True
    while has_more_results:
        print(f'M1. Extracting {sequenceClass} sequences', file=sys.stdout)
        outfile=f'{args.prefix}_{sequenceClass}_{current_pid}'
        sequences = Sequence.select(Sequence.sequenceIdentifier,Sequence.sequence).where(Sequence.sequenceClass == sequenceClass).paginate(page_number, page_size)
        if args.representatives and args.proteins:
            print(f'M2. Extracting only representatives for {sequenceClass} sequences', file=sys.stdout)
            outfile=f'{outfile}_representatives'
            sequences = sequences.join(panTranscriptomeGroup, on=(Sequence.ID == panTranscriptomeGroup.sequenceID)).where(panTranscriptomeGroup.representative == True)
        elif args.representatives and (args.cds or args.transcripts or args.genes):
            print(f'M3. Extracting only representatives for {sequenceClass} sequences', file=sys.stdout)
            outfile=f'{outfile}_representatives'
            subquery = (Sequence
                    .select(Sequence.sequenceIdentifier)
                    .join(panTranscriptomeGroup, on=(Sequence.ID == panTranscriptomeGroup.sequenceID))
                    .where(panTranscriptomeGroup.representative == True))
            query = (Sequence
                .select(Sequence.sequenceIdentifier,Sequence.sequence)
                .where(
                    (Sequence.sequenceClass == sequenceClass) &
                    (Sequence.sequenceIdentifier.in_(subquery))
                )
                .paginate(page_number, page_size))
            sequences = query.execute()
            
        outfile=f'{outfile}.fasta'
        
        if not sequences:
            has_more_results = False
        else:
            with open(outfile, 'a') as f:
                for sequence in sequences:        
                    f.write(f'>{sequence.sequenceIdentifier}\n{sequence.sequence}\n')
            page_number += 1