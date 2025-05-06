# exportSequenceFiles.py
import os
import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from db.model import Sequence, Sequence2Set, SequenceSet, panTranscriptomeGroup
import argparse

#initialize the parser 
parser = argparse.ArgumentParser(description='extract sequences from the database and write them to files')
parser.add_argument('--prefix', type=str, required=True, help='string to use as prefix of output filenames')
args = parser.parse_args()

if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)

sets = SequenceSet.select(SequenceSet.nameSet).distinct()
sequenceClassExt={'transcript':'fna', 'CDS':'fna', 'protein': 'faa'}

for set in sets:
    basename=f"{args.prefix}_{set.nameSet}"
    for ext in sequenceClassExt:
        with open(f"{basename}_{ext}.{sequenceClassExt[ext]}", 'w') as f, open(f"{basename}_representative_{ext}.{sequenceClassExt[ext]}", 'w') as f2:
            for seq in Sequence.select().join(Sequence2Set).join(SequenceSet).where(SequenceSet.nameSet==set.nameSet, Sequence.sequenceClass==ext, Sequence.sequenceVersion==1):
                print(seq.ID)
                og=panTranscriptomeGroup.select(panTranscriptomeGroup.groupID,panTranscriptomeGroup.representative).where(panTranscriptomeGroup.sequenceID==seq.ID).first()
                # print(og)
                # f.write(f">{seq.sequenceIdentifier} OG={og.groupID} representative={og.representative}\n{seq.sequence}\n")
                # if og.representative == 1:
                    # f2.write(f">{seq.sequenceIdentifier} OG={og.groupID} representative={og.representative}\n{seq.sequence}\n")