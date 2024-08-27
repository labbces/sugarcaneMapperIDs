# exportSequenceFiles.py

from model import Sequence, Sequence2Set, SequenceSet, panTranscriptomeGroup
import argparse
import sys

#initialize the parser 
parser = argparse.ArgumentParser(description='Populate the database with sequences and orthogroups')
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
                og=panTranscriptomeGroup.select(panTranscriptomeGroup.groupID,panTranscriptomeGroup.representative).where(panTranscriptomeGroup.sequenceID==seq.ID).first()
                f.write(f">{seq.sequenceIdentifier} OG={og.groupID} representative={og.representative}\n{seq.sequence}\n")
                if og.representative == 1:
                    f2.write(f">{seq.sequenceIdentifier} OG={og.groupID} representative={og.representative}\n{seq.sequence}\n")