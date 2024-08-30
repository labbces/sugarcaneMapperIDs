# populate.py

import sys
import argparse
import re
import os
from Bio import SeqIO

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from db.model import Sequence, Sequence2Set, SequenceSet, panTranscriptomeGroup


#initialize the parser 
parser = argparse.ArgumentParser(description='Populate the database with sequences and orthogroups')
parser.add_argument('--cds', type=argparse.FileType('r'), required=True, help='The file containing the coding sequences')
parser.add_argument('--proteins', type=argparse.FileType('r'), required=True, help='The file containing the protein sequences')
parser.add_argument('--transcripts', type=argparse.FileType('r'), required=True, help='The file containing the transcript sequences')
parser.add_argument('--orthogroup_members', type=argparse.FileType('r'), required=True, help='The file containing the orthogroup members')
parser.add_argument('--orthogroup_representative', type=argparse.FileType('r'), required=True, help='The file containing the sequence IDs of representative for each orthogroup')
args = parser.parse_args()

if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)

def check_existing_sequence(sequenceIdentifier, sequenceClass, sequenceVersion):
    return Sequence.get_or_none(
                        (Sequence.sequenceIdentifier == sequenceIdentifier) &
                        (Sequence.sequenceClass == sequenceClass) &
                        (Sequence.sequenceVersion == sequenceVersion)
                    )

def insert_new_sequence(sequenceIdentifier, sequenceLength, sequenceType, sequence, sequenceClass, sequenceVersion):
    newSeq=Sequence.create(
        sequenceIdentifier=sequenceIdentifier,
        sequenceLength=sequenceLength,
        sequenceType=sequenceType,
        sequence=sequence,
        sequenceClass=sequenceClass,
        sequenceVersion=sequenceVersion
    )
    return newSeq

def check_existing_set(sequenceSet):
    return SequenceSet.get_or_none(
                        (SequenceSet.nameSet == sequenceSet)
                    )

def insert_new_set(sequenceSet):
    return SequenceSet.create(
        nameSet=sequenceSet
    )

def get_set(sequenceIdentifier):
    # Regex patterns updated to capture the content inside parentheses
    pattern1 = re.compile(r'^(.+)_k[23][51].+$')
    pattern2 = re.compile(r'^(SP80-3280).+$')
    pattern3 = re.compile(r'^(R570).+$')
    pattern4 = re.compile(r'^(CC011940).+$')

    # Matching each pattern with the sequenceIdentifier
    match1 = pattern1.match(sequenceIdentifier)
    match2 = pattern2.match(sequenceIdentifier)
    match3 = pattern3.match(sequenceIdentifier)
    match4 = pattern4.match(sequenceIdentifier)

    # Checking which pattern matches and returning the result
    if match1:
        return match1.group(1)
    elif match2:
        return match2.group(1)
    elif match3:
        return match3.group(1)
    elif match4:
        return match4.group(1)
    else:
        return None

def get_set_id(sequenceIdentifier):
    sequenceSet=get_set(sequenceIdentifier)
    existing_set=check_existing_set(sequenceSet)
    if not existing_set:
        # print(f'Inserting new set: {sequenceSet}')
        set = insert_new_set(sequenceSet)
    else:
        set = existing_set
    return set.seID

def insert_new_Sequence2Set(sequenceID, sequenceSetID):
    Sequence2Set.create(
        sequenceID=sequenceID,
        seID=sequenceSetID
    )

def process_sequenceFile(fileObj, sequenceClass, sequenceVersion, sequenceType):
    sequence = ''
    sequenceIdentifier = None
    print(f'Processing {fileObj.name}... sequenceClass: {sequenceClass}, sequenceVersion: {sequenceVersion}, sequenceType: {sequenceType}')
    with open (fileObj.name) as seqHdl:
        for record in SeqIO.parse(seqHdl, "fasta"):
            sequenceIdentifier = record.id
            sequenceLength = len(record.seq)
            sequence = str(record.seq)

            setId=get_set_id(sequenceIdentifier)

            # Check if the sequence is already stored in the database
            existing_sequence = check_existing_sequence(sequenceIdentifier, sequenceClass, sequenceVersion)

            if not existing_sequence:
                # If not found, insert it into the database
                newSeq= insert_new_sequence(sequenceIdentifier, sequenceLength, sequenceType, sequence, sequenceClass, sequenceVersion)
                insert_new_Sequence2Set(newSeq.ID, setId)
        
def process_orthogroups_file(infile):
    representatives={}
    print(f'Processing {infile.name}...')
    for line in args.orthogroup_representative:
        line = line.strip()
        sequenceIdentifier = line.split()[0]
        representatives[sequenceIdentifier]=True

    for line in infile:
        line = line.strip()
        (orthogroup, sequenceIdentifier) = line.split('\t')
        if not orthogroup.startswith('OG'):
            print(f"Invalid orthogroup: {orthogroup}", file=sys.stderr)
        seq = Sequence.get_or_none(Sequence.sequenceIdentifier == sequenceIdentifier, Sequence.sequenceClass == 'protein', Sequence.sequenceVersion == 1)
        if seq:
            if sequenceIdentifier in representatives:
                panTranscriptomeGroup.get_or_create(sequenceID=seq.ID, groupID=orthogroup, representative=1)
            else:
                panTranscriptomeGroup.get_or_create(sequenceID=seq.ID, groupID=orthogroup, representative=0)
        # print(f'{orthogroup}: {sequenceIdentifier}')

#Reading the protein sequences
process_sequenceFile(args.proteins, 'protein', 1, 'Amino acid')
# Read orthogroup members
process_orthogroups_file(args.orthogroup_members)    
# # Read the coding sequences
process_sequenceFile(args.cds, 'CDS', 1, 'DNA')
# #Reading the transcript sequences
process_sequenceFile(args.transcripts, 'transcript', 1, 'DNA')