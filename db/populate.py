# populate.py

from model import Sequence, Sequence2Set, SequenceSet, panTranscriptomeGroup
import argparse
import sys
import re

#initialize the parser 
parser = argparse.ArgumentParser(description='Populate the database with sequences and orthogroups')
parser.add_argument('--cds', type=argparse.FileType('r'), required=True, help='The file containing the coding sequences')
parser.add_argument('--proteins', type=argparse.FileType('r'), required=True, help='The file containing the protein sequences')
parser.add_argument('--transcripts', type=argparse.FileType('r'), required=True, help='The file containing the transcript sequences')
parser.add_argument('--orthogroup_members', type=argparse.FileType('r'), help='The file containing the orthogroup members')
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
    SequenceSet.create(
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
    # print(sequenceSet)
    existing_set=check_existing_set(sequenceSet)
    if not existing_set:
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
    with open (fileObj.name) as cds:
        for line in cds:
            line = line.strip()
            if line.startswith('>'):
                # Reset for the new sequence
                sequenceIdentifier = line.split()[0][1:]
                # If a sequence is already accumulated, insert it before processing the next one
                if sequenceIdentifier and sequence:
                    sequenceLength = len(sequence)
                    # print(sequenceIdentifier)
                    setId=get_set_id(sequenceIdentifier)
                    # print(f'setid:{setId}')

                    # Check if the sequence is already stored in the database
                    existing_sequence = check_existing_sequence(sequenceIdentifier, sequenceClass, sequenceVersion)

                    if not existing_sequence:
                        # If not found, insert it into the database
                        newSeq= insert_new_sequence(sequenceIdentifier, sequenceLength, sequenceType, sequence, sequenceClass, sequenceVersion)
                        # print(f'New seq: {newSeq.ID}; Identifier={sequenceIdentifier}; setID={setId}')
                        insert_new_Sequence2Set(newSeq.ID, setId)
         
                #initialize sequence as blank string when a new identifier is found
                sequence = ''
            else:
                sequence += line

        # Don't forget to insert the last sequence after exiting the loop
        if sequenceIdentifier and sequence:
            sequenceLength = len(sequence)
            setId=get_set_id(sequenceIdentifier)

            # Check if the sequence is already stored in the database
            existing_sequence = check_existing_sequence(sequenceIdentifier, sequenceClass, sequenceVersion)

            if not existing_sequence:
                # If not found, insert it into the database
                newSeq= insert_new_sequence(sequenceIdentifier, sequenceLength, sequenceType, sequence, sequenceClass, sequenceVersion)
                insert_new_Sequence2Set(newSeq.ID, setId)
                
# Read the coding sequences
process_sequenceFile(args.cds, 'ÇDS', 1, 'DNA')
#Reading the transcript sequences
process_sequenceFile(args.transcripts, 'transcript', 1, 'DNA')
#Reading the protein sequences
process_sequenceFile(args.proteins, 'protein', 1, 'Amino acid')