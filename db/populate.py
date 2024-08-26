# populate.py

from model import Sequence, Sequence2Set, SequenceSet, panTranscriptomeGroup
import argparse
import sys
import re

#initialize the parser 
parser = argparse.ArgumentParser(description='Populate the database with sequences and orthogroups')
parser.add_argument('--cds', type=argparse.FileType('r'), help='The file containing the coding sequences')
parser.add_argument('--proteins', type=argparse.FileType('r'), help='The file containing the protein sequences')
parser.add_argument('--orthogroups', type=argparse.FileType('r'), help='The file containing the orthogroups')
parser.add_argument('--orthogroup_members', type=argparse.FileType('r'), help='The file containing the orthogroup members')
parser.add_argument('--sequence_set', type=argparse.FileType('r'), help='The file containing the sequence sets')
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
    Sequence.create(
        sequenceIdentifier=sequenceIdentifier,
        sequenceLength=sequenceLength,
        sequenceType=sequenceType,
        sequence=sequence,
        sequenceClass=sequenceClass,
        sequenceVersion=sequenceVersion
    )

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
    pattern1 = re.compile(r'^>(.+)\(k[23][51]\).+$')
    pattern2 = re.compile(r'^>(SP80-3280)\(.+\)$')
    pattern3 = re.compile(r'^>(R570)\(.+\)$')
    pattern4 = re.compile(r'^>(CC011940)\(.+\)$')

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
# Read the coding sequences
if args.cds:
    sequence = ''
    sequenceIdentifier = None
    sequenceClass = 'CDS'
    sequenceVersion = 1
    sequenceType='DNA'
    with open (args.cds.name) as cds:
        for line in cds:
            line = line.strip()
            if line.startswith('>'):
                # If a sequence is already accumulated, insert it before processing the next one
                if sequenceIdentifier and sequence:
                    sequenceLength = len(sequence)

                    # Check if the sequence is already stored in the database
                    existing_sequence = check_existing_sequence(sequenceIdentifier, sequenceClass, sequenceVersion)

                    if not existing_sequence:
                        # If not found, insert it into the database
                        insert_new_sequence(sequenceIdentifier, sequenceLength, sequenceType, sequence, sequenceClass, sequenceVersion)

                # Reset for the new sequence
                sequenceIdentifier = line.split()[0][1:]
                print(sequenceIdentifier)
                sequenceSet=get_set(sequenceIdentifier)
                existing_set=check_existing_set(sequenceSet)
                if not existing_set:
                    insert_new_set(sequenceSet)
                sequence = ''
            else:
                sequence += line

        # Don't forget to insert the last sequence after exiting the loop
        if sequenceIdentifier and sequence:
            sequenceLength = len(sequence)

            # Check if the sequence is already stored in the database
            existing_sequence = check_existing_sequence(sequenceIdentifier, sequenceClass, sequenceVersion)

            if not existing_sequence:
                # If not found, insert it into the database
                insert_new_sequence(sequenceIdentifier, sequenceLength, sequenceType, sequence, sequenceClass, sequenceVersion)