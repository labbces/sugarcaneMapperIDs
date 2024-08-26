# checkIDs.py
import argparse
import sys
import re

#initialize the parser 
parser = argparse.ArgumentParser(description='Populate the database with sequences and orthogroups')
parser.add_argument('--file', type=argparse.FileType('r'), help='The file containing the sequences')
args = parser.parse_args()

if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)

pattern1 = re.compile(r'^>.+k[23][51].+$')
pattern2 = re.compile(r'^>SP80-3280_.+$')
pattern3 = re.compile(r'^>R570_.+$')
pattern4 = re.compile(r'^>CC011940_.+$')
# Read the coding sequences
if args.file:
    with open(args.file.name) as f:
        for line in f:
            line=line.rstrip()
            if line.startswith('>'):
                line=line.split(' ')
                id=line[0]
                if not pattern1.match(id) and not pattern2.match(id) and not pattern3.match(id) and not pattern4.match(id):
                    print(f"Invalid ID: {id}")