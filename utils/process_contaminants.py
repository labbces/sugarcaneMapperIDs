import argparse
import gzip
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq

def open_maybe_gz(filename, mode='rt'):
    """Open text file or gzipped text file transparently."""
    return gzip.open(filename, mode) if filename.endswith('.gz') else open(filename, mode)

def parse_region(region_str):
    """Converts '1..20' into (1, 20)"""
    start, end = region_str.strip().split('..')
    return int(start), int(end)

def load_contaminations(contaminant_file, offset):
    """
    Loads contaminant regions and determines action (REMOVE_START, REMOVE_END, MASK).
    Returns dict: {sequence_id: [(start, end, action), ...]}
    """
    contamination_dict = defaultdict(list)
    with open(contaminant_file) as f:
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            cols = line.strip().split('\t')
            if len(cols) != 5:
                continue
            if cols[2] != 'ACTION_TRIM':
                continue
            # print(line)
            seq_id, length_str, _, region_str, _ = cols
            length = int(length_str)
            # start, end = parse_region(region_str)

            # if start <= offset:
            #     action = 'REMOVE_START'
            # elif end >= (length - offset + 1):
            #     action = 'REMOVE_END'
            # else:
            #     action = 'MASK'
            regions = region_str.split(',')
            for subregion in regions:
                try:
                    start, end = parse_region(subregion)
                except ValueError:
                    print(f"Skipping malformed region '{subregion}' in line:\n{line}")
                    continue

                if start <= offset:
                    action = 'REMOVE_START'
                elif end >= (length - offset + 1):
                    action = 'REMOVE_END'
                else:
                    action = 'MASK'

                contamination_dict[seq_id].append((start, end, action))

            contamination_dict[seq_id].append((start, end, action))
    return contamination_dict

def process_sequence(seq_record, contaminations):
    """
    Applies trimming or masking to a sequence based on contaminations.
    Returns modified SeqRecord and a summary of actions.
    """
    original_length = len(seq_record.seq)
    sequence = str(seq_record.seq)
    action_types = []

    # Apply from end to start to avoid index shifting
    for start, end, action in sorted(contaminations, key=lambda x: x[0], reverse=True):
        start_idx = start - 1
        end_idx = end
        if action.startswith('REMOVE'):
            sequence = sequence[:start_idx] + sequence[end_idx:]
        elif action == 'MASK':
            sequence = sequence[:start_idx] + 'N' * (end_idx - start_idx) + sequence[end_idx:]
        action_types.append(action)

    seq_record.seq = Seq(sequence)
    summary = {
        'SequenceID': seq_record.id,
        'OriginalLength': original_length,
        'FinalLength': len(sequence),
        'TotalContaminations': len(contaminations),
        'Masked': action_types.count('MASK'),
        'Removed': action_types.count('REMOVE_START') + action_types.count('REMOVE_END'),
        'Actions': ",".join(sorted(set(action_types)))
    }
    return seq_record, summary

def main():
    parser = argparse.ArgumentParser(description="Process sequences to remove/mask contaminants and summarize changes.")
    parser.add_argument('--contaminants', required=True, help="TSV file with contaminant regions")
    parser.add_argument('--fasta', required=True, help="Input FASTA file (can be .gz)")
    parser.add_argument('--output', required=True, help="Output FASTA file (.gz supported)")
    parser.add_argument('--summary', required=True, help="Summary output TSV file")
    parser.add_argument('--offset', type=int, default=20, help="Offset in bp for start/end trimming")
    args = parser.parse_args()

    print(f"Loading contaminants from: {args.contaminants}")
    contaminations = load_contaminations(args.contaminants, args.offset)

    summary_data = []

    print(f"Processing sequences from: {args.fasta}")
    with open_maybe_gz(args.fasta, 'rt') as in_fasta, gzip.open(args.output, 'wt') as out_fasta:
        for record in SeqIO.parse(in_fasta, 'fasta'):
            if record.id in contaminations:
                record, summary = process_sequence(record, contaminations[record.id])
                summary_data.append(summary)
            else:
                summary_data.append({
                    'SequenceID': record.id,
                    'OriginalLength': len(record.seq),
                    'FinalLength': len(record.seq),
                    'TotalContaminations': 0,
                    'Masked': 0,
                    'Removed': 0,
                    'Actions': 'NONE'
                })
            SeqIO.write(record, out_fasta, 'fasta')

    print(f"Saving summary to: {args.summary}")
    with open(args.summary, 'w') as f:
        f.write("SequenceID\tOriginalLength\tFinalLength\tTotalContaminations\tMasked\tRemoved\tActions\n")
        for s in summary_data:
            f.write("{SequenceID}\t{OriginalLength}\t{FinalLength}\t{TotalContaminations}\t{Masked}\t{Removed}\t{Actions}\n".format(**s))

    print("All done.")

if __name__ == '__main__':
    main()
