# analyze_flanking_cli.py

from collections import Counter
import sys
import argparse

def read_fastq(filename):
    """Generator to read FASTQ file and yield (header, sequence, qual)"""
    open_func = open
    if filename.endswith('.gz'):
        import gzip
        open_func = lambda f: gzip.open(f, 'rt')
    
    with open_func(filename) as f:
        while True:
            try:
                header = next(f).strip()
                seq = next(f).strip()
                next(f)  # skip + line
                qual = next(f).strip()
                yield header, seq, qual
            except StopIteration:
                break

def find_flanking_sequences(fastq_file, target_seq, flank_len=20):
    """
    Search for target_seq in reads, extract flank_len bases before and after.
    Returns: (upstream_counter, downstream_counter)
    """
    upstream_list = []
    downstream_list = []
    target_len = len(target_seq)

    print(f"🔍 Searching for '{target_seq}'...", end='', flush=True)

    for header, seq, qual in read_fastq(fastq_file):
        start = 0
        while True:
            pos = seq.find(target_seq, start)
            if pos == -1:
                break
            # Extract upstream
            if pos >= flank_len:
                upstream = seq[pos - flank_len:pos]
                upstream_list.append(upstream)
            # Extract downstream
            end_pos = pos + target_len
            if len(seq) >= end_pos + flank_len:
                downstream = seq[end_pos:end_pos + flank_len]
                downstream_list.append(downstream)
            start = pos + 1  # Find overlapping or multiple occurrences

    upstream_counter = Counter(upstream_list)
    downstream_counter = Counter(downstream_list)
    print(f" found {sum(upstream_counter.values()) + sum(downstream_counter.values())} flanking segments.")

    return upstream_counter, downstream_counter

def main():
    parser = argparse.ArgumentParser(
        description="Find common sequences flanking one or more target motifs in a FASTQ file."
    )
    parser.add_argument('fastq_file', help='Path to the FASTQ file (can be .gz compressed)')
    parser.add_argument('targets', nargs='+', help='One or more target sequences to search for')
    parser.add_argument('--flank', type=int, default=20, help='Number of bases to extract before/after (default: 20)')
    parser.add_argument('--top', type=int, default=10, help='Show top N most common flanking sequences (default: 10)')

    args = parser.parse_args()

    flank_len = args.flank

    for target in args.targets:
        print(f"\n" + "="*60)
        print(f"🎯 TARGET: {target}")
        print("="*60)

        try:
            up_counter, down_counter = find_flanking_sequences(args.fastq_file, target, flank_len)

            total_up = sum(up_counter.values())
            total_down = sum(down_counter.values())

            print(f"\n📌 Upstream ({flank_len} bp before): {total_up} occurrences")
            for seq, cnt in up_counter.most_common(args.top):
                print(f"  {seq} : {cnt}")

            print(f"\n📌 Downstream ({flank_len} bp after): {total_down} occurrences")
            for seq, cnt in down_counter.most_common(args.top):
                print(f"  {seq} : {cnt}")

        except Exception as e:
            print(f"❌ Error processing target '{target}': {e}")

if __name__ == '__main__':
    main()
