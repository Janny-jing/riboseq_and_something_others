def reverse_complement(seq):
    """
    生成 DNA 序列的反向互补序列。
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
                  'a': 't', 'c': 'g', 'g': 'c', 't': 'a'}
    return ''.join(complement[base] for base in reversed(seq))

def process_fasta(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        header = ''
        sequence = ''
        
        for line in infile:
            if line.startswith('>'):
                if header and sequence:
                    # 写入反向互补序列
                    outfile.write(f"{header}\n{reverse_complement(sequence)}\n")
                    sequence = ''
                header = line.strip()
            else:
                sequence += line.strip()

        # 处理最后一个序列
        if header and sequence:
            outfile.write(f"{header}\n{reverse_complement(sequence)}\n")

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 3:
        print("Usage: python generate_reverse_complement.py <input_fasta> <output_fasta>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    process_fasta(input_file, output_file)
