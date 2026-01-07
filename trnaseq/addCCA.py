import sys

def read_fasta(file_path):
    sequences = []
    with open(file_path, 'r') as file:
        name = ''
        seq = ''
        for line in file:
            if line.startswith('>'):  # Header line
                if seq:
                    sequences.append((name, seq))
                    seq = ''
                name = line.strip()
            else:
                seq += line.strip().upper()  # Concatenate the sequence
        if seq:  # Append the last sequence
            sequences.append((name, seq))
    return sequences

def append_CCA_to_sequences(sequences):
    updated_sequences = []
    for name, seq in sequences:
        updated_seq = seq + 'CCA'  # Append CCA to the end of the sequence
        updated_sequences.append((name, updated_seq))
    return updated_sequences

def write_fasta(sequences, file_path):
    with open(file_path, 'w') as file:
        for name, seq in sequences:
            file.write(f"{name}\n{seq}\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python append_CCA.py <input_file> <output_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    # 读取原始序列
    original_sequences = read_fasta(input_file)

    # 给每个序列追加 CCA
    updated_sequences = append_CCA_to_sequences(original_sequences)

    # 将更新后的序列写入新的 FASTA 文件
    write_fasta(updated_sequences, output_file)
