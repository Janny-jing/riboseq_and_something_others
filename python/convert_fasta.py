#!/usr/bin/env python3
import argparse

def convert_to_two_columns(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # 确保文件中的行数是偶数
    if len(lines) % 2 != 0:
        print("警告：文件中的行数不是偶数，最后一行可能无法配对。")

    for i in range(0, len(lines), 2):
        id_line = lines[i].strip()[1:]  # 去掉 '>' 和行尾的换行符
        seq_line = lines[i + 1].strip() if i + 1 < len(lines) else ''
        print(f"{id_line}\t{seq_line}")

def main():
    parser = argparse.ArgumentParser(description='Convert FASTA format to two-column format.')
    parser.add_argument('file_path', type=str, help='Path to the input FASTA file')
    args = parser.parse_args()

    convert_to_two_columns(args.file_path)

if __name__ == "__main__":
    main()
