#!/usr/bin/env python3
import sys
from collections import defaultdict

def fix_fasta_duplicates(input_file, output_file):
    header_count = defaultdict(int)
    headers_seen = set()
    duplicate_count = 0

    with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
        current_header = None
        for line in f_in:
            if line.startswith('>'):
                # 处理序列头
                original_header = line.strip()
                header_count[original_header] += 1

                if header_count[original_header] > 1:
                    # 重复的头，添加唯一后缀
                    new_header = f"{original_header}_dup{header_count[original_header]}"
                    f_out.write(new_header + '\n')
                    duplicate_count += 1
                    print(f"修复重复: {original_header} -> {new_header}")
                else:
                    # 第一次出现的头
                    f_out.write(line)

                current_header = original_header
            else:
                # 序列行
                f_out.write(line)

    print(f"总共修复了 {duplicate_count} 个重复序列标识符")
    return duplicate_count

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("用法: python fix_duplicates.py 输入文件 输出文件")
        sys.exit(1)

    input_fasta = sys.argv[1]
    output_fasta = sys.argv[2]
    fix_fasta_duplicates(input_fasta, output_fasta)
