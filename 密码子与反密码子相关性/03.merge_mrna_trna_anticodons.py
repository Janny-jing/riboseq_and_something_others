import os
import pandas as pd
import argparse

def merge_and_save(input_file, trna_file, output_file):
    # 读取 trna 文件
    trna_data = pd.read_csv(trna_file, sep='\t', header=None, names=['trnaid', 'trna_codons_proportion'])

    # 读取 mrna 文件
    mrna_data = pd.read_csv(input_file, sep='\t', header=None, names=['trnaid', 'mrna_codons_proportion'])

    # 合并数据
    merged_data = pd.merge(trna_data, mrna_data, on='trnaid', how='inner')

    # 保存合并后的数据
    merged_data.to_csv(output_file, sep='\t', index=False)

def process_files(input_directory, trna_file, output_directory):
    # 获取目录中的所有文件
    files = [f for f in os.listdir(input_directory) if os.path.isfile(os.path.join(input_directory, f))]

    # 创建输出目录（如果不存在）
    os.makedirs(output_directory, exist_ok=True)

    for file_name in files:
        input_file_path = os.path.join(input_directory, file_name)
        output_file_path = os.path.join(output_directory, file_name)

        merge_and_save(input_file_path, trna_file, output_file_path)

def main():
    parser = argparse.ArgumentParser(description="Merge two sets of files based on the first column.")
    parser.add_argument("input_directory", help="Path to the directory containing the mrna_codons_proportion files.")
    parser.add_argument("trna_file", help="Path to the trna_codons_proportion file.")
    parser.add_argument("output_directory", help="Path to the directory where the merged files will be saved.")

    args = parser.parse_args()

    # 调用函数
    process_files(args.input_directory, args.trna_file, args.output_directory)

if __name__ == "__main__":
    main()
