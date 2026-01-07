import os
import pandas as pd
import argparse

def read_and_aggregate(file_path, output_path, rename_dict):
    # 读取文件
    df = pd.read_csv(file_path, sep='\t', header=None, names=['code', 'value'])

    # 设置 'code' 作为行索引
    df.set_index('code', inplace=True)

    # 使用 .rename() 方法来更新行索引
    df = df.rename(index=rename_dict)

    # 聚合具有相同行名的值
    df = df.groupby(level=0).sum()

    # 将修改后的 DataFrame 写入新文件
    df.reset_index().to_csv(output_path, sep='\t', header=False, index=False)

def process_files(directory_path, output_directory, rename_dict):
    # 获取目录中的所有文件
    files = [f for f in os.listdir(directory_path) if os.path.isfile(os.path.join(directory_path, f))]

    for file_name in files:
        input_file_path = os.path.join(directory_path, file_name)
        output_file_path = os.path.join(output_directory, f"aggregated_{file_name}")

        read_and_aggregate(input_file_path, output_file_path, rename_dict)

def main():
    parser = argparse.ArgumentParser(description="Batch process text files by renaming row names and aggregating values.")
    parser.add_argument("input_directory", help="Path to the directory containing the input text files.")
    parser.add_argument("output_directory", help="Path to the directory where the output files will be saved.")

    args = parser.parse_args()

    # 创建输出目录（如果不存在）
    os.makedirs(args.output_directory, exist_ok=True)

    # 创建一个字典，其中键是旧行名，值是新行名
    rename_dict = {
        'Ala_GCC': 'Ala-AGC',
        'Ala_GCU': 'Ala-AGC',
        'Ala_GCG': 'Ala-CGC',
        'Ala_GCA': 'Ala-TGC',
        'Arg_CGU': 'Arg-ACG',
        'Arg_CGC': 'Arg-ACG',
        'Arg_CGG': 'Arg-CCG',
        'Arg_AGG': 'Arg-CCT',
        'Arg_CGA': 'Arg-TCG',
        'Arg_AGA': 'Arg-TCT',
        'Asn_AAC': 'Asn-GTT',
        'Asn_AAU': 'Asn-GTT',
        'Asp_GAC': 'Asp-GTC',
        'Asp_GAU': 'Asp-GTC',
        'Cys_UGC': 'Cys-GCA',
        'Cys_UGU': 'Cys-GCA',
        'Gln_CAG': 'Gln-CTG',
        'Gln_CAA': 'Gln-TTG',
        'Glu_GAG': 'Glu-CTC',
        'Glu_GAA': 'Glu-TTC',
        'Gly_GGG': 'Gly-CCC',
        'Gly_GGC': 'Gly-GCC',
        'Gly_GGU': 'Gly-GCC',
        'Gly_GGA': 'Gly-TCC',
        'His_CAC': 'His-GTG',
        'His_CAU': 'His-GTG',
        'Ile_AUU': 'Ile-AAT',
        'Ile_AUC': 'Ile-GAT',
        'Ile_AUA': 'Ile-TAT',
        'Leu_CUC': 'Leu-AAG',
        'Leu_CUU': 'Leu-AAG',
        'Leu_UUG': 'Leu-CAA',
        'Leu_CUG': 'Leu-CAG',
        'Leu_UUA': 'Leu-TAA',
        'Leu_CUA': 'Leu-TAG',
        'Lys_AAG': 'Lys-CTT',
        'Lys_AAA': 'Lys-TTT',
        'Met_AUG': 'Met-CAT',
        'Phe_UUU': 'Phe-GAA',
        'Phe_UUC': 'Phe-GAA',
        'Pro_CCC': 'Pro-AGG',
        'Pro_CCU': 'Pro-AGG',
        'Pro_CCG': 'Pro-CGG',
        'Pro_CCA': 'Pro-TGG',
        'TER_UAA': 'SeC-TCA',
        'TER_UAG': 'SeC-TCA',
        'TER_UGA': 'SeC-TCA',
        'Ser_UCC': 'Ser-AGA',
        'Ser_UCU': 'Ser-AGA',
        'Ser_UCG': 'Ser-CGA',
        'Ser_AGC': 'Ser-GCT',
        'Ser_AGU': 'Ser-GCT',
        'Ser_UCA': 'Ser-TGA',
        'Thr_ACC': 'Thr-AGT',
        'Thr_ACU': 'Thr-AGT',
        'Thr_ACG': 'Thr-CGT',
        'Thr_ACA': 'Thr-TGT',
        'Trp_UGG': 'Trp-CCA',
        'Tyr_UAU': 'Tyr-ATA',
        'Tyr_UAC': 'Tyr-GTA',
        'Val_GUC': 'Val-AAC',
        'Val_GUU': 'Val-AAC',
        'Val_GUG': 'Val-CAC',
        'Val_GUA': 'Val-TAC'
        # 添加更多需要更改的行名
    }

    # 调用函数
    process_files(args.input_directory, args.output_directory, rename_dict)

if __name__ == "__main__":
    main()
