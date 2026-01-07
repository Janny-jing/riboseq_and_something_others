# match_tss_with_stop.py

import pandas as pd
import argparse
import os

def main(tss_file, stop_file, output_file=None):
    # 自动推断输出文件名
    if not output_file:
        base_name = os.path.splitext(os.path.basename(stop_file))[0]
        output_file = f"matched_{base_name}.csv"

    # 读取 TSS 表格
    tss_df = pd.read_csv(tss_file, sep=r'\s+', engine='python')

    # 构建 gene_id -> info 映射
    tss_map = {}
    for _, row in tss_df.iterrows():
        gene_id = row['gene_id']
        chromosome = str(row['chromosome'])
        tss = int(row['TSS'])
        adjusted_tss = int(row['adjusted_TSS'])
        start = min(tss, adjusted_tss)
        end = max(tss, adjusted_tss)
        tss_map[gene_id] = (chromosome, start, end)

    # 读取 stop 文件
    stop_df = pd.read_csv(stop_file, sep=r'\s+', engine='python')

    # 存储符合条件的结果
    matched_rows = []

    # 对每一行进行匹配检查
    for _, row in stop_df.iterrows():
        gene_name = row['gene_name']
        chr_pos = int(row['chr_pos'])
        chr_ev = str(row['chr'])

        if gene_name in tss_map:
            chromosome, start, end = tss_map[gene_name]
            if chromosome == chr_ev and start <= chr_pos <= end:
                matched_rows.append(row)

    # 写出结果
    if matched_rows:
        matched_df = pd.DataFrame(matched_rows)
        matched_df.to_csv(output_file, index=False)
        print(f"✅ 匹配完成，共找到 {len(matched_rows)} 条记录，已保存到 {output_file}")
    else:
        print(f"❌ 未找到匹配记录（{stop_file}）")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="匹配 EV_stop 文件与 TSS 区域")
    parser.add_argument("--tss", required=True, help="TSS 文件路径（CSV 格式）")
    parser.add_argument("--stop", required=True, help="EV_stop 文件路径（TXT 格式）")
    parser.add_argument("--output", help="输出文件路径（CSV）", default=None)

    args = parser.parse_args()

    main(args.tss, args.stop, args.output)
