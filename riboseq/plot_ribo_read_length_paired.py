#!/usr/bin/env python
# plot_ribo_read_length_paired.py
# 功能：统计双端 ribo-seq 数据中 R1 和 R2 的 read 长度分布，并绘制柱状图

import gzip
from collections import Counter
import matplotlib.pyplot as plt
import numpy as np

# ================= 参数设置 =================
# 样本：每个样本包含多个重复，每个重复有 R1 和 R2
SAMPLES = {
    "EV": [
        ("EV-1_clean_R1.fq", "EV-1_clean_R2.fq"),
        ("EV-2_clean_R1.fq", "EV-2_clean_R2.fq")
    ],
    "orf3a": [
        ("orf3a-1_clean_R1.fq", "orf3a-1_clean_R2.fq"),
        ("orf3a-2_clean_R1.fq", "orf3a-2_clean_R2.fq")
    ]
}

MIN_LEN = 25
MAX_LEN = 40
OUTPUT_FILE = "ribo_read_length_distribution_paired.png"
# ============================================

def read_fastq(fastq_file):
    """读取 FASTQ 文件，返回序列长度生成器"""
    open_func = gzip.open if fastq_file.endswith(".gz") else open
    with open_func(fastq_file, "rt") as f:
        for line in f:
            if line.startswith("@"):  # read header
                seq = next(f).strip()
                yield len(seq)
                # skip + and quality
                try:
                    next(f)  # skip +
                    next(f)  # skip quality
                except StopIteration:
                    break

def collect_length_counts(paired_files_list):
    """收集一个样本所有重复的 R1 + R2 长度计数"""
    total_counts = Counter()
    for r1_file, r2_file in paired_files_list:
        print(f"Processing R1: {r1_file}")
        for length in read_fastq(r1_file):
            if MIN_LEN <= length <= MAX_LEN:
                total_counts[length] += 1

        print(f"Processing R2: {r2_file}")
        for length in read_fastq(r2_file):
            if MIN_LEN <= length <= MAX_LEN:
                total_counts[length] += 1

    return total_counts

# ================= 主程序 =================
if __name__ == "__main__":
    plt.figure(figsize=(10, 6))

    bar_width = 0.4
    opacity = 0.8
    index = 0
    colors = ['skyblue', 'salmon']

    for idx, (sample_name, paired_files) in enumerate(SAMPLES.items()):
        length_counter = collect_length_counts(paired_files)

        lengths = sorted(length_counter.keys())
        counts = [length_counter[l] for l in lengths]

        # 错开柱子位置
        x_pos = [x + idx * bar_width for x in lengths]
        plt.bar(x_pos, counts, bar_width, alpha=opacity, color=colors[idx], label=sample_name)

    # 设置横轴刻度在中间
    all_lengths = list(range(MIN_LEN, MAX_LEN + 1))
    plt.xlabel("Read Length (nt)")
    plt.ylabel("Number of Reads")
    plt.title("Ribo-seq Read Length Distribution (R1 + R2, after rRNA removal)")
    plt.xticks([x + bar_width/2 for x in all_lengths], all_lengths, rotation=0)
    plt.legend()
    plt.grid(axis='y', alpha=0.3)
    plt.tight_layout()
    plt.savefig(OUTPUT_FILE, dpi=200)
    plt.show()

    print(f"Plot saved to {OUTPUT_FILE}")
