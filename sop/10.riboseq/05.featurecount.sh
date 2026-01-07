#!/bin/bash
# batch_featurecounts.sh - 批量基因表达定量

# === 设置参数 ===
GTF_FILE="/home/jj2024/00.db/01.hg38/14.gencode_human_rRNA_Trna_ncRNA/gencode.v48.annotation.gtf"
THREADS=8
OUTPUT_DIR="./featurecounts_results"

# 创建输出目录
mkdir -p ${OUTPUT_DIR}

echo "开始批量基因表达定量分析..."

# === 处理所有 *_filter.bam 文件 ===
for bam_file in *_filter.bam; do
    if [ ! -f "$bam_file" ]; then
        echo "未找到filtered BAM文件"
        continue
    fi

    # 提取样本名
    sample_name=$(basename "$bam_file" "_filter.bam")
    echo "处理样本: $sample_name"

    # 使用featureCounts计数，按gene_name分组
    featureCounts \
        -T ${THREADS} \
        -t CDS \                    # 计数CDS区域（推荐用于Ribo-seq）
        --extraAttributes gene_name \              # 按gene_name分组（更易读）
        -a ${GTF_FILE} \
        -o ${OUTPUT_DIR}/${sample_name}_gene_counts.txt \
        $bam_file

    if [ $? -eq 0 ]; then
        echo "✅ $sample_name 计数完成"
    else
        echo "❌ $sample_name 计数失败"
        continue
    fi

    # 生成FPKM标准化文件
    echo "生成FPKM标准化表达量..."
    python3 << EOF
import pandas as pd
import numpy as np

# 读取计数文件
counts_data = pd.read_csv("${OUTPUT_DIR}/${sample_name}_gene_counts.txt",
                         sep='\t', comment='#', header=0)
counts_data = counts_data.dropna()

# 提取基因长度和计数
gene_lengths = counts_data['Length'].values
counts = counts_data.iloc[:, -1].values  # 最后一列是计数

# 计算FPKM
total_counts = counts.sum()
fpkm = (counts * 1e9) / (gene_lengths * total_counts)

# 保存结果
result_df = pd.DataFrame({
    'Geneid': counts_data['Geneid'],
    'gene_name': counts_data['gene_name'],
    'Length': gene_lengths,
    'Counts': counts,
    'FPKM': fpkm
})

result_df.to_csv("${OUTPUT_DIR}/${sample_name}_fpkm.txt",
                sep='\t', index=False)
print("FPKM文件已保存: ${OUTPUT_DIR}/${sample_name}_fpkm.txt")
EOF

done

# === 合并所有样本的计数矩阵 ===
echo "合并所有样本的计数矩阵..."

# 合并原始计数
python3 << EOF
import pandas as pd
import glob

# 收集所有样本的计数
count_files = glob.glob("${OUTPUT_DIR}/*_gene_counts.txt")
count_matrix = None

for file in count_files:
    sample_name = file.split('/')[-1].replace('_gene_counts.txt', '')
    counts_data = pd.read_csv(file, sep='\t', comment='#', header=0)

    if count_matrix is None:
        count_matrix = counts_data[['Geneid', 'gene_name', 'Length']].copy()

    count_matrix[sample_name] = counts_data.iloc[:, -1]

# 保存计数矩阵
count_matrix.to_csv("${OUTPUT_DIR}/combined_count_matrix.txt", sep='\t', index=False)
print("计数矩阵已保存: ${OUTPUT_DIR}/combined_count_matrix.txt")

# 合并FPKM矩阵
fpkm_files = glob.glob("${OUTPUT_DIR}/*_fpkm.txt")
fpkm_matrix = None

for file in fpkm_files:
    sample_name = file.split('/')[-1].replace('_fpkm.txt', '')
    fpkm_data = pd.read_csv(file, sep='\t')

    if fpkm_matrix is None:
        fpkm_matrix = fpkm_data[['Geneid', 'gene_name']].copy()

    fpkm_matrix[sample_name] = fpkm_data['FPKM']

fpkm_matrix.to_csv("${OUTPUT_DIR}/combined_fpkm_matrix.txt", sep='\t', index=False)
print("FPKM矩阵已保存: ${OUTPUT_DIR}/combined_fpkm_matrix.txt")
EOF

echo "🎉 批量featureCounts分析完成！"
echo "输出文件:"
echo "  - 单个样本计数: ${OUTPUT_DIR}/*_gene_counts.txt"
echo "  - 单个样本FPKM: ${OUTPUT_DIR}/*_fpkm.txt"
echo "  - 合并计数矩阵: ${OUTPUT_DIR}/combined_count_matrix.txt"
echo "  - 合并FPKM矩阵: ${OUTPUT_DIR}/combined_fpkm_matrix.txt"
