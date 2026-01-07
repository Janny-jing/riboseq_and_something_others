#!/bin/bash

# 设置参数
THREADS=6
MIN_LEN=20
MAX_LEN=40
ADAPTER_R1="TGGAATTCT"        # 根据实际试剂盒修改！

# 创建输出目录
mkdir -p cleaned reports

# 查找所有 _1.fq.gz 文件，自动匹配 _2
for read1 in *_3.fq.gz; do
    # 提取样本名（如 EV-1, orf3a-1）
    sample=$(basename "$read1" _3.fq.gz)
    fastp \
        -i "$read1" \
        -o "cleaned/${sample}_3.clean.fq.gz" \
        --json "reports/${sample}.fastp.json" \
        --html "reports/${sample}.fastp.html" \
        --report_title "Ribo-seq QC Report - $sample" \
        --adapter_sequence "$ADAPTER_R1" \
        --length_required "$MIN_LEN" \
        --length_limit "$MAX_LEN" \
        --thread "$THREADS" \
	--qualified_quality_phred 15 \
        --unqualified_percent_limit 50 \
        --n_base_limit 10

    if [[ $? -eq 0 ]]; then
        echo "✅ $sample processed successfully."
    else
        echo "❌ fastp failed for $sample"
    fi
done

echo "✅ All samples processed. Clean data in 'cleaned/' and reports in 'reports/'."
