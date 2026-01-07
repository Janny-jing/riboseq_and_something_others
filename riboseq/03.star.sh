#!/bin/bash

# === 设置路径 ===
GENOME_INDEX="/home/jj2024/00.db/01.hg38/12.gencode_human_dna_starindx/index"
GTF_FILE="/home/jj2024/00.db/01.hg38/14.gencode_human_rRNA_Trna_ncRNA/gencode.v48.annotation.gtf"
THREADS=8
OUTPUT_DIR="../ribo_align"  # 可选：指定输出目录

# 创建输出目录
mkdir -p ${OUTPUT_DIR}

# === 循环处理所有 _clean_rmrRNA.fq.gz 文件 ===
for r1 in *_clean_rmrRNA.fq.gz; do
    # 提取样本名
    n=$(basename ${r1} _clean_rmrRNA.fq.gz)

    echo "🔍 处理样本: ${n}"

    # Step 1: STAR全基因组比对
    STAR --runThreadN ${THREADS} \
         --genomeDir ${GENOME_INDEX} \
         --sjdbGTFfile ${GTF_FILE} \
         --readFilesIn ${r1} \
         --readFilesCommand zcat \
         --twopassMode Basic \
         --outFileNamePrefix ${OUTPUT_DIR}/${n}_ribo_ \
         --outFilterType BySJout \
         --alignEndsType Local \
         --outSAMtype BAM SortedByCoordinate \
         --outSAMattributes All \
         --quantMode TranscriptomeSAM GeneCounts 

done

