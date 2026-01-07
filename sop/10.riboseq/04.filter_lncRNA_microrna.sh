#!/bin/bash
# simple_batch_filter.sh - 简化版批量过滤

# 创建蛋白编码基因GTF
#awk '$3 == "gene" && /gene_type "protein_coding"/' /home/jj2024/00.db/01.hg38/12.gencode_human_dna_starindx/gencode.v48.annotation.gtf > protein_coding_genes.gtf

# 批量处理
for bam in *_ribo_Aligned.sortedByCoord.out.bam; do
    sample=$(basename "$bam" "_ribo_Aligned.sortedByCoord.out.bam")
    echo "处理: $sample"
    
    bedtools intersect -a "$bam" -b protein_coding_genes.gtf -wa -u > "${sample}_filter.bam" 2>/dev/null
    samtools index "${sample}_filter.bam"
    
    echo "✅ 完成: ${sample}_filter.bam"
done

echo "🎉 所有样本处理完成！"
