#!/bin/sh

index="/home/jj2024/00.db/01.hg38/14.gencode_human_rRNA_Trna_ncRNA/rfam_human_ncrna/Homo_sapiens_all_ncRNA_index"

for i in *_1.clean.fq.gz; do
    j=${i/_1.clean.fq.gz/}
    echo "Processing sample: $j"

    # 检查输入文件是否存在
    if [ ! -f "${j}_1.clean.fq.gz" ]; then
        echo "Error: Input file ${j}_1.clean.fq.gz not found"
        continue
    fi

    # Step 1: Bowtie alignment against rRNA/ncRNA, output unmapped reads directly
    if ! bowtie2 -x "$index" \
        -U "${j}_1.clean.fq.gz" \
        --un-gz "${j}_clean_rmrRNA.fq.gz" \
        --al-gz "${j}_rRNA_mapped.fq.gz" \
        -S "${j}.sam" 2>"${j}.bowtie2.log"; then
        echo "Error: Bowtie2 alignment failed for $j"
        continue
    fi

    # 可选：计算rRNA去除率用于质控
    original_reads=$(zcat "${j}_1.clean.fq.gz" | echo $((`wc -l`/4)))
    rrna_reads=$(zcat "${j}_rRNA_mapped.fq.gz" 2>/dev/null | echo $((`wc -l`/4)) || echo 0)
    
    if [ $original_reads -gt 0 ]; then
        rrna_percentage=$(echo "scale=2; $rrna_reads * 100 / $original_reads" | bc)
        echo "rRNA removal stats: $rrna_reads/$original_reads reads mapped to rRNA (${rrna_percentage}%)"
    fi

    # 清理中间SAM文件
    rm -f "${j}.sam"

    echo "Completed processing sample: $j"
    echo "rRNA removed reads saved to: ${j}_clean_rmrRNA.fq.gz"

done
