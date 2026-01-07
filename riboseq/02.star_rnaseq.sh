#! /bin/sh

GENOME_INDEX="/home/jj2024/00.db/01.hg38/12.gencode_human_dna_starindx/index"

for i in *_1_val_1.fq.gz;do
    n=${i/_1_val_1.fq.gz/}
    STAR  --runThreadN 8 \
          --genomeDir ${GENOME_INDEX} \
          --readFilesIn ${n}_1_val_1.fq.gz ${n}_2_val_2.fq.gz \
          --readFilesCommand zcat \
          --outFileNamePrefix ${n} \
          --outFilterType BySJout  --alignEndsType EndToEnd --outSAMtype BAM SortedByCoordinate 
done
