#! /bin/bash

index="/home/jj2024/00.db/01.hg38/15.clipper_repeat_rm"

for i in *_1_val_1.fq.gz;do
        n=${i/_1_val_1.fq.gz/}
        STAR --runMode alignReads --runThreadN 8 --genomeDir ${index} --readFilesCommand zcat \
                --genomeLoad LoadAndRemove --readFilesIn ${n}_1_val_1.fq.gz ${n}_2_val_2.fq.gz \
                --outSAMunmapped Within --outFilterMultimapNmax 30 --outFilterMultimapScoreRange 1 \
                --outFileNamePrefix ${n} \
                --outSAMattributes All --outSAMtype BAM Unsorted \
                --outFilterType BySJout --outReadsUnmapped Fastx --outFilterScoreMin 10 --outSAMattrRGline ID:foo --alignEndsType EndToEnd
done

