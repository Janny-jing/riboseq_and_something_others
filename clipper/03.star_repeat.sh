#! /bin/bash

index="/home/jj2024/00.db/01.hg38/15.clipper_repeat_rm"

for i in *_umi_trim.fq;do
        n=${i/_umi_trim.fq/}
        STAR --runMode alignReads --runThreadN 8 --genomeDir ${index} \
                --genomeLoad LoadAndRemove --readFilesIn ${n}_umi_trim.fq \
                --outSAMunmapped Within --outFilterMultimapNmax 30 --outFilterMultimapScoreRange 1 \
                --outFileNamePrefix ${n} \
                --outSAMattributes All --outStd BAM_Unsorted --outSAMtype BAM Unsorted \
                --outFilterType BySJout --outReadsUnmapped Fastx --outFilterScoreMin 10 --outSAMattrRGline ID:foo --alignEndsType EndToEnd
done

