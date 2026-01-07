#! /bin/sh

index="/home/gongweikang/zhangxue_file/db/06.ucsc_hg19"

for i in *.out.mate1;do
	n=${i/.out.mate1/}
	STAR --runMode alignReads --runThreadN 9 --genomeDir ${index}  \
	     --genomeLoad LoadAndRemove --readFilesIn ${n}.out.mate1 \
	     --outSAMunmapped Within --outFilterMultimapNmax 1 --outFilterMultimapScoreRange 1 \
	     --outFileNamePrefix ${n}_unsorted. \
	     --outSAMattributes All --outStd BAM_Unsorted --outSAMtype BAM Unsorted \
	     --outFilterType BySJout --outReadsUnmapped Fastx \
	     --outFilterScoreMin 10 --outSAMattrRGline ID:foo --alignEndsType EndToEnd >${n}_unsorted.bam
done
