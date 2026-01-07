#! /bin/bash

#SBATCH --job-name=03
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --time=4-00:00:00



index="/home/gongweikang/zhangxue_file/db/04.sc2_all"

for i in *.out.mate1;do
	n=${i/.out.mate1/}
	STAR --runMode alignReads --runThreadN 9 --genomeDir ${index}  \
	     --genomeLoad LoadAndRemove --readFilesIn ${n}.out.mate1 \
	     --outSAMunmapped Within --outFilterMultimapNmax 1 --outFilterMultimapScoreRange 1 \
	     --outFileNamePrefix ${n}_unsorted. \
	     --outSAMattributes All --outStd BAM_Unsorted --outSAMtype BAM Unsorted \
	     --outFilterType BySJout --outReadsUnmapped Fastx \
	     --outFilterScoreMin 10 --outSAMattrRGline ID:foo --alignEndsType Local >${n}_unsorted.bam
done
