#! /bin/sh 

index="/home/jj2024/24.riboseq_clipper/01.ref"

for i in *.out.mate1;do
        n=${i/.out.mate1/}
        STAR --runMode alignReads --runThreadN 9 --genomeDir ${index}  \
             --genomeLoad LoadAndRemove --readFilesIn ${n}.out.mate1 ${n}.out.mate2 \
             --outSAMunmapped Within --outFilterMultimapNmax 1 --outFilterMultimapScoreRange 1 \
             --outFileNamePrefix ${n}_unsorted. \
             --outSAMattributes All --outStd BAM_Unsorted --outSAMtype BAM Unsorted \
             --outFilterType BySJout --outReadsUnmapped Fastx \
             --outFilterScoreMin 10 --outSAMattrRGline ID:foo --alignEndsType Local >${n}_unsorted.bam
done

