#! /bin/sh

path="/home/chen/jiangjing/02.nanoRMS/"
ref="/home/chen/jiangjing/00.db/hg38/hg38-tRNAs.fa"
for i in CTR1226CTR1 CTR1226PIC1;do
       bwa mem -W13 -k6 -xont2d -T20 $ref $path/$i/fastq/merge.fastq.gz >$i.sam
       samtools view -b $i.sam >$i_1.bam
       samtools sort $i_1.bam > ${i}.bam
       samtools index ${i}.bam
       rm $i.sam $i_1.bam
done   
