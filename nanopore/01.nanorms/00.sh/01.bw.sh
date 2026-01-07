#! /bin/sh

#ref="/home/chen/jiangjing/00.db/hg38/hg38_trna/hg38-tRNAs.fa"
ref="/home/chen/jiangjing/00.db/mm39/GtRNAdb_fa/mm39-tRNAs.fa"
for i in *.fastq.gz;do
       n=${i/.fastq.gz/}
       bwa mem -W13 -k6 -xont2d -T20 $ref $n.fastq.gz >$n.sam
       samtools view -b $n.sam >$n_1.bam
       samtools sort $n_1.bam > $n.bam
       samtools index $n.bam
       rm $n.sam $n_1.bam
done   
