#! /bin/sh

# ref="/home/chen/jiangjing/00.db/mm39/GtRNAdb_fa_dataprocessing/mm39_rmdup_add.fasta"
ref="/home/chen/jiangjing/00.db/hg38/hg38_trna_dataprocessing/hg38_rmdup_add.fasta"
for i in *.fastq;do
       n=${i/.fastq/}
       bwa mem -W13 -k6 -xont2d -T20 $ref $n.fastq >$n.sam
       samtools view -m80 -F4 -b $n.sam >$n_1.bam
       samtools sort $n_1.bam > $n.bam
       samtools index $n.bam
       rm $n.sam  $n_1.bam
       samtools idxstats ${n}.bam >${n}_trna_stat.txt
done   
