#! /bin/sh
# ref="/home/jj2024/00.db/01.hg38/02.hg38_trna_mt/human_tRNA_mt.fa"
ref="/home/jj2024/00.db/01.hg38/08.hg38_trna_nomt_3adapter/human_tRNA_nomt_3adapter.fa"

for i in *.fq.gz;do
       n=${i/.fq.gz/}
       bwa mem -W13 -k6 -xont2d -T20 $ref ${n}.fq.gz |\
       samtools view -q1 -F4 -b - |\
       samtools sort - > ${n}_nomt.bam
       samtools index ${n}_nomt.bam
       samtools idxstats ${n}_nomt.bam >${n}_trna_stat_q1_nomt.txt
done


