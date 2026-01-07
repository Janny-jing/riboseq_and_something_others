#! /bin/sh

ref="/home/chen/jiangjing/00.db/mm39/GtRNAdb_fa_dataprocessing_mt/mm39_rmdup_add_mt.fa"
for i in *.fq.gz;do
       n=${i/.fq.gz/}
       bwa mem -W13 -k6 -xont2d -T20 $ref ${n}.fq.gz |\
       samtools view -m80 -F4 -b - |\
       samtools sort - > ${n}_mt.bam
       samtools index ${n}_mt.bam
       samtools idxstats ${n}_mt.bam >${n}_trna_mt_stat.txt
done   

ref="/home/chen/jiangjing/00.db/mm39/GtRNAdb_fa_dataprocessing/mm39_rmdup_add.fa"
for i in *.fq.gz;do
       n=${i/.fq.gz/}
       bwa mem -W13 -k6 -xont2d -T20 $ref ${n}.fq.gz |\
       samtools view -m80 -F4 -b - |\
       samtools sort - > ${n}.bam
       samtools index ${n}.bam
       samtools idxstats ${n}.bam >${n}_trna_stat.txt
done

mkdir 02.mm39_mt_80 03.mm39_80
mv *_mt.bam* *_trna_mt_stat.txt 02.mm39_mt_80
mv *.bam* *_trna_stat.txt 03.mm39_80

ref="/home/chen/jiangjing/00.db/mm39/GtRNAdb_fa_dataprocessing/mm39_rmdup_add.fa"
for i in *.fq.gz;do
       n=${i/.fq.gz/}
       bwa mem -W13 -k6 -xont2d -T20 $ref ${n}.fq.gz |\
       samtools view -q1 -m80 -F4 -b - |\
       samtools sort - > ${n}_1.bam
       samtools index ${n}_1.bam
       samtools idxstats ${n}_1.bam >${n}_trna_stat_q0.txt
done

ref="/home/chen/jiangjing/00.db/mm39/GtRNAdb_fa_dataprocessing_mt/mm39_rmdup_add_mt.fa"
for i in *.fq.gz;do
       n=${i/.fq.gz/}
       bwa mem -W13 -k6 -xont2d -T20 $ref ${n}.fq.gz |\
       samtools view -q1 -m80 -F4 -b - |\
       samtools sort - > ${n}_mt_1.bam
       samtools index ${n}_mt_1.bam
       samtools idxstats ${n}_mt_1.bam >${n}_trna_mt_stat_q0.txt
done

ref="/home/chen/jiangjing/00.db/mm39/GtRNAdb_fa_dataprocessing/mm39_rmdup_add.fa"
for i in *.fq.gz;do
       n=${i/.fq.gz/}
       bwa mem -W13 -k6 -xont2d -T20 $ref $n.fq.gz |\
       samtools view -q5 -m80 -F4 -b - |\
       samtools sort - > ${n}_5.bam
       samtools index ${n}_5.bam
       samtools idxstats ${n}_5.bam >${n}_trna_stat_q5.txt
done

ref="/home/chen/jiangjing/00.db/mm39/GtRNAdb_fa_dataprocessing_mt/mm39_rmdup_add_mt.fa"
for i in *.fq.gz;do
       n=${i/.fq.gz/}
       bwa mem -W13 -k6 -xont2d -T20 $ref ${n}.fq.gz |\
       samtools view -q5 -m80 -F4 -b - |\
       samtools sort - > ${n}_mt_5.bam
       samtools index ${n}_mt_5.bam
       samtools idxstats ${n}_mt_5.bam >${n}_trna_mt_stat_q5.txt
done

mkdir 04.mm39_mt_q0 05.mm39_q0 06.mm39_mt_q5 07.mm39_q5
mv *_mt_1.bam* *_trna_mt_stat_q0.txt 04.mm39_mt_q0
mv *_1.bam* *_trna_stat_q0.txt 05.mm39_q0
mv *_mt_5.bam* *_trna_mt_stat_q5.txt 06.mm39_mt_q5
mv *_5.bam* *_trna_stat_q5.txt 07.mm39_q5

