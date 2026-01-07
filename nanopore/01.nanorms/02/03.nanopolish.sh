#! /bin/sh

ref="/home/chen/jiangjing/00.db/hg38/hg38-tRNAs.fa"
for d in /home/chen/jiangjing/02.nanoRMS/CTR1226CTR1 /home/chen/jiangjing/02.nanoRMS/CTR1226PIC1; do
  nanopolish index -d $d/fast5 $d/fastq/merge.fastq.gz 2> /dev/null;
  minimap2 -ax map-ont $ref $d/fastq/merge.fastq.gz 2> /dev/null | samtools sort -o $d/fastq/fq.gz.bam;
  samtools index $d/fastq/fq.gz.bam;
  nanopolish eventalign -n -t 6 -q 10 --progress --signal-index --scale-events --reads $d/fastq/merge.fastq.gz --bam $d/fastq/fq.gz.bam --genome $ref | gzip > $d/fastq/fq.gz.bam.events.gz;
done
