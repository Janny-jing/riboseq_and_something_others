#! /bin/sh

ref="/home/chen/jiangjing/00.db/hg38/hg38-tRNAs.fa"
for i in *.bam;do
	python3 /home/chen/jiangjing/soft/nanoRMS/epinano_RMS/epinano_rms.py -R $ref -b $i -s /home/chen/jiangjing/soft/nanoRMS/epinano_RMS/sam2tsv.jar
done
