#! /bin/sh

#ref="/home/chen/jiangjing/00.db/hg38/hg38_trna/hg38-tRNAs.fa"
ref="/home/chen/jiangjing/00.db/mm39/GtRNAdb_fa/mm39-tRNAs.fa"

for i in *.bam;do
	python3 /home/chen/jiangjing/soft/nanoRMS/epinano_RMS/epinano_rms.py -R $ref -b $i -s /home/chen/jiangjing/soft/nanoRMS/epinano_RMS/sam2tsv.jar
done
