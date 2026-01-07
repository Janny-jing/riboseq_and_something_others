#! /bin/sh

ref="/home/chen/jiangjing/00.db/mm39/GtRNAdb_fa_dataprocessing/mm39_rmdup_add.fasta"

for i in *.bam;do
	python3 /home/chen/jiangjing/soft/nanoRMS/epinano_RMS/epinano_rms.py -R $ref -b $i -s /home/chen/jiangjing/soft/nanoRMS/epinano_RMS/sam2tsv.jar
done
