#! /bin/sh

ref="/home/chen/jiangjing/00.db/hg38/hg38-tRNAs.fa"
/home/chen/jiangjing/soft/nanoRMS/per_read/get_features.py --rna -f $ref -t 6 -i /home/chen/jiangjing/02.nanoRMS/CTR*/fast5/*.fast5 
