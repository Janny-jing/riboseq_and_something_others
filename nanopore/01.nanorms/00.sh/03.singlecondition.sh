#! /bin/sh

for i in *.per.site.baseFreq.csv;do
	Rscript --vanilla /home/chen/jiangjing/02.nanoRMS/00.sh/Pseudou_prediction_singlecondition.R -f $i -p /home/chen/jiangjing/02.nanoRMS/00.sh/single.tsv
done
