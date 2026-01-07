#! /bin/sh

for i in *.fq;do
	n=${i/.fq/}
	mkdir ${n}_nanoplot
	NanoPlot -t4 --fastq ${n}.fq -o ${n}_nanoplot
done
