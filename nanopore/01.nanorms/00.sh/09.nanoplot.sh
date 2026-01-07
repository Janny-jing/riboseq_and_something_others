#! /bin/sh

for i in *.fastq;do
	n=${i/.fastq/}
	mkdir ${n}_nanoplot
	NanoPlot -t4 --fastq ${n}.fastq -o ${n}_nanoplot
done
