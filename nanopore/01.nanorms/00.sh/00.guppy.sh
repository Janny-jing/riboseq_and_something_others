#! /bin/sh

/home/chen/jiangjing/soft/ont-guppy-cpu/bin/guppy_basecaller -i fast5 -s fastq -r --min_qscore 0 --flowcell FLO-MIN106 --kit SQK-RNA002
cat ./*fastq >all.fastq
seqkit seq all.fastq --rna2dna > alldna.fastq

