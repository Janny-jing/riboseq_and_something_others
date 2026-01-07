#! /bin/sh

for i in IP*_umi.fq;do
        n=${i/_umi.fq/}
        cutadapt --match-read-wildcards --quality-cutoff 6 --times 1 -e 0.1 -O 1 -m 18 -a NNNNNAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -g CTTCCGATCTACAAGTT -g CTTCCGATCTTGGTCCT -o ${n}_umi_trim.fq ${n}_umi.fq > ${n}.metrics
done

