#! /bin/sh

for i in *.fq;do
        n=${i/.fq/}
        seqkit seq -Q 6 -m 80 ${n}.fq | gzip -c >${n}.fq.gz
done

