#! /bin/sh

for i in *_rev.fa;do
    n=${i/_rev.fa/}
    seqkit seq -i ${n}_rev.fa | awk 'NR%2==0 {print substr($0, 2)} NR%2==1' >${n}_rev_1.fa
done
