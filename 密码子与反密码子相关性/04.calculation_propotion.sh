#!/bin/bash

file="A549_RNASEQ_hg38_nameid_modify_ctl_calculation.txt"

awk '{
    sum[$1] += $2
    totalsum += $2
}
END {
    for (geneid in sum) {
        printf("%s\t%g\n", geneid, sum[geneid]/totalsum)
    }
}' "$file" |sort > codons_propotion.txt

awk '{sum[$1] += $2} END {for (geneid in sum){printf("%s\t%g\n",geneid,sum[geneid])}}' pic.txt >pic_sum.txt