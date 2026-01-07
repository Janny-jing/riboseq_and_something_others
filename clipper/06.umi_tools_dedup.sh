#! /bin/bash

for i in *_sorted.bam;do
        n=${i/_sorted.bam/}
        umi_tools dedup -I ${n}_sorted.bam --output-stats=deduplicatd -S ${n}_rmRep_sorted.bam
done

