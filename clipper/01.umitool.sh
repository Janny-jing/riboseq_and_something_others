#! /bin/bash

for i in *_1.fq.gz;do
        n=${i/_1.fq.gz/}
        umi_tools extract --stdin ${n}_1.fq.gz --bc-pattern NNNNNNNNNN --log ${n}_processed.log --stdout ${n}_umi.fq
done

