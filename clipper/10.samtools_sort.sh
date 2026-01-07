#! /bin/sh
for i in *_unsorted.bam;do
        n=${i/_unsorted.bam/}
        samtools view -bF4 ${n}_unsorted.bam |\
        samtools sort - >${n}_mapped.bam
        samtools index ${n}_mapped.bam
done

