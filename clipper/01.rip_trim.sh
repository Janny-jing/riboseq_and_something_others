#! /bin/sh

for i in *_1.fq.gz;do
    n=${i/_1.fq.gz/}
	 trim_galore --paired --fastqc --adapter 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC' --adapter2 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT' --quality 20 ${n}_1.fq.gz ${n}_2.fq.gz
done

