#! /bin/bash

#SBATCH --job-name=06
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --time=2-00:00:00



for i in *_sorted.bam;do
	n=${i/_sorted.bam/}
        umi_tools dedup -I ${n}_sorted.bam --output-stats=deduplicatd -S ${n}_rmRep_sorted.bam
done
