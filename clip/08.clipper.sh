#! /bin/bash

#SBATCH --job-name=08
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --time=2-00:00:00



for i in IP*_unsorted_rmRep_gatksorted.bam;do
	n=${i/_rmRep_gatksorted.bam/}
        clipper -b ${n}_rmRep_gatksorted.bam -s hg19 -o ${n}.peaks.bed --processors 10
done
