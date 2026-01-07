#! /bin/bash

#SBATCH --job-name=09
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --time=2-00:00:00


for i in *_c3.0_cond1.bed;do
    n=${i/_c3.0_cond1.bed/}
    /home/gongweikang/.conda/envs/chipseq/share/homer/bin/annotatePeaks.pl ${n}_c3.0_cond1.bed hg38 > ${n}_annotation.xls 
done
