#! /bin/bash

#SBATCH --job-name=getfasta
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=2:00:00


for i in *.bed;do
    n=${i/.bed/}
    bedtools getfasta -fi /home/gongweikang/zhangxue_file/db/Homo_sapiens.GRCh38.dna.toplevel.fa -bed ${n}.bed -fo ${n}.fa -nameOnly
done
