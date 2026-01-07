#!/bin/bash
#SBATCH -c 4
#SBATCH -J cellranger
#SBATCH --mem=64G
#SBATCH -p intel-sc3,amd-ep2
#SBATCH -q normal

source /home/liuxiaodongLab/jiangjing/miniconda3/etc/profile.d/conda.sh
conda activate dnbc4tools

dnbc4tools rna multi --list sample.tsv --genomeDir --threads 4
