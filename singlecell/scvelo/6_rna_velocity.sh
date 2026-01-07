#!/bin/bash
#SBATCH -c 4
#SBATCH -J velocyto
#SBATCH --mem=128G
#SBATCH -p intel-sc3,amd-ep2
#SBATCH -q normal

source /home/liuxiaodongLab/jiangjing/miniconda3/etc/profile.d/conda.sh
conda activate diffusionmap

python rna_velocity_1.py
