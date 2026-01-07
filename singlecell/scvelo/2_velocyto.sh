#!/bin/bash
#SBATCH -c 16
#SBATCH -J velocyto
#SBATCH --mem=64G
#SBATCH -p intel-sc3,amd-ep2
#SBATCH -q normal

source /home/liuxiaodongLab/jiangjing/miniconda3/etc/profile.d/conda.sh
conda activate scvelo

cellranger_ref="/storage/liuxiaodongLab/jiangjing/00.db/02.human/01.GRCh38_cellranger_ref/refdata-gex-GRCh38-2024-A//genes/genes.gtf"
cellranger_outdir="/storage/liuxiaodongLab/jiangjing/Projects/YutingFu/PD_YutingFu/scvelo/201129A_P1_A" #cellranger output

velocyto run10x  $cellranger_outdir $cellranger_ref
velocyto run10x /storage/liuxiaodongLab/jiangjing/Projects/YutingFu/PD_YutingFu/scvelo/201129A_P1_B /storage/liuxiaodongLab/jiangjing/00.db/02.human/01.GRCh38_cellranger_ref/refdata-gex-GRCh38-2024-A//genes/genes.gtf
