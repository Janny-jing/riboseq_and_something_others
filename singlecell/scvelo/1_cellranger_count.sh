#!/bin/bash
#SBATCH -c 16
#SBATCH -J cellranger
#SBATCH --mem=64G
#SBATCH -p intel-sc3,amd-ep2
#SBATCH -q normal

source /home/liuxiaodongLab/jiangjing/miniconda3/etc/profile.d/conda.sh
conda activate scvelo

module load cellranger/7.1.0

for i in 201129A_P1*_fastqs;do
	n=${i/_fastqs/}
        cellranger count --id ${n} --transcriptome /storage/liuxiaodongLab/jiangjing/00.db/02.human/01.GRCh38_cellranger_ref/refdata-gex-GRCh38-2024-A/ --fastqs /storage/liuxiaodongLab/jiangjing/Projects/YutingFu/PD_YutingFu/scvelo/${n}_fastqs/  --sample ${n}  
done

