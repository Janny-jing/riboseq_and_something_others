#!/usr/bin/env bash

#SBATCH -J scPoli
#SBATCH -o scPoli.%j.out
#SBATCH -e scPoli.%j.err
#SBATCH -p intel-sc3,amd-ep2
#SBATCH -q normal
#SBATCH -c 8           #这里需要填上CPU数量
#SBATCH --mem=64G       #这里需要填上内存需求


source /home/liuxiaodongLab/jiangjing/miniconda3/etc/profile.d/conda.sh
conda activate /storage2/liuxiaodongLab/jiangjing/miniconda3/envs/agent_new
module load gcc/11.2.0

python integrate_lineage_to_h5ad_adata.py --csv_folder /storage2/liuxiaodongLab/jiangjing/Projects/XueyingFan/PD_XueyingFan/20251106_fetal_adult_merge_name_scdevelopment/output/model_validate_csv/lineage \
	--original_adata /storage2/liuxiaodongLab/jiangjing/Projects/XueyingFan/PD_XueyingFan/20251106_fetal_adult_merge_name_scdevelopment/output/validation_100k_with_lineage.h5ad \
	--output_adata /storage2/liuxiaodongLab/jiangjing/Projects/XueyingFan/PD_XueyingFan/20251106_fetal_adult_merge_name_scdevelopment/output/model_validate_csv/lineage/validation_100k_with_lineage.h5ad
