#!/usr/bin/env bash

#SBATCH -J scPoli
#SBATCH -o scPoli.%j.out
#SBATCH -e scPoli.%j.err
#SBATCH -p intel-fat
#SBATCH -q hmem
#SBATCH -c 8           #这里需要填上CPU数量
#SBATCH --mem=480G       #这里需要填上内存需求


source /home/liuxiaodongLab/jiangjing/miniconda3/etc/profile.d/conda.sh
conda activate /storage2/liuxiaodongLab/jiangjing/miniconda3/envs/agent_new

# Run label transfer
python label_transfer_color_v2.py \
      --file_path /storage2/liuxiaodongLab/jiangjing/Projects/XueyingFan/PD_XueyingFan/20251106_fetal_adult_merge_name_scdevelopment/output/validation_100k_without_lineage_preprocessed.h5ad \
      --figures_folder /storage2/liuxiaodongLab/jiangjing/Projects/XueyingFan/PD_XueyingFan/20251106_fetal_adult_merge_name_scdevelopment/label_transfer_model_text/reanno_embryo_hvg4000/ \
      --custom_model_dir /storage2/liuxiaodongLab/fanxueying/developmental_atlas/code/20251113_scpoli_brain_all_reanno_training/label_transfer_reanno_20251113_v3/scpoli_model_reanno_hvg4000 \
      --custom_adata_path /storage2/liuxiaodongLab/fanxueying/developmental_atlas/code/20251113_scpoli_brain_all_reanno_training/label_transfer_reanno_20251113_v3/scpoli_model_reanno_hvg4000/adata.h5ad \
      --model_type cell_type      

echo "Label transfer completed"
