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
      --file_path /storage2/liuxiaodongLab/jiangjing/Projects/XueyingFan/PD_XueyingFan/20251117_bgi_second_dataset/code/combined_with_first_second_ips_dpc_npc_hvg2000.h5ad \
      --figures_folder /storage2/liuxiaodongLab/jiangjing/Projects/XueyingFan/PD_XueyingFan/20251117_bgi_second_dataset/label_transfer_hvg2000/ \
      --custom_model_dir /storage2/liuxiaodongLab/jiangjing/Projects/XueyingFan/PD_XueyingFan/20251106_fetal_adult_merge_name_scdevelopment/model_defined_number/without_invitro_embryo_fetal_hvg2000_merge_adult_dims50_ep50/scpoli_model_lineage_hvg2000  \
      --custom_adata_path /storage2/liuxiaodongLab/jiangjing/Projects/XueyingFan/PD_XueyingFan/20251106_fetal_adult_merge_name_scdevelopment/model_defined_number/without_invitro_embryo_fetal_hvg2000_merge_adult_dims50_ep50/scpoli_model_lineage_hvg2000/adata.h5ad \
      #--model_type cell_type

echo "Label transfer completed"
