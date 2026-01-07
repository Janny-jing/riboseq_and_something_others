#!/usr/bin/env bash

#SBATCH -J scPoli
#SBATCH -o scPoli.%j.out
#SBATCH -e scPoli.%j.err
#SBATCH -p intel-fat
#SBATCH -q hmem
#SBATCH -c 8           #这里需要填上CPU数量
#SBATCH --mem=480G      #这里需要填上内存需求


source /home/liuxiaodongLab/jiangjing/miniconda3/etc/profile.d/conda.sh
conda activate /storage2/liuxiaodongLab/jiangjing/miniconda3/envs/agent_new
module load gcc/11.2.0

# Run label transfer
python scPoli_without_invitro_combineddata_examine_count_lineage.py --cell_type_key lineage_pred --n_top_genes 2000 \
	--model_dir /storage2/liuxiaodongLab/jiangjing/Projects/XueyingFan/PD_XueyingFan/20251205_hypoblast_in_house_d20_d26/model \
        --data_path /storage2/liuxiaodongLab/jiangjing/Projects/XueyingFan/PD_XueyingFan/20251205_hypoblast_in_house_d20_d26/code/adata_human_qc_processed.h5ad
