#!/usr/bin/env bash

#SBATCH -J preprocess
#SBATCH -o preprocess.%j.out
#SBATCH -e preprocess.%j.err
#SBATCH -p intel-fat
#SBATCH -q hmem
#SBATCH -c 2           #这里需要填上CPU数量
#SBATCH --mem=480G       #这里需要填上内存需求

source /home/liuxiaodongLab/jiangjing/miniconda3/etc/profile.d/conda.sh
conda activate /storage2/liuxiaodongLab/jiangjing/miniconda3/envs/agent_new
module load gcc/11.2.0

# Run label transfer
 python preprocess_data.py /storage2/liuxiaodongLab/jiangjing/Projects/XueyingFan/PD_XueyingFan/20251106_fetal_adult_merge_name_scdevelopment/data/3sample_smalldataset_test_script_model_nolineage.h5ad \
        --output_dir /storage2/liuxiaodongLab/jiangjing/Projects/XueyingFan/PD_XueyingFan/20251106_fetal_adult_merge_name_scdevelopment/output 
        
echo "preprocess completed"
