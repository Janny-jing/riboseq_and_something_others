#!/bin/bash
#SBATCH -c 1
#SBATCH -J cellranger
#SBATCH --mem=64G
#SBATCH -p intel-sc3,amd-ep2
#SBATCH -q normal

source /home/liuxiaodongLab/jiangjing/miniconda3/etc/profile.d/conda.sh
conda activate diffusionmap

module load cellranger/7.1.0


for main_dir in /storage/liuxiaodongLab/jiangjing/Projects/YutingFu/PD_YutingFu/scvelo/W*; do
    for subdir in "$main_dir"/*_L01_QD*; do
	    # 遍历每个子目录中的fq.gz文件
            for fq1 in "$subdir"/*_1.fq.gz; do
                # 配对对应的R2文件
                fq2="${fq1/_1.fq.gz/_2.fq.gz}"

                # 检查文件是否真实存在
                if [ -f "$fq2" ]; then
                    # 提取样本编号（例如67/49/40等）
                    sample_id=$(basename "$fq1" | awk -F'_' '{print $3}')

                    # 打印文件路径用于验证
                    echo "Processing sample $sample_id:"
                    echo "  R1: $fq1"
                    echo "  R2: $fq2"

                    # 这里可以添加实际处理命令，例如：
                    # fastqc -o output_dir "$fq1" "$fq2"
                else
                    echo "WARNING: Missing R2 file for $fq1"
                fi
            done
      done
done
