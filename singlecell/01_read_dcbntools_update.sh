#!/bin/bash
#SBATCH -c 2
#SBATCH -J dnbctools_main
#SBATCH --mem=10G
#SBATCH -p intel-sc3,amd-ep2
#SBATCH -q normal

source /storage2/liuxiaodongLab/jiangjing/02.soft/dnbc4tools2.1.3/sourceC4.bash

index="/storage2/liuxiaodongLab/jiangjing/00.db/02.human/03.GRCH38_dnbc4tools_ref_10x"

# 默认路径（当未指定命令行参数时使用）
default_data_dir="/storage2/liuxiaodongLab/jiangjing/Projects/XueyingFan/PD_XueyingFan/20251217_bgi_forth_dataset/data"
default_output_base="/storage2/liuxiaodongLab/jiangjing/Projects/XueyingFan/PD_XueyingFan/20251217_bgi_forth_dataset/output"

# 解析命令行参数
while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--input)
            data_dir="$2"
            shift 2
            ;;
        -o|--output)
            output_base="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            echo "Usage: $0 [-i|--input data_dir] [-o|--output output_dir]"
            exit 1
            ;;
    esac
done

# 使用命令行参数或默认值
data_dir="${data_dir:-$default_data_dir}"
output_base="${output_base:-$default_output_base}"
scripts_dir="${output_base}/job_scripts"

# 检查数据目录是否存在
if [[ ! -d "$data_dir" ]]; then
    echo "Error: Data directory does not exist: $data_dir"
    exit 1
fi

# 创建输出目录和脚本目录
mkdir -p "$output_base" "$scripts_dir"

echo "======================================"
echo "Data directory: $data_dir"
echo "Output directory: $output_base"
echo "======================================"

echo "Scanning for CDNA and Oligo sample pairs in: $data_dir"

# 显示目录结构
echo "=== Directory structure ==="
find "$data_dir" -maxdepth 1 -type d | sort

echo ""
echo "=== Finding sample pairs (case-insensitive) ==="

# 查找所有CDNA和Oligo目录（不区分大小写）
declare -A directories

# 首先收集所有相关目录
for dir in "$data_dir"/*; do
    if [[ -d "$dir" ]]; then
        dir_name=$(basename "$dir")
        # 转换为小写进行匹配
        dir_lower=$(echo "$dir_name" | tr '[:upper:]' '[:lower:]')
        
        if [[ "$dir_lower" == *cdna ]]; then
            # 提取样本基础名称（去掉后缀，不区分大小写）
            sample_base=$(echo "$dir_name" | sed 's/[-_]cdna$//i')
            directories["${sample_base}_cdna"]="$dir"
        elif [[ "$dir_lower" == *oligo ]]; then
            # 提取样本基础名称（去掉后缀，不区分大小写）
            sample_base=$(echo "$dir_name" | sed 's/[-_]oligo$//i')
            directories["${sample_base}_oligo"]="$dir"
        fi
    fi
done

# 配对所有找到的目录
declare -A sample_pairs
matched_samples=0

# 查找所有唯一的样本基础名称
declare -A unique_samples
for key in "${!directories[@]}"; do
    sample_base=$(echo "$key" | sed 's/_[^_]*$//')
    unique_samples["$sample_base"]=1
done

# 为每个样本基础名称查找配对的目录
for sample_base in "${!unique_samples[@]}"; do
    cdna_key="${sample_base}_cdna"
    oligo_key="${sample_base}_oligo"
    
    cdna_dir="${directories[$cdna_key]}"
    oligo_dir="${directories[$oligo_key]}"
    
    if [[ -n "$cdna_dir" && -n "$oligo_dir" ]]; then
        # 获取目录的实际名称（保留原始大小写）
        cdna_name=$(basename "$cdna_dir")
        oligo_name=$(basename "$oligo_dir")
        
        sample_pairs["$sample_base"]="$cdna_dir:$oligo_dir"
        ((matched_samples++))
        echo "Found pair: $sample_base"
        echo "  CDNA: $cdna_name"
        echo "  Oligo: $oligo_name"
    elif [[ -n "$cdna_dir" ]]; then
        echo "Warning: CDNA directory found but no matching Oligo directory for: $sample_base"
        echo "  CDNA: $(basename "$cdna_dir")"
    elif [[ -n "$oligo_dir" ]]; then
        echo "Warning: Oligo directory found but no matching CDNA directory for: $sample_base"
        echo "  Oligo: $(basename "$oligo_dir")"
    fi
done

if [ ${#sample_pairs[@]} -eq 0 ]; then
    echo "Error: No complete sample pairs found!"
    echo "Available directories in $data_dir:"
    find "$data_dir" -maxdepth 1 -type d | sort
    echo ""
    echo "Expected directory naming pattern:"
    echo "  SampleName-CDNA/  and  SampleName-Oligo/"
    echo "  or"
    echo "  SampleName-cdna/  and  SampleName-oligo/"
    exit 1
fi

echo ""
echo "Found ${#sample_pairs[@]} complete sample pairs"
echo "Generating job scripts..."

# 为每个样本生成独立的作业脚本
job_scripts=()

for sample_base in "${!sample_pairs[@]}"; do
    echo "======================================"
    echo "Generating script for sample: $sample_base"

    # 解析目录路径
    IFS=':' read -r CDNA_dir Oligo_dir <<< "${sample_pairs[$sample_base]}"

    # 直接在目录中查找文件
    CDNA_r1=""
    CDNA_r2=""
    Oligo_r1=""
    Oligo_r2=""

    # 扫描CDNA目录中的文件（不区分_fq.gz或.fastq.gz）
    for file in "$CDNA_dir"/*; do
        if [[ -f "$file" ]]; then
            filename=$(basename "$file")
            # 转换为小写进行匹配
            filename_lower=$(echo "$filename" | tr '[:upper:]' '[:lower:]')
            
            if [[ "$filename_lower" == *_1.fq.gz || "$filename_lower" == *_1.fastq.gz || 
                  "$filename_lower" == *r1*.fq.gz || "$filename_lower" == *r1*.fastq.gz ]]; then
                CDNA_r1="$file"
            elif [[ "$filename_lower" == *_2.fq.gz || "$filename_lower" == *_2.fastq.gz ||
                    "$filename_lower" == *r2*.fq.gz || "$filename_lower" == *r2*.fastq.gz ]]; then
                CDNA_r2="$file"
            fi
        fi
    done

    # 扫描Oligo目录中的文件（不区分_fq.gz或.fastq.gz）
    for file in "$Oligo_dir"/*; do
        if [[ -f "$file" ]]; then
            filename=$(basename "$file")
            # 转换为小写进行匹配
            filename_lower=$(echo "$filename" | tr '[:upper:]' '[:lower:]')
            
            if [[ "$filename_lower" == *_1.fq.gz || "$filename_lower" == *_1.fastq.gz ||
                  "$filename_lower" == *r1*.fq.gz || "$filename_lower" == *r1*.fastq.gz ]]; then
                Oligo_r1="$file"
            elif [[ "$filename_lower" == *_2.fq.gz || "$filename_lower" == *_2.fastq.gz ||
                    "$filename_lower" == *r2*.fq.gz || "$filename_lower" == *r2*.fastq.gz ]]; then
                Oligo_r2="$file"
            fi
        fi
    done

    # 检查文件是否存在
    if [[ -z "$CDNA_r1" || -z "$CDNA_r2" ]]; then
        echo "Error: Missing CDNA files for $sample_base in $CDNA_dir, skipping..."
        echo "  Available files:"
        ls -1 "$CDNA_dir"/*.fq.gz "$CDNA_dir"/*.fastq.gz 2>/dev/null || echo "    No fq.gz or fastq.gz files found"
        continue
    fi
    
    if [[ -z "$Oligo_r1" || -z "$Oligo_r2" ]]; then
        echo "Error: Missing Oligo files for $sample_base in $Oligo_dir, skipping..."
        echo "  Available files:"
        ls -1 "$Oligo_dir"/*.fq.gz "$Oligo_dir"/*.fastq.gz 2>/dev/null || echo "    No fq.gz or fastq.gz files found"
        continue
    fi

    # 创建样本特定的作业脚本
    job_script="${scripts_dir}/${sample_base}_dnbc4tools.sh"

    cat > "$job_script" << EOF
#!/bin/bash
#SBATCH -c 10
#SBATCH -J dnb_${sample_base}
#SBATCH --mem=256G
#SBATCH -p intel-sc3,amd-ep2
#SBATCH -q normal
#SBATCH -o ${output_base}/${sample_base}_%j.out
#SBATCH -e ${output_base}/${sample_base}_%j.err

source /storage2/liuxiaodongLab/jiangjing/02.soft/dnbc4tools2.1.3/sourceC4.bash

echo "Starting processing for sample: $sample_base"
echo "Start time: \$(date)"
echo "Data directory: $data_dir"
echo "Output directory: $output_base"

# 创建输出目录
sample_output="${output_base}/${sample_base}_sample"
mkdir -p "\$sample_output"

echo "Input files:"
echo "  CDNA R1: $(basename "$CDNA_r1")"
echo "  CDNA R2: $(basename "$CDNA_r2")"
echo "  Oligo R1: $(basename "$Oligo_r1")"
echo "  Oligo R2: $(basename "$Oligo_r2")"

# 运行dnbc4tools
dnbc4tools rna run \\
    --cDNAfastq1 "$CDNA_r1" \\
    --cDNAfastq2 "$CDNA_r2" \\
    --oligofastq1 "$Oligo_r1" \\
    --oligofastq2 "$Oligo_r2" \\
    --genomeDir $index \\
    --name "${sample_base}_sample" \\
    --threads 10 \\
    --outdir "\$sample_output"

exit_code=\$?
echo "End time: \$(date)"

if [ \$exit_code -eq 0 ]; then
    echo "Successfully completed processing sample: $sample_base"
else
    echo "Error processing sample: $sample_base (exit code: \$exit_code)"
fi

exit \$exit_code
EOF

    # 使脚本可执行
    chmod +x "$job_script"
    job_scripts+=("$job_script")
    echo "Generated job script: $job_script"
done

echo ""
echo "Generated ${#job_scripts[@]} job scripts"

if [ ${#job_scripts[@]} -eq 0 ]; then
    echo "Error: No valid job scripts were generated!"
    exit 1
fi

# 提交所有作业
echo ""
echo "Submitting jobs to SLURM..."
job_ids=()

for job_script in "${job_scripts[@]}"; do
    job_name=$(basename "$job_script" .sh)
    echo "Submitting job: $job_name"
    job_id=$(sbatch "$job_script" | awk '{print $4}')
    if [[ -n "$job_id" ]]; then
        job_ids+=("$job_id")
        echo "  Submitted with job ID: $job_id"
    else
        echo "  Failed to submit job: $job_name"
    fi
done

echo ""
echo "======================================"
echo "Submitted ${#job_ids[@]} jobs successfully!"
echo "Job IDs: ${job_ids[*]}"
echo ""
echo "You can check job status using:"
echo "  squeue -u \$USER"
echo ""
echo "Job outputs will be saved in: $output_base"
echo "Individual sample outputs in: ${output_base}/{sample}_sample/"
echo "Job scripts are saved in: $scripts_dir"

# 生成一个监控脚本
monitor_script="${output_base}/monitor_jobs.sh"
cat > "$monitor_script" << 'EOF'
#!/bin/bash
echo "=== Job Status Monitor ==="
echo "Current time: $(date)"
echo ""
echo "Active jobs:"
squeue -u $USER -o "%.10i %.20j %.10u %.8T %.10M %.6D %.20R %b" | grep -E "(dnb_|JOBID)"
echo ""
echo "Completed jobs:"
if [[ -n "$(sacct -u $USER --format=JobID,JobName,State,ExitCode,Elapsed -X | grep dnb_)" ]]; then
    sacct -u $USER --format=JobID,JobName,State,ExitCode,Elapsed -X | grep dnb_
else
    echo "No completed dnbctools jobs found"
fi
EOF

chmod +x "$monitor_script"
echo ""
echo "Monitor script created: $monitor_script"
echo "Run './$monitor_script' to check job status"

# 生成使用说明
echo ""
echo "======================================"
echo "USAGE SUMMARY:"
echo "======================================"
echo "To run with custom directories:"
echo "  sbatch $0 -i /path/to/data -o /path/to/output"
echo ""
echo "To run with default directories:"
echo "  sbatch $0"
echo ""
echo "Default data directory: $default_data_dir"
echo "Default output directory: $default_output_base"
