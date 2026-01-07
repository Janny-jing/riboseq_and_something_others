#!/bin/bash
# 保存为 filter_bams.sh，然后执行 chmod +x filter_bams.sh

# 设置工作目录
cd /home/jj2024/22.riboseq_three/fram_graph

# 创建标准染色体数组
standard_chr=($(cat standard_chromosomes.txt))

# 查找所有BAM文件（排除索引文件）
for bam_file in *.bam; do
    # 跳过索引文件
    if [[ "$bam_file" == *.bai ]] || [[ "$bam_file" == *.csi ]]; then
        continue
    fi
    
    echo "正在处理: $bam_file"
    
    # 生成输出文件名
    base_name="${bam_file%.bam}"
    output_bam="${base_name}_standardChr.bam"
    
    # 统计原始reads数
    original_count=$(samtools view -c "$bam_file")
    
    # 执行过滤
    samtools view -b -h "$bam_file" "${standard_chr[@]}" > "$output_bam"
    
    # 建立索引
    samtools index "$output_bam"
    
    # 统计过滤后reads数
    filtered_count=$(samtools view -c "$output_bam")
    
    # 计算保留比例
    retention=$(echo "scale=2; $filtered_count * 100 / $original_count" | bc)
    
    # 输出统计信息
    echo "  原始reads数: $original_count"
    echo "  过滤后reads数: $filtered_count"
    echo "  保留比例: ${retention}%"
    echo ""
done

echo "批量过滤完成！"
