#!/bin/bash
# 脚本名称：riboseq_frame_complete.sh
# 完整的Ribo-seq frame分析，基于CDS起始位置

set -e  # 出错时退出

if [ $# -lt 2 ]; then
    echo "用法: $0 <gtf_file> <bam_file1> [bam_file2 ...]"
    echo "示例: $0 gencode.gtf *.bam"
    exit 1
fi

GTF_FILE="$1"
shift
BAM_FILES=("$@")

echo "=================================================="
echo "      RIBO-SEQ FRAME 分析工具"
echo "=================================================="
echo "GTF文件: $(basename "$GTF_FILE")"
echo "样本数: ${#BAM_FILES[@]}"
echo "开始时间: $(date)"
echo "=================================================="

# 检查必需工具
command -v samtools >/dev/null 2>&1 || { echo "错误: 需要samtools"; exit 1; }
command -v bedtools >/dev/null 2>&1 || { echo "错误: 需要bedtools"; exit 1; }

# 创建输出目录
OUTPUT_DIR="riboseq_frame_results_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$OUTPUT_DIR"
echo "输出目录: $OUTPUT_DIR"

# 创建临时目录
TEMP_DIR="$OUTPUT_DIR/temp"
mkdir -p "$TEMP_DIR"

# 输出文件
SUMMARY_CSV="$OUTPUT_DIR/summary.csv"
DETAILED_CSV="$OUTPUT_DIR/detailed_results.csv"
LOG_FILE="$OUTPUT_DIR/analysis.log"

# 写入日志
exec > >(tee -a "$LOG_FILE") 2>&1

# ============================================================================
# 第一步：从GTF提取蛋白质编码基因的CDS起始信息
# ============================================================================
echo ""
echo "步骤1: 提取蛋白质编码基因的CDS起始信息"
echo "----------------------------------------"

# 提取所有蛋白质编码基因的起始密码子
echo "提取起始密码子位置..."
awk 'BEGIN {OFS="\t"; FS="\t"}
# 跳过注释行
/^#/ {next}

# 提取基因信息（只保留蛋白质编码基因）
$3 == "gene" {
    gene_id = ""
    gene_type = ""
    gene_name = ""
    
    # 解析属性字段
    split($9, attrs, ";")
    for(i=1; i<=length(attrs); i++) {
        gsub(/^[[:space:]]+|[[:space:]]+$/, "", attrs[i])
        if(attrs[i] ~ /^gene_id/) {
            split(attrs[i], parts, "\"")
            gene_id = parts[2]
        }
        if(attrs[i] ~ /^gene_type/) {
            split(attrs[i], parts, "\"")
            gene_type = parts[2]
        }
        if(attrs[i] ~ /^gene_name/) {
            split(attrs[i], parts, "\"")
            gene_name = parts[2]
        }
    }
    
    # 只保存蛋白质编码基因
    if(gene_type == "protein_coding" && gene_id != "") {
        genes[gene_id] = $1 "\t" $4 "\t" $5 "\t" $7 "\t" gene_name
    }
}

# 提取起始密码子
$3 == "start_codon" {
    gene_id = ""
    transcript_id = ""
    
    split($9, attrs, ";")
    for(i=1; i<=length(attrs); i++) {
        gsub(/^[[:space:]]+|[[:space:]]+$/, "", attrs[i])
        if(attrs[i] ~ /^gene_id/) {
            split(attrs[i], parts, "\"")
            gene_id = parts[2]
        }
        if(attrs[i] ~ /^transcript_id/) {
            split(attrs[i], parts, "\"")
            transcript_id = parts[2]
        }
    }
    
    # 如果是蛋白质编码基因且有转录本ID
    if(gene_id in genes && transcript_id != "") {
        split(genes[gene_id], gene_info, "\t")
        chr = gene_info[1]
        gene_start = gene_info[2]
        gene_end = gene_info[3]
        strand = gene_info[4]
        gene_name = gene_info[5]
        
        start_codon_pos = $4
        
        # 对于每个转录本，我们取第一个起始密码子
        if(!(transcript_id in processed)) {
            # 计算相对于基因起始的offset（如果需要）
            if(strand == "+") {
                offset_from_gene_start = start_codon_pos - gene_start
            } else {
                offset_from_gene_start = gene_end - start_codon_pos
            }
            
            # 输出: chr, start_codon, end_codon, gene_name, strand
            # 起始密码子长度通常是3nt
            end_codon_pos = start_codon_pos + 2
            
            print chr "\t" start_codon_pos-1 "\t" end_codon_pos "\t" \
                  gene_name ":" transcript_id "\t0\t" strand
            
            processed[transcript_id] = 1
            count++
        }
    }
}

END {
    print "提取了 " count " 个起始密码子" > "/dev/stderr"
}' "$GTF_FILE" > "$TEMP_DIR/start_codons.bed"

START_COUNT=$(wc -l < "$TEMP_DIR/start_codons.bed")
echo "成功提取 $START_COUNT 个起始密码子"

if [ "$START_COUNT" -eq 0 ]; then
    echo "警告: 未找到起始密码子，尝试提取CDS区域作为替代..."
    # 提取CDS区域
    awk 'BEGIN {OFS="\t"} $3 == "CDS" {
        # 提取第一个CDS（通常是起始位置）
        gene_name = "unknown"
        split($9, attrs, ";")
        for(i=1; i<=length(attrs); i++) {
            gsub(/^[[:space:]]+|[[:space:]]+$/, "", attrs[i])
            if(attrs[i] ~ /^gene_name/) {
                split(attrs[i], parts, "\"")
                gene_name = parts[2]
                break
            }
        }
        print $1 "\t" $4-1 "\t" $5 "\t" gene_name "\t0\t" $7
    }' "$GTF_FILE" | head -10000 > "$TEMP_DIR/start_codons.bed"
    
    START_COUNT=$(wc -l < "$TEMP_DIR/start_codons.bed")
    echo "使用CDS区域作为替代: $START_COUNT 个区域"
fi

# ============================================================================
# 第二步：处理每个BAM文件
# ============================================================================
echo ""
echo "步骤2: 处理BAM文件"
echo "----------------------------------------"

# 输出文件头
echo "sample,total_reads,cds_reads,frame0,frame1,frame2,frame0_pct,frame1_pct,frame2_pct,enrichment,quality" > "$SUMMARY_CSV"
echo "sample,gene,strand,read_pos,cds_start,offset,frame" > "$DETAILED_CSV"

for BAM in "${BAM_FILES[@]}"; do
    if [ ! -f "$BAM" ]; then
        echo "警告: 跳过不存在的文件 $BAM"
        continue
    fi
    
    SAMPLE_NAME=$(basename "$BAM" .bam | sed 's/\.bam//; s/_filter//; s/_sorted//; s/_aligned//')
    echo ""
    echo "处理样本: $SAMPLE_NAME"
    echo "--------------------------------"
    
    # 检查BAM文件
    if ! samtools quickcheck "$BAM"; then
        echo "错误: BAM文件损坏或格式不正确"
        continue
    fi
    
    # 统计总reads
    TOTAL_READS=$(samtools view -c -F 0x4 "$BAM")
    echo "总比对reads: $(printf "%'d" $TOTAL_READS)"
    
    # 提取reads的5'端位置（P-site）
    echo "提取reads 5'端位置（P-site）..."
    SAMPLE_TEMP_DIR="$TEMP_DIR/$SAMPLE_NAME"
    mkdir -p "$SAMPLE_TEMP_DIR"
    
    # 方法1：使用固定偏移（对于Ribo-seq，通常5'端+12nt或15nt是P-site）
    # 这里简化处理：使用5'端位置
    OFFSET_FOR_PSITE=12  # 可根据实际数据调整
    
    samtools view -F 0x4 "$BAM" | \
    awk -v offset="$OFFSET_FOR_PSITE" -v out="$SAMPLE_TEMP_DIR/reads_psite.bed" '
    BEGIN {
        OFS = "\t"
        count = 0
    }
    {
        chr = $3
        pos = $4  # 1-based起始位置
        flag = $2
        cigar = $6
        seq_len = length($10)
        
        # 跳过非常规染色体
        if (!(chr ~ /^chr[0-9XYM]/ || chr ~ /^[0-9XYM]/)) next
        
        # 确定链方向
        if (flag >= 16) {
            strand = "-"
            # 反向read：5端在3端
            # 计算比对长度（简化：取CIGAR中的M）
            aligned_len = seq_len
            if (cigar ~ /^[0-9]+M/) {
                match(cigar, /^[0-9]+/)
                aligned_len = substr(cigar, RSTART, RLENGTH)
            }
            five_prime = pos + aligned_len - 1  # 5端位置
            # P-site位置：5端 - offset（因为反向）
            psite_pos = five_prime - offset
        } else {
            strand = "+"
            five_prime = pos  # 5端位置
            # P-site位置：5端 + offset
            psite_pos = five_prime + offset
        }
        
        # 确保P-site位置有效
        if (psite_pos < 1) next
        
        # 转换为0-based BED格式
        bed_start = psite_pos - 1
        
        print chr, bed_start, bed_start+1, "read_" NR, 0, strand > out
        count++
        
        if (count % 1000000 == 0) {
            printf "  已处理 %dM reads\n", count/1000000 > "/dev/stderr"
        }
    }
    END {
        print "提取了 " count " 个reads的P-site位置" > "/dev/stderr"
    }'
    
    # 过滤起始密码子区域的reads
    echo "过滤起始密码子区域的reads..."
    bedtools intersect -a "$SAMPLE_TEMP_DIR/reads_psite.bed" \
                       -b "$TEMP_DIR/start_codons.bed" \
                       -s -wa -wb > "$SAMPLE_TEMP_DIR/reads_in_start.bed" 2>/dev/null
    
    READS_IN_START=$(wc -l < "$SAMPLE_TEMP_DIR/reads_in_start.bed" 2>/dev/null || echo "0")
    echo "起始密码子区域reads: $(printf "%'d" $READS_IN_START)"
    
    if [ "$READS_IN_START" -lt 100 ]; then
        echo "警告: 起始密码子区域reads太少，可能数据质量有问题"
        # 使用CDS区域作为备选
        echo "尝试使用CDS区域..."
        awk 'BEGIN {OFS="\t"} $3 == "CDS" {
            print $1 "\t" $4-1 "\t" $5 "\tCDS\t0\t" $7
        }' "$GTF_FILE" | head -50000 > "$TEMP_DIR/cds_regions.bed"
        
        bedtools intersect -a "$SAMPLE_TEMP_DIR/reads_psite.bed" \
                           -b "$TEMP_DIR/cds_regions.bed" \
                           -s -wa -wb > "$SAMPLE_TEMP_DIR/reads_in_cds.bed" 2>/dev/null
        
        READS_IN_CDS=$(wc -l < "$SAMPLE_TEMP_DIR/reads_in_cds.bed" 2>/dev/null || echo "0")
        echo "CDS区域reads: $(printf "%'d" $READS_IN_CDS)"
        
        ANALYSIS_FILE="$SAMPLE_TEMP_DIR/reads_in_cds.bed"
        READS_FOR_ANALYSIS=$READS_IN_CDS
    else
        ANALYSIS_FILE="$SAMPLE_TEMP_DIR/reads_in_start.bed"
        READS_FOR_ANALYSIS=$READS_IN_START
    fi
    
    # 分析frame分布
    echo "分析frame分布..."
    if [ "$READS_FOR_ANALYSIS" -gt 100 ]; then
        awk -v sample="$SAMPLE_NAME" '
        BEGIN {
            FS = OFS = "\t"
            f0 = f1 = f2 = 0
            total = 0
            
            # 用于详细输出
            detail_out = "'"$SAMPLE_TEMP_DIR"'/frame_detail.txt"
            print "gene,strand,read_pos,cds_start,offset,frame" > detail_out
        }
        {
            # read信息: $1-$6
            # start_codon/CDS信息: $7-$12
            read_chr = $1
            read_start = $2  # 0-based P-site位置
            read_strand = $6
            
            cds_chr = $7
            cds_start = $8   # 0-based CDS起始
            cds_end = $9
            cds_name = $10
            cds_strand = $12
            
            # 确保染色体和链匹配
            if (read_chr != cds_chr || read_strand != cds_strand) next
            
            # 计算相对于CDS起始的offset
            if (cds_strand == "+") {
                # 正向链：从CDS起始开始计算
                offset = read_start - cds_start
            } else {
                # 反向链：从CDS结束开始计算（反向）
                # 注意：对于反向链，cds_start实际上是BED格式的起始（较小值）
                # 但CDS的真实起始在右侧
                cds_real_start = cds_end - 3  # 起始密码子长度3nt
                offset = cds_real_start - read_start
            }
            
            # 只考虑合理的offset（0-1000，避免过大偏移）
            if (offset >= 0 && offset < 1000) {
                frame = offset % 3
                
                # 计数
                if (frame == 0) f0++
                else if (frame == 1) f1++
                else if (frame == 2) f2++
                
                total++
                
                # 保存详细数据
                print cds_name, cds_strand, read_start+1, cds_start+1, offset, frame >> detail_out
            }
        }
        END {
            if (total > 0) {
                # 计算比例
                p0 = f0 / total
                p1 = f1 / total
                p2 = f2 / total
                
                # 计算Frame 0富集度
                enrichment = p0 / 0.333333
                
                # 质量评估
                if (p0 > 0.50) quality = "Excellent"
                else if (p0 > 0.45) quality = "Good"
                else if (p0 > 0.40) quality = "Fair"
                else if (p0 > 0.35) quality = "Marginal"
                else quality = "Poor"
                
                # 输出摘要
                printf "%s,%d,%d,%d,%d,%d,%.3f,%.3f,%.3f,%.2f,%s\n",
                    sample, '"$TOTAL_READS"', total, f0, f1, f2,
                    p0, p1, p2, enrichment, quality
                    
                # 打印结果
                printf "\n分析结果:\n"
                printf "  有效reads数: %d\n", total
                printf "  Frame分布:\n"
                printf "    Frame 0: %d (%.1f%%)\n", f0, p0*100
                printf "    Frame 1: %d (%.1f%%)\n", f1, p1*100
                printf "    Frame 2: %d (%.1f%%)\n", f2, p2*100
                printf "  Frame 0富集度: %.2f倍\n", enrichment
                printf "  质量评估: %s\n", quality
            } else {
                printf "%s,%d,0,0,0,0,0.000,0.000,0.000,0.00,No_data\n",
                    sample, '"$TOTAL_READS"'
                print "  无有效数据"
            }
        }' "$ANALYSIS_FILE" >> "$SUMMARY_CSV"
        
        # 合并详细数据
        if [ -f "$SAMPLE_TEMP_DIR/frame_detail.txt" ]; then
            tail -n +2 "$SAMPLE_TEMP_DIR/frame_detail.txt" | \
            awk -v sample="$SAMPLE_NAME" '{print sample "," $0}' >> "$DETAILED_CSV"
        fi
        
        # 为每个样本生成报告
        cat > "$OUTPUT_DIR/${SAMPLE_NAME}_report.txt" << EOF
========================================
Ribo-seq Frame分析报告
========================================
样本: $SAMPLE_NAME
分析时间: $(date)
总reads数: $(printf "%'d" $TOTAL_READS)

关键指标:
---------
有效分析reads数: $(printf "%'d" $READS_FOR_ANALYSIS)

Frame分布:
---------
$(awk -F',' -v sample="$SAMPLE_NAME" '$1==sample {
    printf "Frame 0: %d (%.1f%%)\n", $4, $7*100
    printf "Frame 1: %d (%.1f%%)\n", $5, $8*100
    printf "Frame 2: %d (%.1f%%)\n", $6, $9*100
}' "$SUMMARY_CSV")

质量评估:
---------
$(awk -F',' -v sample="$SAMPLE_NAME" '$1==sample {
    printf "Frame 0富集度: %.2f倍\n", $10
    printf "评估结果: %s\n", $11
    print ""
    if($11 == "Excellent" || $11 == "Good") {
        print "✓ 数据质量良好，符合Ribo-seq特征"
        print "✓ Frame 0偏好明显，表明有效的翻译暂停"
    } else if($11 == "Fair" || $11 == "Marginal") {
        print "○ 数据质量一般，Frame 0偏好较弱"
        print "○ 可能需要优化实验条件或数据分析参数"
    } else {
        print "✗ 数据可能存在问题"
        print "✗ Frame分布接近随机，可能不是Ribo-seq数据"
    }
}' "$SUMMARY_CSV")

建议:
------
1. 高质量Ribo-seq数据通常Frame 0 > 45%
2. 如果Frame 0 < 40%，建议检查：
   - 样本处理（环己酰亚胺处理时间）
   - 核糖体消化条件
   - 文库构建质量
3. 可考虑调整P-site偏移参数重新分析

EOF
    else
        echo "有效reads数不足 ($READS_FOR_ANALYSIS)，跳过详细分析"
        echo "$SAMPLE_NAME,$TOTAL_READS,$READS_FOR_ANALYSIS,0,0,0,0.000,0.000,0.000,0.00,Insufficient_data" >> "$SUMMARY_CSV"
    fi
    
    # 生成可视化数据
    echo "生成可视化数据..."
    awk -F',' -v sample="$SAMPLE_NAME" '$1==sample {
        printf "Frame\tCount\tPercentage\n"
        printf "0\t%d\t%.1f\n", $4, $7*100
        printf "1\t%d\t%.1f\n", $5, $8*100
        printf "2\t%d\t%.1f\n", $6, $9*100
    }' "$SUMMARY_CSV" > "$OUTPUT_DIR/${SAMPLE_NAME}_plot_data.tsv"
done

# ============================================================================
# 第三步：生成汇总报告
# ============================================================================
echo ""
echo "步骤3: 生成汇总报告"
echo "----------------------------------------"

# 汇总统计
cat > "$OUTPUT_DIR/summary_report.txt" << EOF
========================================
Ribo-seq Frame分析汇总报告
========================================
分析时间: $(date)
样本数量: ${#BAM_FILES[@]}
GTF文件: $(basename "$GTF_FILE")

样本结果汇总:
------------
$(awk -F',' 'NR==1 {printf "%-20s %12s %12s %8s %8s %8s %10s %10s\n", 
    "Sample", "Total_Reads", "CDS_Reads", "Frame0%", "Frame1%", "Frame2%", "Enrichment", "Quality"}
NR>1 {printf "%-20s %12d %12d %8.1f %8.1f %8.1f %10.2f %10s\n", 
    $1, $2, $3, $7*100, $8*100, $9*100, $10, $11}' "$SUMMARY_CSV")

质量统计:
---------
$(awk -F',' 'NR>1 {
    if($11 == "Excellent") excellent++
    else if($11 == "Good") good++
    else if($11 == "Fair") fair++
    else if($11 == "Marginal") marginal++
    else if($11 == "Poor") poor++
    total++
}
END {
    print "总样本数: " total
    print "Excellent (Frame0 > 50%): " excellent
    print "Good (Frame0 45-50%): " good
    print "Fair (Frame0 40-45%): " fair
    print "Marginal (Frame0 35-40%): " marginal
    print "Poor (Frame0 < 35%): " poor
}' "$SUMMARY_CSV")

Frame分布统计:
-------------
$(awk -F',' 'NR>1 {
    frame0_avg += $7
    frame1_avg += $8
    frame2_avg += $9
    count++
}
END {
    if(count > 0) {
        print "平均Frame分布:"
        printf "  Frame 0: %.1f%%\n", frame0_avg/count*100
        printf "  Frame 1: %.1f%%\n", frame1_avg/count*100
        printf "  Frame 2: %.1f%%\n", frame2_avg/count*100
    }
}' "$SUMMARY_CSV")

数据质量建议:
------------
$(awk -F',' 'NR>1 && $7 > 0.45 {good_samples++} 
END {
    if(good_samples > 0) {
        print "✓ " good_samples " 个样本显示良好的Ribo-seq特征"
    } else {
        print "⚠ 所有样本Frame 0比例均低于45%，可能需要:"
        print "  1. 检查实验条件（环己酰亚胺处理）"
        print "  2. 验证是否为真正的Ribo-seq数据"
        print "  3. 调整P-site偏移参数重新分析"
    }
}' "$SUMMARY_CSV")

生成的文件:
----------
1. summary.csv - 汇总结果
2. detailed_results.csv - 详细数据
3. <sample>_report.txt - 各样本详细报告
4. <sample>_plot_data.tsv - 绘图数据
5. analysis.log - 分析日志

EOF

# 显示简要结果
echo ""
echo "=================================================="
echo "分析完成！"
echo "=================================================="
echo "主要结果:"
echo "----------"
awk -F',' 'NR==1 {printf "%-15s %10s %10s %8s %8s %8s\n", 
    "Sample", "Total", "CDS", "Frame0%", "Frame1%", "Frame2%"}
NR>1 {printf "%-15s %10d %10d %8.1f %8.1f %8.1f\n", 
    $1, $2, $3, $7*100, $8*100, $9*100}' "$SUMMARY_CSV" | head -10

echo ""
echo "数据质量统计:"
awk -F',' 'NR>1 {
    if($11 == "Excellent" || $11 == "Good") good++
    else if($11 == "Fair" || $11 == "Marginal") medium++
    else poor++
} END {
    print "  良好: " good " 个样本"
    print "  一般: " medium " 个样本"
    print "  较差: " poor " 个样本"
}' "$SUMMARY_CSV"

echo ""
echo "所有结果已保存到: $OUTPUT_DIR"
echo "查看完整报告: cat $OUTPUT_DIR/summary_report.txt"

# 清理临时文件（可选）
# echo "清理临时文件..."
# rm -rf "$TEMP_DIR"
