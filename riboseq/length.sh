#!/bin/bash
# riboseq_length_distribution_fixed.sh - 修复版长度分布分析

# === 设置参数 ===
OUTPUT_DIR="./length_distribution"
mkdir -p ${OUTPUT_DIR}

echo "分析Ribo-seq reads长度分布..."

# 处理所有BAM文件
for bam_file in *.bam; do
    if [ ! -f "$bam_file" ]; then
        continue
    fi
    
    sample_name=$(basename "$bam_file" .bam)
    echo "处理样本: $sample_name"
    
    # 提取reads长度并统计
    samtools view "$bam_file" | \
    awk '{print length($10)}' | \
    sort | uniq -c | \
    sort -n -k2 > ${OUTPUT_DIR}/${sample_name}_length_stats.txt
    
    echo "✅ 长度统计: ${OUTPUT_DIR}/${sample_name}_length_stats.txt"
done

# 使用Python绘制柱状图（修复语法错误）
python3 << EOF
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import os

# 设置中文字体和样式
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False
sns.set_style("whitegrid")

# 读取所有样本的长度数据
data_list = []
for file in glob.glob("${OUTPUT_DIR}/*_length_stats.txt"):
    sample_name = os.path.basename(file).replace('_length_stats.txt', '')
    df = pd.read_csv(file, sep='\\s+', header=None, names=['count', 'length'])
    df['sample'] = sample_name
    df['percentage'] = df['count'] / df['count'].sum() * 100
    data_list.append(df)

# 合并数据
combined_data = pd.concat(data_list)

# 绘制多样本长度分布图
plt.figure(figsize=(12, 8))
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']  # 定义颜色
sample_names = combined_data['sample'].unique()

for i, sample in enumerate(sample_names):
    sample_data = combined_data[combined_data['sample'] == sample]
    color = colors[i % len(colors)]
    plt.bar(sample_data['length'], sample_data['percentage'], 
            alpha=0.7, label=sample, width=0.8, color=color)

plt.xlabel('Read Length (nt)')
plt.ylabel('Percentage (%)')
plt.title('Ribo-seq Reads Length Distribution')
plt.legend()
plt.xticks(range(15, 36))  # 常见的Ribo-seq长度范围
plt.tight_layout()
plt.savefig('${OUTPUT_DIR}/riboseq_length_distribution.png', dpi=300, bbox_inches='tight')
plt.close()

# 分别绘制每个样本的分布图
for sample in combined_data['sample'].unique():
    sample_data = combined_data[combined_data['sample'] == sample]
    
    plt.figure(figsize=(10, 6))
    bars = plt.bar(sample_data['length'], sample_data['percentage'], 
            color='steelblue', alpha=0.8, width=0.8)
    
    # 标记主要峰（修复f-string语法）
    main_peak_idx = sample_data['percentage'].idxmax()
    main_peak_length = sample_data.loc[main_peak_idx, 'length']
    main_peak_percentage = sample_data.loc[main_peak_idx, 'percentage']
    
    # 使用传统字符串格式化代替f-string
    peak_label = 'Main peak: {}nt ({:.1f}%)'.format(main_peak_length, main_peak_percentage)
    plt.axvline(x=main_peak_length, color='red', linestyle='--', 
                alpha=0.8, label=peak_label)
    
    plt.xlabel('Read Length (nt)')
    plt.ylabel('Percentage (%)')
    plt.title('Ribo-seq Reads Length Distribution - {}'.format(sample))
    plt.legend()
    plt.xticks(range(15, 36))
    plt.tight_layout()
    plt.savefig('${OUTPUT_DIR}/{}_length_distribution.png'.format(sample), dpi=300, bbox_inches='tight')
    plt.close()

print("所有图表已保存到 ${OUTPUT_DIR}/ 目录")

# 输出统计信息
print("\\nReads长度分布统计:")
for sample in combined_data['sample'].unique():
    sample_data = combined_data[combined_data['sample'] == sample]
    total_reads = sample_data['count'].sum()
    main_peak = sample_data.loc[sample_data['percentage'].idxmax()]
    print("{}:".format(sample))
    print("  总reads数: {}".format(total_reads))
    print("  主要峰: {}nt ({:.1f}%)".format(main_peak['length'], main_peak['percentage']))
    print("  长度范围: {}-{}nt".format(sample_data['length'].min(), sample_data['length'].max()))
EOF

echo "🎉 Ribo-seq reads长度分布分析完成！"
echo "📊 结果文件在: ${OUTPUT_DIR}/"
