#!/bin/bash

# =============================================
# 从Rfam数据库提取人类非编码RNA序列的完整脚本
# 修复Python脚本的序列头匹配问题
# =============================================

set -e  # 遇到任何错误立即退出

# 配置参数
WORK_DIR="./rfam_human_ncrna"
RFAM_MYSQL_HOST="mysql-rfam-public.ebi.ac.uk"
RFAM_MYSQL_PORT="4497"
RFAM_MYSQL_USER="rfamro"
RFAM_DATABASE="Rfam"

# 检查依赖
check_dependencies() {
    echo "检查必要工具..."
    local missing=()
    
    ! command -v mysql &> /dev/null && missing+=("mysql-client")
    ! command -v bowtie2 &> /dev/null && missing+=("bowtie2")
    
    if [ ${#missing[@]} -gt 0 ]; then
        echo "错误: 缺少以下工具: ${missing[*]}"
        echo "请运行: conda install -c bioconda mysql-client bowtie2"
        exit 1
    fi
    echo "所有必要工具已安装"
}

# 检查MySQL连接
check_mysql_connection() {
    echo "检查MySQL连接..."
    if mysql -u "$RFAM_MYSQL_USER" -h "$RFAM_MYSQL_HOST" -P "$RFAM_MYSQL_PORT" -e "USE $RFAM_DATABASE;" 2>/dev/null; then
        echo "MySQL连接成功"
        return 0
    else
        echo "错误: 无法连接到Rfam MySQL数据库"
        echo "请检查网络连接或使用VPN"
        return 1
    fi
}

# 主程序
echo "=================================================="
echo "开始提取人类非编码RNA序列"
echo "工作目录: $WORK_DIR"
echo "=================================================="

# 检查依赖
check_dependencies

# 创建工作目录
mkdir -p "$WORK_DIR"
cd "$WORK_DIR"

# 步骤1: 检查Rfam序列文件
echo "步骤1/6: 检查Rfam序列文件..."
if [ ! -f "Rfam.fa" ]; then
    echo "错误: Rfam.fa 文件不存在"
    echo "请确保已在当前目录下载Rfam.fa文件"
    exit 1
fi

echo "Rfam.fa 文件存在，文件大小: $(du -h Rfam.fa | cut -f1)"

# 步骤2: 检查MySQL连接
echo "步骤2/6: 检查数据库连接..."
if ! check_mysql_connection; then
    echo "无法连接数据库，退出脚本"
    exit 1
fi

# 步骤3: 创建SQL查询文件
echo "步骤3/6: 创建SQL查询文件..."

# 创建提取所有人类ncRNA的SQL查询
cat > query_all_human_ncrna.sql << 'EOF'
SELECT CONCAT(fr.rfamseq_acc, '/', fr.seq_start, '-', fr.seq_end) as sequence_id
FROM full_region fr 
JOIN genseq gs ON gs.rfamseq_acc = fr.rfamseq_acc
WHERE fr.is_significant = 1
AND fr.type = 'full'
AND gs.upid = 'UP000005640';
EOF

# 创建提取人类snoRNA的SQL查询
cat > query_human_snorna.sql << 'EOF'
SELECT CONCAT(fr.rfamseq_acc, '/', fr.seq_start, '-', fr.seq_end) as sequence_id
FROM full_region fr 
JOIN genseq gs ON gs.rfamseq_acc = fr.rfamseq_acc
JOIN family f ON f.rfam_acc = fr.rfam_acc
WHERE fr.is_significant = 1
AND fr.type = 'full'
AND gs.upid = 'UP000005640'
AND f.type LIKE '%snoRNA%';
EOF

# 创建提取人类miRNA的SQL查询
cat > query_human_mirna.sql << 'EOF'
SELECT CONCAT(fr.rfamseq_acc, '/', fr.seq_start, '-', fr.seq_end) as sequence_id
FROM full_region fr 
JOIN genseq gs ON gs.rfamseq_acc = fr.rfamseq_acc
JOIN family f ON f.rfam_acc = fr.rfam_acc
WHERE fr.is_significant = 1
AND fr.type = 'full'
AND gs.upid = 'UP000005640'
AND f.type LIKE '%miRNA%';
EOF

# 创建提取人类rRNA的SQL查询
cat > query_human_rrna.sql << 'EOF'
SELECT CONCAT(fr.rfamseq_acc, '/', fr.seq_start, '-', fr.seq_end) as sequence_id
FROM full_region fr 
JOIN genseq gs ON gs.rfamseq_acc = fr.rfamseq_acc
JOIN family f ON f.rfam_acc = fr.rfam_acc
WHERE fr.is_significant = 1
AND fr.type = 'full'
AND gs.upid = 'UP000005640'
AND f.type LIKE '%rRNA%';
EOF

# 步骤4: 执行SQL查询获取序列ID
echo "步骤4/6: 执行SQL查询获取序列ID..."

execute_sql_query() {
    local query_file=$1
    local output_file=$2
    local description=$3
    
    echo "提取${description}序列ID..."
    if mysql -u "$RFAM_MYSQL_USER" -h "$RFAM_MYSQL_HOST" -P "$RFAM_MYSQL_PORT" \
          --skip-column-names --database "$RFAM_DATABASE" \
          < "$query_file" > "$output_file" 2>/dev/null; then
        
        count=$(wc -l < "$output_file" 2>/dev/null || echo 0)
        echo "  → 成功找到 ${count} 个${description}序列"
        
        if [ "$count" -gt 0 ]; then
            echo "  → 前3个序列ID示例:"
            head -3 "$output_file" | sed 's/^/      /'
        fi
        return 0
    else
        echo "  → 错误: 提取${description}失败"
        return 1
    fi
}

execute_sql_query "query_all_human_ncrna.sql" "all_human_ncrna_accessions.txt" "所有人类ncRNA"
execute_sql_query "query_human_snorna.sql" "human_snorna_accessions.txt" "人类snoRNA"
execute_sql_query "query_human_mirna.sql" "human_mirna_accessions.txt" "人类miRNA"
execute_sql_query "query_human_rrna.sql" "human_rrna_accessions.txt" "人类rRNA"

# 步骤5: 使用修复的Python脚本提取序列
echo "步骤5/6: 使用修复的Python脚本提取FASTA序列..."

# 创建修复的Python提取脚本
cat > extract_sequences_fixed.py << 'PYTHON_EOF'
#!/usr/bin/env python3
import sys
import re

def debug_sequence_matching(fasta_file, accession_file):
    """调试序列匹配情况"""
    print("=== 调试序列匹配 ===")
    
    # 读取accession文件
    with open(accession_file, 'r') as f:
        target_accessions = [line.strip() for line in f if line.strip()]
    
    print(f"需要匹配的accession数量: {len(target_accessions)}")
    
    # 测试前几个accession
    test_accessions = target_accessions[:3]
    for i, acc in enumerate(test_accessions):
        print(f"测试accession {i+1}: '{acc}'")
        
        # 在FASTA文件中搜索
        with open(fasta_file, 'r') as f:
            found = False
            for line in f:
                if line.startswith('>'):
                    header = line[1:].strip()  # 去掉'>'符号
                    if acc in header:
                        print(f"  找到匹配! 完整序列头: {header}")
                        found = True
                        break
            if not found:
                print(f"  未找到匹配")
    
    print("=== 调试结束 ===\n")

def extract_sequences_fixed(fasta_file, accession_file, output_file):
    # 首先进行调试
    debug_sequence_matching(fasta_file, accession_file)
    
    # 读取所有要提取的accession
    with open(accession_file, 'r') as f:
        target_accessions = set(line.strip() for line in f if line.strip())
    
    print(f"需要提取 {len(target_accessions)} 个序列")
    
    # 读取FASTA文件并提取序列
    with open(fasta_file, 'r') as fin, open(output_file, 'w') as fout:
        current_header = None
        current_sequence = []
        found_count = 0
        line_count = 0
        
        for line in fin:
            line_count += 1
            if line_count % 1000000 == 0:
                print(f"已处理 {line_count} 行，找到 {found_count} 个序列")
                
            if line.startswith('>'):
                # 处理前一个序列
                if current_header is not None:
                    header_content = current_header[1:].strip()  # 去掉'>'符号
                    
                    # 检查这个头是否匹配任何target accession
                    for acc in target_accessions:
                        if acc in header_content:
                            # 找到匹配，写入序列
                            fout.write(current_header)
                            fout.writelines(current_sequence)
                            found_count += 1
                            if found_count % 1000 == 0:
                                print(f"已找到 {found_count} 个序列")
                            break
                
                # 开始新序列
                current_header = line
                current_sequence = []
            else:
                current_sequence.append(line)
        
        # 处理最后一个序列
        if current_header is not None:
            header_content = current_header[1:].strip()
            for acc in target_accessions:
                if acc in header_content:
                    fout.write(current_header)
                    fout.writelines(current_sequence)
                    found_count += 1
                    break
        
        print(f"处理完成: 总共扫描了 {line_count} 行")
        print(f"总共找到 {found_count} 个匹配序列")
        
        # 输出一些统计信息
        if found_count > 0:
            # 重新读取输出文件显示一些匹配示例
            with open(output_file, 'r') as f:
                headers_found = []
                current_header = None
                for line in f:
                    if line.startswith('>'):
                        current_header = line.strip()
                        headers_found.append(current_header)
                    if len(headers_found) >= 3:
                        break
            
            print("匹配的序列头示例:")
            for i, header in enumerate(headers_found[:3]):
                print(f"  {i+1}. {header}")
        else:
            print("错误: 没有找到任何匹配的序列")
            print("可能的原因:")
            print("  1. accession格式不匹配")
            print("  2. FASTA文件格式问题") 
            print("  3. 文件编码问题")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("用法: python extract_sequences_fixed.py fasta_file accession_file output_file")
        sys.exit(1)
    
    extract_sequences_fixed(sys.argv[1], sys.argv[2], sys.argv[3])
PYTHON_EOF

# 提取序列的函数
extract_sequences_fixed() {
    local accession_file=$1
    local output_file=$2
    local description=$3
    
    if [ ! -s "$accession_file" ]; then
        echo "警告: $accession_file 为空，跳过${description}提取"
        return 1
    fi
    
    echo "提取${description}序列..."
    echo "  → 使用修复的Python脚本..."
    
    # 首先手动测试一个accession
    test_acc=$(head -1 "$accession_file")
    echo "  → 测试第一个accession: '$test_acc'"
    matches=$(grep -c ">$test_acc" Rfam.fa 2>/dev/null || echo 0)
    echo "  → 在Rfam.fa中找到 $matches 个完全匹配"
    
    if [ "$matches" -eq 0 ]; then
        echo "  → 警告: 完全匹配失败，将使用包含匹配"
    fi
    
    python3 extract_sequences_fixed.py Rfam.fa "$accession_file" "$output_file"
    
    local count=$(grep -c "^>" "$output_file" 2>/dev/null || echo 0)
    if [ "$count" -gt 0 ]; then
        echo "  → ${description}提取完成: $output_file (${count} 条序列)"
    else
        echo "  → 错误: 未成功提取任何${description}序列"
        # 尝试最后的备用方法
        extract_sequences_final_backup "$accession_file" "$output_file" "$description"
    fi
}

# 最后的备用方法：使用grep和sed
extract_sequences_final_backup() {
    local accession_file=$1
    local output_file=$2
    local description=$3
    
    echo "  → 使用最终备用方法（grep+sed）..."
    
    # 为每个accession创建搜索模式
    > "$output_file"
    total=$(wc -l < "$accession_file")
    processed=0
    found=0
    
    while IFS= read -r acc; do
        if [ -n "$acc" ]; then
            # 使用grep找到包含该accession的头行
            header_line=$(grep -n ">$acc" Rfam.fa | head -1 | cut -d: -f1)
            if [ -n "$header_line" ]; then
                # 使用sed提取从这个头行开始到下一个头行或文件结束的序列
                next_header_line=$(grep -n "^>" Rfam.fa | awk -F: -v line="$header_line" '$1 > line {print $1; exit}')
                
                if [ -n "$next_header_line" ]; then
                    sed -n "${header_line},$((next_header_line - 1))p" Rfam.fa >> "$output_file"
                else
                    sed -n "${header_line},\$p" Rfam.fa >> "$output_file"
                fi
                ((found++))
            fi
            ((processed++))
            
            if [ $((processed % 1000)) -eq 0 ]; then
                echo "   进度: $processed/$total, 找到: $found"
            fi
        fi
    done < "$accession_file"
    
    echo "  → 最终方法完成: 处理 $processed, 找到 $found"
}

# 提取各种ncRNA类型
extract_sequences_fixed "all_human_ncrna_accessions.txt" "Homo_sapiens_all_ncRNA.fa" "所有人类ncRNA"
extract_sequences_fixed "human_snorna_accessions.txt" "Homo_sapiens_snoRNA.fa" "人类snoRNA"
extract_sequences_fixed "human_mirna_accessions.txt" "Homo_sapiens_miRNA.fa" "人类miRNA"
extract_sequences_fixed "human_rrna_accessions.txt" "Homo_sapiens_rRNA.fa" "人类rRNA"

# 清理临时文件
rm -f extract_sequences_fixed.py

# 步骤6: 为Bowtie2构建索引
echo "步骤6/6: 构建Bowtie2索引..."

if [ -f "Homo_sapiens_all_ncRNA.fa" ] && [ -s "Homo_sapiens_all_ncRNA.fa" ]; then
    sequence_count=$(grep -c "^>" "Homo_sapiens_all_ncRNA.fa")
    file_size=$(du -h "Homo_sapiens_all_ncRNA.fa" | cut -f1)
    echo "为所有人类ncRNA序列构建Bowtie2索引..."
    echo "  → 输入文件: $sequence_count 条序列, 大小: $file_size"
    
    if [ "$sequence_count" -gt 0 ]; then
        if bowtie2-build "Homo_sapiens_all_ncRNA.fa" "Homo_sapiens_all_ncRNA_index" 2>&1 | tee bowtie2_build.log; then
            echo "Bowtie2索引构建成功!"
            echo "生成的索引文件:"
            ls -la Homo_sapiens_all_ncRNA_index*.bt2 | head -5
        else
            echo "错误: Bowtie2索引构建失败"
            echo "请查看日志文件: bowtie2_build.log"
        fi
    else
        echo "错误: Homo_sapiens_all_ncRNA.fa 中没有序列"
    fi
else
    echo "警告: Homo_sapiens_all_ncRNA.fa 不存在或为空，跳过索引构建"
fi

# 生成最终报告
echo "=================================================="
echo "提取完成！最终结果:"
echo "=================================================="

for file in Homo_sapiens_*.fa; do
    if [ -f "$file" ]; then
        count=$(grep -c "^>" "$file" 2>/dev/null || echo 0)
        size=$(du -h "$file" 2>/dev/null | cut -f1 || echo "0")
        echo "  $(basename "$file"): $count 条序列, 大小: $size"
    fi
done

echo ""
echo "使用说明:"
if [ -f "Homo_sapiens_all_ncRNA_index.1.bt2" ]; then
    echo "Bowtie2过滤命令:"
    echo "  bowtie2 -x Homo_sapiens_all_ncRNA_index -1 reads_1.fq -2 reads_2.fq --un-conc-gz filtered_%.fq.gz"
else
    echo "Bowtie2索引未成功构建，请检查FASTA文件"
fi
echo "=================================================="
