#!/bin/bash

# 定义要处理的文件模式（例如 *.txt 或 *_no2ndcol.txt）
file_pattern="*.txt"  # 可以根据需要调整文件模式

# 使用 awk 处理并求和，针对每个文件
for file in $file_pattern; do
    if [[ -f "$file" ]]; then # 确保是文件
        output_file="${file%.txt}_sum.txt"
        awk '
        BEGIN {
            FS=OFS="\t";  # 设置输入和输出字段分隔符为制表符

            # 初始化所有可能的密码子到 all_codons 数组中
            nuc="TCAG";
            split(nuc, nucleotides, "");  # 将核苷酸字符串分割成数组
            for (i in nucleotides)
                for (j in nucleotides)
                    for (k in nucleotides)
                        all_codons[toupper(nucleotides[i] nucleotides[j] nucleotides[k])] = 0;
        }

        {
            codon = $1;
            count = $2;

            if (codon in all_codons) {
                all_codons[codon] += count;
            } else {
                print "Warning: Codon", codon, "is not a valid codon and will be ignored." > "/dev/stderr"
            }
        }

        END {
            # 输出所有密码子及其对应的总和
            for (codon in all_codons) {
                print codon, all_codons[codon];
            }
        }' "$file" | sort > "$output_file"
        echo "Processed $file and saved results to $output_file"
    fi
done
