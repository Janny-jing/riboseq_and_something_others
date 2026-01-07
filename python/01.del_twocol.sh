#!/bin/bash

# 遍历当前目录下的所有 .txt 文件
for file in *.txt; do
    # 检查是否为文件，防止意外操作非文件项
    if [ -f "$file" ]; then
        # 创建一个新的输出文件名，例如在原文件名前加上 'modified_'
        output="modified_$file"
        
        # 使用 awk 删除第二列并保存到新的文件中
        awk -F'\t' 'BEGIN{OFS="\t"} { $2=""; gsub(/^[\t]+|[\t]+$/, ""); print }' "$file" > "$output"
        
        # 如果你需要覆盖原来的文件，可以取消下面这行的注释
        mv "$output" "$file"
        
        echo "Processed $file and saved as $output"
    fi
done
