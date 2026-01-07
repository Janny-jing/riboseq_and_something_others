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
    # 对于每一对列，以第一列为键累积计数
    for (pair=1; pair<=6; pair++) {  # 因为有6对列
        col1 = pair * 2 - 1;         # 第一、三、五...列
        col2 = pair * 2;             # 第二、四、六...列

        if (col2 <= NF) {
            codon = $col1;
            if (codon in all_codons) {
                counts[pair, codon] += 1;
            }
        }
    }
}

END {
    # 输出结果，为每组创建一个单独的输出块
    for (pair=1; pair<=6; pair++) {
        printf "Group %d:\n", pair;

        # 遍历所有可能的密码子
        for (codon in all_codons) {
            # 获取当前密码子的计数，如果不存在则默认为0
            count = (pair SUBSEP codon) in counts ? counts[pair, codon] : 0;
            printf "%s\t%d\n", codon, count;
        }
        print "";  # 空行分隔不同组的结果
    }
}
