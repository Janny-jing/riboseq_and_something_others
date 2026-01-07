#! /bin/sh

for i in *.txt;do
	n=${i/.txt/}
        awk '
        BEGIN {
               FS=OFS="\t";  # 设置输入和输出字段分隔符为制表符
        }

        {
             # 对于每一对列，以第一列为键累积计数
        for (pair=1; pair<=6; pair++) {  # 因为有6对列
             col1 = pair * 2 - 1;         # 第一、三、五...列
             col2 = pair * 2;             # 第二、四、六...列

             if (col2 <= NF) {
               # 使用数组来累积每个键的出现次数
                counts[pair, $col1] += 1;
                 }
          }
        }

         END {
        # 输出结果，为每组创建一个单独的输出块
           for (pair=1; pair<=6; pair++) {
              printf "Group %d:\n", pair;
               for (key in counts) {
                 split(key, key_parts, SUBSEP);  # 分割复合键
                    group = key_parts[1];
                   if (group == pair) {
                       printf "%s\t%d\n", key_parts[2], counts[key];
                    }
                }
              print "";  # 空行分隔不同组的结果
          }
        }' ${n}.txt > ${n}_sum.txt
done
