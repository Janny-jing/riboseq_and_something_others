clipper——microRNA分析

1.microRNA数据库数据下载：https://www.mirbase.org/download/[hairpin.fa](https://www.mirbase.org/download/hairpin.fa)

2.由于这个是包含很多物种的文件，我们需要根据需求筛选出自己物种的fa文件（这里我们需要hsa），因此根据我们的需求写了一个python脚本

```
import argparse
import re
import sys

def filter_sequences(input_file, output_file, keyword):
    # 简化正则表达式，直接查找关键词
    pattern = re.compile(re.escape(keyword))

    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        current_header = ''
        current_sequence = []
        
        for line in infile:
            line = line.strip()
            if line.startswith('>'):
                # 如果当前序列不为空，则检查是否包含关键词
                if current_header and keyword in current_header:
                    print(f"Matched header: {current_header}", file=sys.stderr)
                    outfile.write(current_header + '\n' + ''.join(current_sequence) + '\n')
                
                current_header = line
                current_sequence = []
                print(f"Processing header: {current_header}", file=sys.stderr)
            else:
                current_sequence.append(line)
        
        # 检查最后一个序列
        if current_header and keyword in current_header:
            print(f"Matched header: {current_header}", file=sys.stderr)
            outfile.write(current_header + '\n' + ''.join(current_sequence) + '\n')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Filter FASTA sequences based on a keyword in the header.')
    
    parser.add_argument('-i', '--input', required=True, help='Input FASTA file path')
    parser.add_argument('-o', '--output', required=True, help='Output FASTA file path')
    parser.add_argument('-k', '--keyword', default='has', help='Keyword to search for in the sequence headers (default: has)')
    
    args = parser.parse_args()

    filter_sequences(args.input, args.output, args.keyword)
```

3.取出序列中含有hsa的行及其对应的序列

```
python fa_keyword.py -i hairpin.fa -o microRNA_sequences.fa -k hsa
```

4.对序列进行bowtie2建立索引

```
bowtie2-build -f <name>.fa <name>
```

5.对序列进行比对(前期以对序列进行了质控)

```
#! /bin/bash

#SBATCH --job-name=01
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --time=2-00:00:00

index="/home/gongweikang/zhangxue_file/db/08.microRNA/microRNA_sequences"

for i in *_umi_trim.fq;do
    n=${i/_umi_trim.fq/}
        bowtie2 -p 8 --local --very-sensitive-local --no-unal --no-mixed --no-discordant -x ${index} -U ${n}_umi_trim.fq | \
        samtools view -b -F 4 - | \
        samtools sort -@ 8 - | \
        samtools rmdup -s -  ${n}.bam >err 2>log
        samtools index ${n}.bam
        samtools idxstats ${n}.bam > ${n}.txt
done
```

