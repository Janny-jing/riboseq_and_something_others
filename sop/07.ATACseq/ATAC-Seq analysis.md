ATAC-Seq analysis 参数：

1. 质控：

Trim Galore （ 0.6.10）using cutadapt（with parameters：--illumina  -paired）

其中corresponding to cutadapt parameters: -e 0.1 -q 20 -O 1是默认参数

```
trim_galore --illumina -paired ${n}_1.fq.gz ${n}_2.fq.gz
```

2. 比对：

Bowtie2 （2.5.4）（with parameters：-p 8 --end-to-end --very-sensitive --dovetail  --no-mixed --no-discordant --fr -k 20 -X 1000）

PCR duplicates in the resulting BAM files were marked and removed with SAMtools , using samtools rmdup. MarkDuplicates and filtered with SAMtools（1.20）（with parameters：-q 30 -cat-F 256 -F 4）. retaining only properly-paired primary alignments, removing low-quality alignments and
those that fall within blacklisted regions（[ENCFF356LFX – ENCODE (encodeproject.org)](https://www.encodeproject.org/files/ENCFF356LFX/)） /（[Blacklist/lists at master · Boyle-Lab/Blacklist · GitHub](https://github.com/Boyle-Lab/Blacklist/tree/master/lists)）and the mitochondrial chromosome

```
#! /bin/sh

index="/home/jj2024/00.db/01.hg38/01.hg38_dna_bow2idx/hg38_dna"

for i in *_1_val_1.fq.gz;do
    j=${i/_1_val_1.fq.gz/}
        bowtie2 -p 8 --end-to-end --very-sensitive --dovetail  --no-mixed --no-discordant --fr -k 20 -X 1000  -x ${index} -1 ${j}_1_val_1.fq.gz -2 ${j}_2_val_2.fq.gz  | \
        samtools view -h -q 30 -f 2 -F 256 -F 4 -b - | \
        samtools sort -@ 8 -n - | \
        samtools rmdup - ${j}.bam >err 2>log
        samtools index ${j}.bam
done
```

详解：

```bash
--end-to-end 比对是将整个read和参考序列进行比对. 该模式--ma的值为0. 该模式为默认模式。
--very-sensitive Same as: -D 20 -R 3 -N 0 -L 20 -i S,1,0.50
--dovetail: 当配对读段相互延伸时视为一致。
--no-mixed       默认设置下, 一对reads不能成对比对到参考序列上, 则单独对每个read进行比对. 该选项则阻止此行为.
--no-discordant  默认设置下, 一对reads不能和谐比对(concordant alignment,即满足-I, -X, --fr/--rf/--ff的条件)到参考序列上, 则搜寻其不和谐比对(disconcordant alignment, 即两条reads都能独一无二地比对到参考序列上, 但是不满足-I,-X,--fr/--rf/--ff的条件). 该选项阻止此行为.--fr/--rf/--ff: -1, -2 配对读段对齐的方向。
-k <int>: 报告每个读段最多 <int> 个对齐结果
-X/--maxins <int>: 最大片段长度（默认为 500）。
```

3. call peak：

Genrich（0.6.1）(-E GRCh38_blacklist.bed -a 200 -m 30 -q 0.05 -j -r -e chrM,chrY -v）

详解：`-a <float>`: 峰值的最小 AUC 阈值（默认为 200.0）。

​            `-m <int>`: 保留对齐的最小 MAPQ 分数（默认为 0）。

​            `-q <float>`: 最大 q 值阈值（FDR-adjusted p-value; 默认为 1）。

​            `-j`: 使用 ATAC-seq 模式（默认为关闭）。

​            `-r`: 移除 PCR 重复。

​            `-e <arg>`: 排除的染色体列表，用逗号分隔。

​            `-v`: 打印状态更新和计数到标准错误输出（stderr）。

```
Genrich -t A549_control_r1.bam -o A549_control_r1.narrowpeak -E /home/jj2024/00.db/03.blacklist/GRCh38_blacklist.bed -a 200 -m 30 -q 0.05 -j -r -e chrM,chrY -v

##macs2
 macs2 callpeak -f BAMPE -t  ${n}.bam  -n ${n}  --outdir macs2  -g hs --keep-dup all -q 0.05 --fe-cutoff 2
```

