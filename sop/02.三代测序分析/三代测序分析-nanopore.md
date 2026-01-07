三代测序分析-nanopore

1.软件环境

```
# packages in environment at /home/jj2024/anaconda3/envs/nanopore:
#
# Name                    Version                   Build  Channel
_libgcc_mutex             0.1                 conda_forge    conda-forge
_openmp_mutex             4.5                       2_gnu    conda-forge
bwa                       0.7.18               he4a0461_1    bioconda
bzip2                     1.0.8                h4bc722e_7    conda-forge
c-ares                    1.33.0               ha66036c_0    conda-forge
ca-certificates           2024.7.4             hbcca054_0    conda-forge
dos2unix                  7.5.2                ha770c72_3    conda-forge
gettext                   0.22.5               he02047a_3    conda-forge
gettext-tools             0.22.5               he02047a_3    conda-forge
homer                     4.11            pl5262h4ac6f70_9    bioconda
htslib                    1.20                 h5efdd21_2    bioconda
keyutils                  1.6.1                h166bdaf_0    conda-forge
krb5                      1.21.3               h659f571_0    conda-forge
ld_impl_linux-64          2.40                 hf3520f5_7    conda-forge
libasprintf               0.22.5               he8f35ee_3    conda-forge
libasprintf-devel         0.22.5               he8f35ee_3    conda-forge
libcurl                   8.8.0                hca28451_1    conda-forge
libdeflate                1.21                 h4bc722e_0    conda-forge
libedit                   3.1.20191231         he28a2e2_2    conda-forge
libev                     4.33                 hd590300_2    conda-forge
libffi                    3.4.2                h7f98852_5    conda-forge
libgcc-ng                 14.1.0               h77fa898_0    conda-forge
libgettextpo              0.22.5               he02047a_3    conda-forge
libgettextpo-devel        0.22.5               he02047a_3    conda-forge
libgomp                   14.1.0               h77fa898_0    conda-forge
libidn2                   2.3.7                hd590300_0    conda-forge
libnghttp2                1.58.0               h47da74e_1    conda-forge
libnsl                    2.0.1                hd590300_0    conda-forge
libsqlite                 3.46.0               hde9e2c9_0    conda-forge
libssh2                   1.11.0               h0841786_0    conda-forge
libstdcxx-ng              14.1.0               hc0a3c3a_0    conda-forge
libunistring              0.9.10               h7f98852_0    conda-forge
libuuid                   2.38.1               h0b41bf4_0    conda-forge
libxcrypt                 4.4.36               hd590300_1    conda-forge
libzlib                   1.2.13               h4ab18f5_6    conda-forge
ncurses                   6.5                  h59595ed_0    conda-forge
openssl                   3.3.1                hb9d3cd8_3    conda-forge
perl                      5.32.1          7_hd590300_perl5    conda-forge
pip                       24.2               pyhd8ed1ab_0    conda-forge
python                    3.11.0          he550d4f_1_cpython    conda-forge
readline                  8.2                  h8228510_1    conda-forge
samtools                  1.20                 h50ea8bc_1    bioconda
seqkit                    2.8.2                h9ee0642_1    bioconda
setuptools                72.1.0             pyhd8ed1ab_0    conda-forge
subread                   2.0.6                he4a0461_2    bioconda
tk                        8.6.13          noxft_h4845f30_101    conda-forge
tzdata                    2024a                h0c530f3_0    conda-forge
unzip                     6.0                  h7f98852_3    conda-forge
wget                      1.21.4               hda4d442_0    conda-forge
wheel                     0.44.0             pyhd8ed1ab_0    conda-forge
xz                        5.2.6                h166bdaf_0    conda-forge
zlib                      1.2.13               h4ab18f5_6    conda-forge
zstd                      1.5.6                ha6fb4c9_0    conda-forge
```

2.对三代下机数据进行pod5文件或者fast5类型转换为fastq

```
##guppy碱基识别
./guppy/bin/guppy_basecaller -i /dir/fast5/ -s /dir/fastq/ -r --min_qscore 0 --flowcell FLO-MIN106 --kit SQK-RNA002 --num_callers 1 --cpu_threads_per_caller 8
##dorado碱基识别
dorado basecaller /home/jj2024/.dorado/models/rna004_130bps_hac@v5.1.0  pod5/ --emit-fastq > 20241107_sample5.fastq
```

3.RNADNA转换

```
seqkit seq ./rna.fq --rna2dna > dna.fastq
```

4.碱基质量筛选（去除q<6的）

```
seqkit seq -Q 6 -m 80 ${n}.fq | gzip -c >${n}.fq.gz
```

5.比对tRNA序列（正常到这一步，linux上的就结束了）

```
#! /bin/sh
ref="/home/jj2024/00.db/01.hg38/08.hg38_trna_nomt_3adapter/human_tRNA_nomt_3adapter.fa"

for i in *.fq;do
       n=${i/.fq/}
       bwa mem -W13 -k6 -xont2d -T20 $ref ${n}.fq |\
       samtools view -q1 -F4 -b - |\
       samtools sort - > ${n}.bam
       samtools index ${n}.bam
       samtools idxstats ${n}.bam >${n}_trna_stat_q1.txt
done
```



####接下来的步骤是用于探寻自建库数据，未比对上tRNA参考基因组的序列是否比对到了全基因组或者adapter上面

1.提取未必对序列数据

```
#! /bin/sh
ref="/home/jj2024/00.db/01.hg38/02.hg38_trna_mt/human_tRNA_mt.fa"

for i in *.fastq;do
       n=${i/.fastq/}
       bwa mem -W13 -k6 -xont2d -T20 $ref ${n}.fastq |\
       samtools view -f4 -b - |\
       samtools sort - > ${n}_mt.bam
       samtools index ${n}_mt.bam
done

for i in *bam;do
        n=${i/.bam/}
        samtools fastq -n ${n}.bam >${n}_unmapped.fq
done

```

2.比对全基因组

```
#! /bin/sh

ref="/home/jj2024/00.db/01.hg38/06.hg38_dna_bwaidx/hg38_dna.fa"

for i in *.fq;do
       n=${i/.fq/}
       bwa mem -W13 -k6 -xont2d -T20 $ref ${n}.fq |\
       samtools view -F4 -b - |\
       samtools sort - > ${n}.bam
       samtools index ${n}.bam
done

bedtools genomecov -ibam ${n}.bam -bg > coverage.bedgraph
```

3.比对adapter

```
#! /bin/sh
ref="/home/jj2024/00.db/01.hg38/07.adapter_3/adapted_3.fa"

for i in *.fq;do
       n=${i/.fq/}
       bwa mem -W13 -k6 -xont2d -T20 -M $ref ${n}.fq |\
       samtools view -q0 -F2308 -b - |\
       samtools sort - > ${n}_adapter.bam
       samtools index ${n}_adapter.bam
       samtools idxstats ${n}_adapter.bam >${n}_adapter.txt
done

```

4.查看比对到adapter上的序列比对上了多少tRNA

```
#! /bin/sh

for i in *_adapter.bam;do
        n=${i/_adapter.bam/}
        samtools view -h ${n}_adapter.bam adapter1_adapter2_1|samtools view -bS - |samtools fastq - >${n}_adapter12.fq
        samtools view -h ${n}_adapter.bam adapter1_1|samtools view -bS - |samtools fastq - >${n}_adapter1.fq
        samtools view -h ${n}_adapter.bam adapter2_1|samtools view -bS - |samtools fastq - >${n}_adapter2.fq
done

#! /bin/sh
ref="/home/jj2024/00.db/01.hg38/03.hg38_trna_nomt/human_tRNA_nomt.fa"

for i in *.fq;do
       n=${i/.fq/}
       bwa mem -W13 -k6 -xont2d -T20 -M $ref ${n}.fq |\
       samtools view -q0 -F2308 -b - |\
       samtools sort - > ${n}_nomt.bam
       samtools index ${n}_nomt.bam
       samtools idxstats ${n}_nomt.bam >${n}_trna_q0_nomt.txt
done

```

