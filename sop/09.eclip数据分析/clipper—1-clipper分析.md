clipper—1-clipper分析

1.clipper环境构建,之前详细过程在我的知乎上写了，配上网址[！！！血书clipper安装过程 - 知乎](https://zhuanlan.zhihu.com/p/7046358695)

![屏幕截图 2025-01-08 134225](C:\Users\jiang\Desktop\屏幕截图 2025-01-08 134225.png)

2.其他软件的安装借助conda

```
##所需安装软件：
cutadapt
umi_tools
STAR
gatk
wf_get_reproducible_eclip_peaks
```

3.clipper软件环境及版本号

```
# packages in environment at /home/gongweikang/.conda/envs/clipper3:

#

# Name                    Version                   Build  Channel

_libgcc_mutex             0.1                        main    conda-forge
alabaster                 0.7.13                   pypi_0    pypi
babel                     2.14.0                   pypi_0    pypi
bedtools                  2.29.2               hc088bd4_0    bioconda
blas                      1.0                         mkl  
bowtie2                   2.4.2            py37h8270d21_1    bioconda
bz2file                   0.98             py37h06a4308_1  
bzip2                     1.0.8                h7b6447c_0  
ca-certificates           2024.11.26           h06a4308_0  
certifi                   2022.12.7        py37h06a4308_0  
charset-normalizer        3.4.1                    pypi_0    pypi
clipper                   2.0.0                    pypi_0    pypi
curl                      7.71.0               hbc83047_0  
cutadapt                  2.6              py37h516909a_0    bioconda
cwl                       0.0.1                    pypi_0    pypi
cycler                    0.10.0                   py37_0  
cython                    0.29.20          py37he6710b0_0  
dbus                      1.13.16              hb2f20db_0  
deeptools                 3.5.4.post1              pypi_0    pypi
deeptoolsintervals        0.1.9            py37hf01694f_2    bioconda
dnaio                     0.3              py37h14c3975_1    bioconda
docutils                  0.19                     pypi_0    pypi
expat                     2.2.9                he6710b0_2  
fontconfig                2.13.0               h9420a91_0  
fonttools                 4.38.0                   pypi_0    pypi
freetype                  2.10.2               h5ab3b9f_0  
future                    0.18.3           py37h06a4308_0  
gatk                      3.6                           5    bioconda
gatk4                     4.0.5.1                       0    bioconda
glib                      2.65.0               h3eb4bd4_0  
gst-plugins-base          1.14.0               hbbd80ab_1  
gstreamer                 1.14.0               hb31296c_0  
htseq                     0.11.3           py37hb3f55d8_0    bioconda
icu                       58.2                 he6710b0_3  
idna                      3.10                     pypi_0    pypi
imagesize                 1.4.1                    pypi_0    pypi
importlib-metadata        6.7.0                    pypi_0    pypi
intel-openmp              2020.1                      217  
java-jdk                  8.0.92                        1    bioconda
jinja2                    3.1.5                    pypi_0    pypi
joblib                    0.15.1                     py_0    conda-forge
jpeg                      9b                   h024ee3a_2  
kiwisolver                1.2.0            py37hfd86e86_0  
krb5                      1.18.2               h173b8e3_0  
ld_impl_linux-64          2.33.1               h53a641e_7    conda-forge
libcurl                   7.71.0               h20c2e04_0  
libdeflate                1.0                  h14c3975_1    bioconda
libedit                   3.1.20191231         h7b6447c_0  
libffi                    3.3                  he6710b0_1  
libgcc                    7.2.0                h69d50b8_2    conda-forge
libgcc-ng                 9.1.0                hdf63c60_0  
libgfortran-ng            7.3.0                hdf63c60_0  
libpng                    1.6.37               hbc83047_0  
libssh2                   1.9.0                h1ba5d50_1  
libstdcxx-ng              9.1.0                hdf63c60_0  
libuuid                   1.0.3                h1bed415_2  
libxcb                    1.13                 h1bed415_1  
libxml2                   2.9.10               he19cac6_1  
markupsafe                2.1.5                    pypi_0    pypi
matplotlib                3.5.3                    pypi_0    pypi
mkl                       2020.1                      217  
mkl-service               2.3.0            py37he904b0f_0  
mkl_fft                   1.1.0            py37h23d657b_0  
mkl_random                1.1.1            py37h0573a6f_0  
ncurses                   6.2                  he6710b0_1  
numpy                     1.18.5           py37ha1c710e_0  
numpy-base                1.18.5           py37hde5b4d6_0  
numpydoc                  1.5.0                    pypi_0    pypi
openjdk                   8.0.152              h7b6447c_3  
openssl                   1.1.1w               h7f8727e_0  
packaging                 24.0                     pypi_0    pypi
pandas                    1.0.5            py37h0573a6f_0  
pcre                      8.44                 he6710b0_0  
perl                      5.26.2               h14c3975_0  
pigz                      2.6                  h27cfd23_0  
pillow                    9.5.0                    pypi_0    pypi
pip                       20.1.1                   py37_1  
plotly                    5.9.0            py37h06a4308_0  
py2bit                    0.3.0            py37hf01694f_4    bioconda
pybedtools                0.8.1            py37h4ef193e_1    bioconda
pybigwig                  0.3.17           py37hfb91e9a_1    bioconda
pygments                  2.17.2                   pypi_0    pypi
pyparsing                 2.4.7                      py_0  
pyqt                      5.9.2            py37h05f1152_2  
pysam                     0.15.3           py37hda2845c_1    bioconda
python                    3.7.7                hcff3b4d_5  
python-dateutil           2.8.1                      py_0    conda-forge
python_abi                3.7                     1_cp37m    conda-forge
pytz                      2020.1                     py_0  
qt                        5.9.7                h5867ecd_1  
readline                  8.0                  h7b6447c_0  
regex                     2022.3.15        py37h7f8727e_0  
requests                  2.31.0                   pypi_0    pypi
samtools                  1.7                           1    bioconda
scikit-learn              0.23.1           py37h423224d_0  
scipy                     1.5.0            py37h0b6359f_0  
setuptools                47.3.1                   py37_0  
sip                       4.19.8           py37hf484d3e_0  
six                       1.15.0                     py_0  
snowballstemmer           2.2.0                    pypi_0    pypi
sphinx                    5.3.0                    pypi_0    pypi
sphinxcontrib-applehelp   1.0.2                    pypi_0    pypi
sphinxcontrib-devhelp     1.0.2                    pypi_0    pypi
sphinxcontrib-htmlhelp    2.0.0                    pypi_0    pypi
sphinxcontrib-jsmath      1.0.1                    pypi_0    pypi
sphinxcontrib-qthelp      1.0.3                    pypi_0    pypi
sphinxcontrib-serializinghtml 1.1.5                    pypi_0    pypi
sqlite                    3.32.3               h62c20be_0  
star                      2.7.10b              h9ee0642_0    bioconda
tbb                       2020.3               hfd86e86_0  
tenacity                  8.0.1            py37h06a4308_1  
threadpoolctl             2.1.0              pyh5ca1d4c_0    conda-forge
tk                        8.6.10               hbc83047_0  
tornado                   6.0.4            py37h7b6447c_1  
typing-extensions         4.7.1                    pypi_0    pypi
umi_tools                 1.0.1            py37hf01694f_2    bioconda
urllib3                   2.0.7                    pypi_0    pypi
wheel                     0.34.2                   py37_0    conda-forge
xopen                     0.7.3                      py_0    bioconda
xz                        5.2.5                h7b6447c_0  
zipp                      3.15.0                   pypi_0    pypi
zlib                      1.2.11               h7b6447c_3  
```

4. 从 FASTQ 文件中提取 Unique Molecular Identifiers (UMIs)。UMIs 是在文库制备过程中添加到每个原始  DNA 或 RNA 分子上的短的随机核苷酸序列，用于区分 PCR 扩增产生的重复序列和原始分子。通过使用  UMIs，可以在后续的数据分析中识别并去除 PCR 引入的偏差，从而提高测序数据的质量和准确性。

   ```
   umi_tools extract --stdin ${n}_1.fq.gz --bc-pattern NNNNNNNNNN --log ${n}_processed.log --stdout ${n}_umi.fq
   ```

5. 去除特定的接头（adapter）和引物序列，并根据设定的质量和长度标准来过滤读段

   ```
   cutadapt --match-read-wildcards --quality-cutoff 6 --times 1 -e 0.1 -O 1 -m 18 -a NNNNNAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -g CTTCCGATCTACAAGTT -g CTTCCGATCTTGGTCCT -o ${n}_umi_trim.fq ${n}_umi.fq > ${n}.metrics
   ```

6. 比对重复序列原件，取出未必对上的读段

   	STAR --runMode alignReads --runThreadN 8 --genomeDir ${index} \
   		--genomeLoad LoadAndRemove --readFilesIn ${n}_umi_trim.fq \
   		--outSAMunmapped Within --outFilterMultimapNmax 30 --outFilterMultimapScoreRange 1 \
   		--outFileNamePrefix ${n} \
   		--outSAMattributes All --outStd BAM_Unsorted --outSAMtype BAM Unsorted \
   		--outFilterType BySJout --outReadsUnmapped Fastx --outFilterScoreMin 10 --outSAMattrRGline ID:foo --alignEndsType EndToEnd

7. 未必对上的的读段比对全基因组进行重新比对（因为我们还要特定病毒序列，因此这一步要分两步，同步进行

   ```
   #比对全基因组序列（由于后续clipper软件原因需要下载ucsc上的hg19序列，其他的不满足要求，会出错无法进行）
   	STAR --runMode alignReads --runThreadN 9 --genomeDir ${index}  \
   	     --genomeLoad LoadAndRemove --readFilesIn ${n}.out.mate1 \
   	     --outSAMunmapped Within --outFilterMultimapNmax 1 --outFilterMultimapScoreRange 1 \
   	     --outFileNamePrefix ${n}_unsorted. \
   	     --outSAMattributes All --outStd BAM_Unsorted --outSAMtype BAM Unsorted \
   	     --outFilterType BySJout --outReadsUnmapped Fastx \
   	     --outFilterScoreMin 10 --outSAMattrRGline ID:foo --alignEndsType EndToEnd >${n}_unsorted.bam
   #比对病毒基因组序列
   for i in *.out.mate1;do
   	n=${i/.out.mate1/}
   	STAR --runMode alignReads --runThreadN 9 --genomeDir ${index}  \
   	     --genomeLoad LoadAndRemove --readFilesIn ${n}.out.mate1 \
   	     --outSAMunmapped Within --outFilterMultimapNmax 1 --outFilterMultimapScoreRange 1 \
   	     --outFileNamePrefix ${n}_unsorted. \
   	     --outSAMattributes All --outStd BAM_Unsorted --outSAMtype BAM Unsorted \
   	     --outFilterType BySJout --outReadsUnmapped Fastx \
   	     --outFilterScoreMin 10 --outSAMattrRGline ID:foo --alignEndsType Local >${n}_unsorted.bam
   ```

8. 对第7步每一步都要做一下处理，首先是对bam进行排序

   ```
   gatk --java-options "-XX:+UseParallelOldGC -XX:ParallelGCThreads=4 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" SortSam --INPUT=${n}.bam --OUTPUT=${n}_sorted.bam --VALIDATION_STRINGENCY=SILENT -SO=coordinate --CREATE_INDEX=true
   ```

9. 其次是利用umitools对bam进行去除PCR重复

   ```
   umi_tools dedup -I ${n}_sorted.bam --output-stats=deduplicatd -S ${n}_rmRep_sorted.bam
   ```

10. 并再次进行bam文件排序

    ```
    gatk --java-options "-XX:+UseParallelOldGC -XX:ParallelGCThreads=4 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" SortSam --INPUT=${n}_rmRep_sorted.bam --OUTPUT=${n}_rmRep_gatksorted.bam --VALIDATION_STRINGENCY=SILENT -SO=coordinate --CREATE_INDEX=true
    ```

11.将排序好的bam文件利用clipper进行call peak

```
 clipper -b ${n}_rmRep_gatksorted.bam -s hg19 -o ${n}.peaks.bed --processors 10
```

12.利用软件对peak文件进行注释，可以用homer或者R语言，以下列出homer的使用

```
/home/gongweikang/.conda/envs/chipseq/share/homer/bin/annotatePeaks.pl ${n}_c3.0_cond1.bed hg38 > ${n}_annotation.xls 
```









