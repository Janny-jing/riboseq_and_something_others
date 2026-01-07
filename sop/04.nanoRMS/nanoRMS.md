# nanoRMS：同时定位不同RNA修饰

## 安装

​        wget https://github.com/novoalab/nanoRMS/archive/refs/tags/v1.0.tar.gz   &&  tar -zxvf v1.0.tar.gz 

## 环境配置

​       ###docker build -f /home/chen/jiangjing/soft/nanoRMS-1.0/docker/Dockerfile .  （不可用，内置有问题） 

​       ###linux配置所需依赖包

​        conda create -n nanoRMS

​        conda activate nanoRMS

​        conda install -c bioconda ont-tombo

​        conda install conda-forge::hdf5

​        conda install minimap2

​        conda install dask

​        conda install -c bioconda nanopolish

​        pip install "numpy<1.20" ont_fast5_api matplotlib mappy numba pysam pandas seaborn scikit-learn jupyterlab openpyxl mappy scipy

​        ###linux R配置包

​        ###环境：当前环境R版本低，需安装版本大于3.5，但不能大于3.7会和tombo冲突，安装依赖包

​       ①conda uninstall r-base

​       ②conda install conda-forge::r-base

​       ③install.packages(c(“optparse”,”stringr”,”dplyr”,”plyr”,”VennDiagram”,”RColorBrewer”,”data.table”))

## 使用

###         1.   prediction of RNA modified sits

####           1.1 Extract base-calling features using Epinano-RMS

​              ① 对reference fasta构建索引，并且同一目录下需要有fai文件（samtools faidx）

```
java -jar epinano_RMS/picard.jar CreateSequenceDictionary REFERENCE=reference.fasta OUTPUT= reference.fasta.dict
```

​               脚本示例：java -jar /home/chen/jiangjing/soft/nanoRMS/epinano_RMS/picard.jar  CreateSequenceDictionary REFERENCE=yeast_all_rRNA.fa OUTPUT= reference.fasta.dict

​              ② 计算碱基位点频率（！！！如果跑第二遍要确定是否删掉了隐藏文件（.r1.sorted.tmp_splitted_base_freq.done_splitting））

```
python3 epinano_RMS/epinano_rms.py -R <reference_file> -b <bam_file> -s epinano_RMS/sam2tsv
```

​               脚本：python3 /home/chen/jiangjing/soft/nanoRMS/epinano_RMS/epinano_rms.py -R  /home/chen/jiangjing/soft/nanoRMS/predict_rna_mod/references/yeast_all_rRNA.fa -b /home/chen/jiangjing/soft/nanoRMS/predict_rna_mod/test_data/test.bam -s /home/chen/jiangjing/soft/nanoRMS/epinano_RMS/sam2tsv.jar

​               提醒会报错，显示无某个文件夹，直接建一个，后面删掉就行：mkdir -p ./home/chen/jiangjing/soft/nanoRMS/predict_rna_mod/test_data/test.tmp_splitted_base_freq.done_splitting

​                结果文件：包含每个碱基位点频率的.csv文件（test.per.site.baseFreq.csv）

![image-20240413180318310](C:\Users\admin\AppData\Roaming\Typora\typora-user-images\image-20240413180318310.png)

####          1.2 Predict RNA modifications

​              ① Single sample RNA modification prediction (i.e. "de novo" prediction)

​               ##predicting pseudouridine RNA modifications in mitochondrial rRNAs, validating 2 out of the 2 sites that were predicted in all 3 biological replicates.Specifically, **pseudouridine causes strong mismatch signatures** in the modified position, largely **in the form of C-to-U mismatches**

```
Rscript --vanilla Pseudou_prediction_singlecondition.R [options] -f <epinano_file1> (-s <epinano_file2> -t <epinano_file3>)
```


Options:

        -h, --help
                    Show this help message and exit
                    
        -f CHARACTER, --epinano1=CHARACTER
                first epinano file
    
        -s CHARACTER, --epinano2=CHARACTER
                second epinano file
    
        -t CHARACTER, --epinano3=CHARACTER
                third epinano file
    
        -p CHARACTER, --modpos=CHARACTER
                mod positions file [default= positions/RNA_Mod_Positions_rRNAYeast.tsv]
    
        -m NUMBER, --misfreq=NUMBER
                Mismatch frequency threshold [default= 0.137]
    
        -c NUMBER, --Cfreq=NUMBER
                C mismatch frequency threshold [default= 0.578]
​               脚本示例：Rscript --vanilla /home/chen/jiangjing/soft/nanoRMS/predict_rna_mod/Pseudou_prediction_singlecondition.R -f WT_rRNA_Epinano.csv -s sn34KO_rRNA_Epinano.csv -t sn36KO_rRNA_Epinano.csv -p ../positions/RNA_Mod_Positions_rRNAYeast.tsv

​               ##csv文件是1.1的结果文件

​               ##其中mod_position文件是依据已知的假尿嘧啶修饰位点相关数据库和文献，自建的

​               ##mod_position文件相关解释：https://github.com/novoalab/nanoRMS/issues/32

​               ##找到的修饰位点相关数据库：https://rna.sysu.edu.cn/rmbase3/pseudosearch.php

​               结果文件：（1）threesample_predicted_reproducible_y_sites.tsv（三个重复样本交集位点）

​                                  （2）三个重复样本交集venn图：Threesamples_venn_diagramm_mitochondrial.tiff

​                                                       ![image-20240413181334314](C:\Users\admin\AppData\Roaming\Typora\typora-user-images\image-20240413181334314.png)   <img src="C:\Users\admin\AppData\Roaming\Typora\typora-user-images\image-20240413181427338.png" alt="image-20240413181427338" style="zoom:50%;" />



修改：Rscript --vanilla Pseudou_prediction_singlecondition.R -f CTR1226CTR1.per.site.baseFreq.csv -p single.tsv



​              ① Paired sample RNA modification prediction (i.e. "differential-error"-based prediction)

​               ##Pseudouridine is not always present in high stoichiometries (e.g. rRNAs), but can also be present in low stoichiometries (e.g. in mRNAs). Please note that pseudouridines in mRNAs cannot be accurately predicted using "de novo" mode, because the background nanopore 'error' is too similar to the 'error' caused by the presence of pseudouridine. For such cases, we can predict differentially pseudouridylated sites by identifying which sites show pseudouridine differential error signatures between two conditions, as shown below. This type of pairwise comparison can be done for WT-KO, or between two conditions (e.g. normal-heat stress). In the image below, some examples of heat-responsive sites that were identified using this script are shown. Differential-error based prediction can be applied to any type of RNA (mRNA, snoRNAs, snRNAs, rRNAs etc).

​                  1）For Transcriptome mapped reads (Reads from one strand)

```
Rscript --vanilla Pseudou_prediction_pairedcondition_transcript.R [options] -f <epinano_file1> -s <epinano_file2> 
```

Options:

        -h, --help
                    Show this help message and exit
                    
        -f CHARACTER, --epinano1=CHARACTER
                first epinano file
    
        -s CHARACTER, --epinano2=CHARACTER
                second epinano file
    
        -p CHARACTER, --modpos=CHARACTER
                mod positions file    
                [default=positions/RNA_Mod_Positions_ncRNAYeast_HeatSensitive.tsv]
    
        -m NUMBER, --misfreq=NUMBER
                Mismatch frequency threshold [default= 0.137]
    
        -c NUMBER, --Cfreq=NUMBER
                C mismatch frequency threshold [default= 0.578]
    
        -d NUMBER, --diff=NUMBER
                mismatch frequency difference threshold [default= 0.1]

​               脚本示例：Rscript --vanilla ../Pseudou_prediction_pairedcondition_transcript.R -f WT_ncRNA_Normal_Rep1_Epinano.csv -s WT_ncRNA_HeatShock_Rep1_Epinano.csv -p ../positions/RNA_Mod_Positions_ncRNAYeast_HeatSensitive.tsv

​               结果文件：（1）Paired_comparison_altering_sites_pU_predictions.bed 

​                                  （2）Paired_comparison_altering_sites_pU_predictions.tsv

​                ![image-20240413181702273](C:\Users\admin\AppData\Roaming\Typora\typora-user-images\image-20240413181702273.png)



​                  2）For Genome mapped reads (Reads from both strands)

​                        ①csv转bed

```
./Epinano_to_BED.sh <epinano_file1>
```

​                           脚本示例：./Epinano_to_BED.sh WT_mRNA_Normal_Epinano.csv

​                       ②Convert GTF (Only CDS) output into BED

```
Rscript --vanilla GTF_to_BED.R <GTF_File>
```

​                           脚本示例：Rscript --vanilla ../GTF_to_BED.R Saccer64.gtf

​                       ③Intersect Epinano_BED file and GTF_BED file

```
bedtools intersect -a <Epinano_Bed> -b <GTF_Bed> -wa -wb > output.bed
```

​                           脚本示例：bedtools intersect -a ../test_data/WT_mRNA_Normal_Epinano.csv.bed -b Saccer3.bed -wa -wb > WT_mRNA_Normal_Epinano_final.bed

​                        ④mapped

​                          官网给的是错的：

​                          ~~Rscript --vanilla ../Pseudou_prediction_pairedcondition_genome.R -f WT_rRNA_Epinano.csv -s sn36KO_rRNA_Epinano.csv~~

​                           第一步：修改脚本：Pseudou_prediction_pairedcondition_genome.R

![image-20240413145139189](C:\Users\admin\AppData\Roaming\Typora\typora-user-images\image-20240413145139189.png)

​                             mod.file去掉#

​                           第二步：脚本示例：

​                            Rscript --vanilla ../Pseudou_prediction_pairedcondition_genome.R -f WT_mRNA_Normal_Epinano_final.bed -s WT_mRNA_HeatShock_Epinano_final.bed -p ../positions/RNA_Mod_Positions_mRNAYeast_HeatSensitive.tsv

​                            结果文件：（1）Paired_comparison_altering_sites_pU_predictions.bed 

​                                               （2）Paired_comparison_altering_sites_pU_predictions.tsv

​                             ![image-20240413182242511](C:\Users\admin\AppData\Roaming\Typora\typora-user-images\image-20240413182242511.png)



### ~~2. RNA modification stoichiometry estimation using Nanopolish resquiggling（官网不推荐）~~



###         3. RNA modification stoichiometry estimation using Tombo resquiggling

####           3.1 第一步下载相关数据的fast5及fq文件

####           3.2 retrieve per-read features from all sample（fast5转bam）

```
./get_features.py --rna -f reference.fa -t 6 -i .fast5文件路径
```

结果文件：![image-20240413150351578](C:\Users\admin\AppData\Roaming\Typora\typora-user-images\image-20240413150351578.png)



####           3.3 Run nanopolish if you wish to compare it with tombo. You can skip this step if you don't want to compare the tools.

```
ref=reference.fa 
for d in guppy3.0.3.hac/*_{snR*,wt}/; do
 if [ ! -s $d/fq.gz.bam.events.gz ]; then
  echo `date` $d;
  cat $d/*.fastq.gz > $d/fq.gz;
  nanopolish index -d $d $d/fq.gz 2> /dev/null;
  minimap2 -ax map-ont $ref $d/fq.gz 2> /dev/null | samtools sort -o $d/fq.gz.bam;
  samtools index $d/fq.gz.bam;
  nanopolish eventalign -n -t 6 -q 10 --progress --signal-index --scale-events --reads $d/fq.gz --bam $d/fq.gz.bam --genome $ref | gzip > $d/fq.gz.bam.events.gz;
 fi;
done; date
```

结果文件：![image-20240413150545230](C:\Users\admin\AppData\Roaming\Typora\typora-user-images\image-20240413150545230.png)

####           3.4 Estimate modification frequency difference between two samples

####Note, you'll need to provide candidate positions that are likely modified. Those were identified earlier -- please see above section 1.2. Predict RNA modifications. so here we'll just generate BED file from existing candidate file.

```
# prepare BED - this is no longer needed step!
f=per_read/results/predictions_ncRNA_WT30C_WT45C.tsv.gz
zgrep -v X.Ref $f |awk -F'\t' 'BEGIN {OFS = FS} {print $1,$2-1,$2,".",100,"+"}' > $f.bed
```

bed格式：

![image-20240413150934959](C:\Users\admin\AppData\Roaming\Typora\typora-user-images\image-20240413150934959.png)

```
# calculate modification frequency difference between control and sample of interest
per_read/get_freq.py -f $ref -b $f.bed -o $f.bed.tsv.gz -1 per_read/guppy3.0.3.hac/*WT30C/workspace/*.fast5.bam -2 per_read/guppy3.0.3.hac/*WT45C/workspace/*.fast5.bam
```

usage: get_freq.py

option：

    -h/--help
                Show this help message and exit
                
    --version/-v 
    
    -1 control...
    
    -2 sample...
    
    -o output
    
    -b BED文件
    
    -f/--fasta 
    
    -m mincov (required minimum number of reads,but it has to be at least 5 due to KNN requirements.)
结果文件：predictions_ncRNA_WT30C_WT45C.bed.tsv.gz

(KNN 的最终目的是分类，而 K-means 的目的是给所有距离相近的点分配一个类别，也就是聚类)

![image-20240413151828320](C:\Users\admin\AppData\Roaming\Typora\typora-user-images\image-20240413151828320.png)

###         4. Visualization of per-read current intensities at individual sites

####           4.1 Data pre-processing

​              ###Firstly, generate a collapsed Nanopolish event align output(3.3结果文件), by collapsing all the multiple observations for a given position from a same read.

```
python3 per_read_mean.py <event_align_file>
```

​               脚本示例：python3 per_read_mean.py test_data/data1_eventalign_output.txt

​               输出文件：data1_eventalign_output.txt_processed_perpos_mean.csv

​                ![image-20240413161024642](C:\Users\admin\AppData\Roaming\Typora\typora-user-images\image-20240413161024642.png)

​              ###Secondly, create 15-mer windows of per-read current intensities centered in positions of interest

```
Rscript --vanilla nanopolish_window.R positions_file <input_table> <label>
```

​               脚本示例：Rscript --vanilla nanopolish_window.R test_data/positions.tsv test_data/data1_eventalign_output.txt_processed_perpos_mean.csv data1

​               输出文件：data1_window_file.tsv：

​                ![image-20240413161209559](C:\Users\admin\AppData\Roaming\Typora\typora-user-images\image-20240413161209559.png)

####           4.2. Visualization of current intensity information:

#####             4.2.1 Density plots

```
Rscript --vanilla density_nanopolish.R <window_file1> <window_file2> <window_file3(optional)> <window_file4(optional)>
```

​               脚本示例：Rscript --vanilla nanopolish_density_plot.R bc2_window_file.tsv bc4_window_file.tsv

​               输出文件：

![image-20240413161348221](C:\Users\admin\AppData\Roaming\Typora\typora-user-images\image-20240413161348221.png)

#####            4.2.2. Mean current intensity plots centered in the modified sites

```
Rscript --vanilla nanopolish_meanlineplot.R <window_file1> <window_file2> <window_file3(optional)> <window_file4(optional)>
```

​               脚本示例：Rscript --vanilla nanopolish_meanlineplot.R bc2_window_file.tsv bc4_window_file.tsv

​               输出文件：

​                 ![image-20240413161513432](C:\Users\admin\AppData\Roaming\Typora\typora-user-images\image-20240413161513432.png)

#####            4.2.3 Per-read current intensity plots centered in the modified sites

```
Rscript --vanilla nanopolish_perreadlineplot.R <window_file1> <window_file2> <window_file3(optional)> <window_file4(optional)>
```

​               脚本示例：Rscript --vanilla nanopolish_perreadlineplot.R bc2_window_file.tsv bc4_window_file.tsv

​               输出文件：

​                ![image-20240413161559274](C:\Users\admin\AppData\Roaming\Typora\typora-user-images\image-20240413161559274.png)

#####            4.2.4 PCA plots from the per-read 15-mer current intensity data

```
Rscript --vanilla nanopolish_pca.R <window_file1.tsv> <window_file2.tsv> <window_file3.tsv(optional)> <window_file4.tsv(optional)>
```

​               脚本示例：Rscript --vanilla nanopolish_pca.R bc2_window_file.tsv bc4_window_file.tsv

​               输出文件：

​               ![image-20240413161645509](C:\Users\admin\AppData\Roaming\Typora\typora-user-images\image-20240413161645509.png)





示例：

01.bw.sh

#! /bin/sh

path="/home/chen/jiangjing/02.nanoRMS/"
ref="/home/chen/jiangjing/00.db/hg38/hg38-tRNAs.fa"
for i in CTR1226CTR1 CTR1226PIC1;do
       bwa mem -W13 -k6 -xont2d -T20 $ref $path/$i/fastq/merge.fastq.gz >$i.sam
       samtools view -b $i.sam >$i_1.bam
       samtools sort $i_1.bam > ${i}.bam
       samtools index ${i}.bam
       rm $i.sam $i_1.bam
done





02.nanorms.sh

#! /bin/sh

ref="/home/chen/jiangjing/00.db/hg38/hg38-tRNAs.fa"
for i in *.bam;do
        python3 /home/chen/jiangjing/soft/nanoRMS/epinano_RMS/epinano_rms.py -R $ref -b $i -s /home/chen/jiangjing/soft/nanoRMS/epinano_RMS/sam2tsv.jar
done



03.nanopolish.sh

#! /bin/sh

ref="/home/chen/jiangjing/00.db/hg38/hg38-tRNAs.fa"
for d in /home/chen/jiangjing/02.nanoRMS/CTR1226CTR1 /home/chen/jiangjing/02.nanoRMS/CTR1226PIC1; do
  nanopolish index -d $d/fast5 $d/fastq/merge.fastq.gz 2> /dev/null;
  minimap2 -ax map-ont $ref $d/fastq/merge.fastq.gz 2> /dev/null | samtools sort -o $d/fastq/fq.gz.bam;
  samtools index $d/fastq/fq.gz.bam;
  nanopolish eventalign -n -t 6 -q 10 --progress --signal-index --scale-events --reads $d/fastq/merge.fastq.gz --bam $d/fastq/fq.gz.bam --genome $ref | gzip > $d/fastq/fq.gz.bam.events.gz;
done





04.fast5_bam.sh

#! /bin/sh

ref="/home/chen/jiangjing/00.db/hg38/hg38-tRNAs.fa"
/home/chen/jiangjing/soft/nanoRMS/per_read/get_features.py --rna -f $ref -t 6 -i /home/chen/jiangjing/02.nanoRMS/CTR*/fast5/*.fast5