#on SFTP
查看服务器上的软件
module list 
启动anaconda
module load anaconda
查看目前有的环境
conda env list
conda info --evn
选择一个conda环境启用
conda activate tf#tf：文件名称
#看环境里装有哪些软件
conda list
#有umitools等

#去除adapter可以直接在base环境即可
conda actibate base

cutadapt -f fastq --match-read-wildcards --quality-cutoff 6 --times 1 -e 0.1 -O 1 -m 18 -a NNNNNAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -g CTTCCGATCTACAAGTT -g CTTCCGATCTTGGTCCT -o PTBP1-IP1.trim.read1.fastq PTBP1-IP1_L1_1.fq.gz >PTBP1-IP1.metrics

cutadapt -f fastq --match-read-wildcards --quality-cutoff 6 --times 1 -e 0.1 -O 1 -m 18 -a NNNNNAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -g CTTCCGATCTACAAGTT -g CTTCCGATCTTGGTCCT -o PTBP1-IP2.trim.read1.fastq PTBP1-IP2_L1_1.fq.gz >PTBP1-IP2.metrics

cutadapt -f fastq --match-read-wildcards --quality-cutoff 6 --times 1 -e 0.1 -O 1 -m 18 -a NNNNNAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -g CTTCCGATCTACAAGTT -g CTTCCGATCTTGGTCCT -o PTBP1-Input1.trim.read1.fastq P-Input1_L1_1.fq.gz >PTBP1-Input1.metrics

cutadapt -f fastq --match-read-wildcards --quality-cutoff 6 --times 1 -e 0.1 -O 1 -m 18 -a NNNNNAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -g CTTCCGATCTACAAGTT -g CTTCCGATCTTGGTCCT -o PTBP1-Input2.trim.read1.fastq P-Input2_L1_1.fq.gz >PTBP1-Input2.metrics

# umi treatment
#需要回到有umi_tools的环境 CLIP_env
umi_tools extract --stdin=PTBP1-IP1.trim.read1.fastq --bc-pattern=NNNNNNNNNN --log=./processed_PTBP1-IP_trimed_read1.log --stdout ./processed_PTBP1-IP_trimed_read1.fastq.gz

umi_tools extract --stdin=PTBP1-IP2.trim.read1.fastq --bc-pattern=NNNNNNNNNN --log=./processed_PTBP1-IP_trimed_read2.log --stdout ./processed_PTBP1-IP_trimed_read2.fastq.gz

umi_tools extract --stdin=PTBP1-Input1.trim.read1.fastq --bc-pattern=NNNNNNNNNN --log=./processed_PTBP1-Input_trimed_read1.log --stdout ./processed_PTBP1-Input_trimed_read1.fastq.gz

umi_tools extract --stdin=PTBP1-Input2.trim.read1.fastq --bc-pattern=NNNNNNNNNN --log=./processed_PTBP1-Input_trimed_read2.log --stdout ./processed_PTBP1-Input_trimed_read2.fastq.gz


#module load star

nohup STAR --runMode alignReads --runThreadN 8 --genomeDir /cluster/facility/cxiao/pd/rep_index/ --genomeLoad LoadAndRemove --readFilesIn processed_PTBP1-IP_trimed_read1.fastq.gz --outSAMunmapped Within --outFilterMultimapNmax 30 --outFilterMultimapScoreRange 1 --outFileNamePrefix PTBP1-IP1_trimed_R1.rep.bam --readFilesCommand zcat --outSAMattributes All --outStd BAM_Unsorted --outSAMtype BAM Unsorted --outFilterType BySJout --outReadsUnmapped Fastx --outFilterScoreMin 10 --outSAMattrRGline ID:foo --alignEndsType EndToEnd > processed.PTBP1_IP1_trimed_R1.rep.bam

nohup STAR --runMode alignReads --runThreadN 8 --genomeDir /cluster/facility/cxiao/pd/rep_index/ --genomeLoad LoadAndRemove --readFilesIn processed_PTBP1-IP_trimed_read2.fastq.gz --outSAMunmapped Within --outFilterMultimapNmax 30 --outFilterMultimapScoreRange 1 --outFileNamePrefix PTBP1-IP2_trimed_R1.rep.bam --readFilesCommand zcat --outSAMattributes All --outStd BAM_Unsorted --outSAMtype BAM Unsorted --outFilterType BySJout --outReadsUnmapped Fastx --outFilterScoreMin 10 --outSAMattrRGline ID:foo --alignEndsType EndToEnd > processed.PTBP1_IP2_trimed_R1.rep.bam

STAR --runMode alignReads --runThreadN 8 --genomeDir /cluster/facility/cxiao/pd/rep_index/ --genomeLoad LoadAndRemove --readFilesIn processed_PTBP1-Input_trimed_read1.fastq.gz --outSAMunmapped Within --outFilterMultimapNmax 30 --outFilterMultimapScoreRange 1 --outFileNamePrefix PTBP1-Input1_trimed_R1.rep.bam --readFilesCommand zcat --outSAMattributes All --outStd BAM_Unsorted --outSAMtype BAM Unsorted --outFilterType BySJout --outReadsUnmapped Fastx --outFilterScoreMin 10 --outSAMattrRGline ID:foo --alignEndsType EndToEnd > processed.PTBP1_Input1_trimed_R1.rep.bam

STAR --runMode alignReads --runThreadN 8 --genomeDir /cluster/facility/cxiao/pd/rep_index/ --genomeLoad LoadAndRemove --readFilesIn processed_PTBP1-Input_trimed_read2.fastq.gz --outSAMunmapped Within --outFilterMultimapNmax 30 --outFilterMultimapScoreRange 1 --outFileNamePrefix PTBP1-Input2_trimed_R1.rep.bam --readFilesCommand zcat --outSAMattributes All --outStd BAM_Unsorted --outSAMtype BAM Unsorted --outFilterType BySJout --outReadsUnmapped Fastx --outFilterScoreMin 10 --outSAMattrRGline ID:foo --alignEndsType EndToEnd > processed.PTBP1_Input2_trimed_R1.rep.bam

#align genome
#新建一个文件夹Align
#mkdir Align
nohup STAR --runMode alignReads --runThreadN 9 --genomeDir /cluster/database/Index_old/STAR-Index/mm10/ --genomeLoad LoadAndRemove --readFilesIn PTBP1-IP1_trimed_R1.rep.bamUnmapped.out.mate1 --outSAMunmapped Within --outFilterMultimapNmax 1 --outFilterMultimapScoreRange 1 --outFileNamePrefix ./Align/PTBP1-IP1.align.rmRep.bam --outSAMattributes All --outStd BAM_Unsorted --outSAMtype BAM Unsorted --outFilterType BySJout --outReadsUnmapped Fastx --outFilterScoreMin 10 --outSAMattrRGline ID:foo --alignEndsType EndToEnd > ./Align/PTBP1-IP1.align.rmRep.bam

STAR --runMode alignReads --runThreadN 9 --genomeDir /cluster/database/Index_old/STAR-Index/mm10/ --genomeLoad LoadAndRemove --readFilesIn PTBP1-IP2_trimed_R1.rep.bamUnmapped.out.mate1 --outSAMunmapped Within --outFilterMultimapNmax 1 --outFilterMultimapScoreRange 1 --outFileNamePrefix ./Align/PTBP1-IP2.align.rmRep.bam --outSAMattributes All --outStd BAM_Unsorted --outSAMtype BAM Unsorted --outFilterType BySJout --outReadsUnmapped Fastx --outFilterScoreMin 10 --outSAMattrRGline ID:foo --alignEndsType EndToEnd > ./Align/PTBP1-IP2.align.rmRep.bam

STAR --runMode alignReads --runThreadN 9 --genomeDir /cluster/database/Index_old/STAR-Index/mm10/ --genomeLoad LoadAndRemove --readFilesIn PTBP1-Input1_trimed_R1.rep.bamUnmapped.out.mate1 --outSAMunmapped Within --outFilterMultimapNmax 1 --outFilterMultimapScoreRange 1 --outFileNamePrefix ./Align/PTBP1-Input1.align.rmRep.bam --outSAMattributes All --outStd BAM_Unsorted --outSAMtype BAM Unsorted --outFilterType BySJout --outReadsUnmapped Fastx --outFilterScoreMin 10 --outSAMattrRGline ID:foo --alignEndsType EndToEnd > ./Align/PTBP1-Input1.align.rmRep.bam

STAR --runMode alignReads --runThreadN 9 --genomeDir /cluster/database/Index_old/STAR-Index/mm10/ --genomeLoad LoadAndRemove --readFilesIn PTBP1-Input2_trimed_R1.rep.bamUnmapped.out.mate1 --outSAMunmapped Within --outFilterMultimapNmax 1 --outFilterMultimapScoreRange 1 --outFileNamePrefix ./Align/PTBP1-Input2.align.rmRep.bam --outSAMattributes All --outStd BAM_Unsorted --outSAMtype BAM Unsorted --outFilterType BySJout --outReadsUnmapped Fastx --outFilterScoreMin 10 --outSAMattrRGline ID:foo --alignEndsType EndToEnd > ./Align/PTBP1-Input2.align.rmRep.bam

#module load samtools+htslib
#module load anaconda/201903
#mkdir sorted
/cluster/facility/cxiao/pd/gatk-4.1.3.0/gatk --java-options "-Xmx2048m -XX:+UseParallelOldGC -XX:ParallelGCThreads=4 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" SortSam --INPUT=PTBP1-IP1.align.rmRep.bam --OUTPUT=./sorted/PTBP1-IP1.align.rmRep.sorted.bam --VALIDATION_STRINGENCY=SILENT -SO=coordinate --CREATE_INDEX=true

/cluster/facility/cxiao/pd/gatk-4.1.3.0/gatk --java-options "-Xmx2048m -XX:+UseParallelOldGC -XX:ParallelGCThreads=4 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" SortSam --INPUT=PTBP1-IP2.align.rmRep.bam --OUTPUT=./sorted/PTBP1-IP2.align.rmRep.sorted.bam --VALIDATION_STRINGENCY=SILENT -SO=coordinate --CREATE_INDEX=true

/cluster/facility/cxiao/pd/gatk-4.1.3.0/gatk --java-options "-Xmx2048m -XX:+UseParallelOldGC -XX:ParallelGCThreads=4 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" SortSam --INPUT=PTBP1-Input1.align.rmRep.bam --OUTPUT=./sorted/PTBP1-Input1.align.rmRep.sorted.bam --VALIDATION_STRINGENCY=SILENT -SO=coordinate --CREATE_INDEX=true

/cluster/facility/cxiao/pd/gatk-4.1.3.0/gatk --java-options "-Xmx2048m -XX:+UseParallelOldGC -XX:ParallelGCThreads=4 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" SortSam --INPUT=PTBP1-Input2.align.rmRep.bam --OUTPUT=./sorted/PTBP1-Input2.align.rmRep.sorted.bam --VALIDATION_STRINGENCY=SILENT -SO=coordinate --CREATE_INDEX=true

#需要回到有umi_tools的环境 CLIP_env

umi_tools dedup -I PTBP1-IP1.align.rmRep.sorted.bam --output-stats=deduplicatd -S deduplicated.PTBP1-IP1.align.rmRep.sorted.bam

umi_tools dedup -I PTBP1-IP2.align.rmRep.sorted.bam --output-stats=deduplicatd -S deduplicated.PTBP1-IP2.align.rmRep.sorted.bam

umi_tools dedup -I PTBP1-Input1.align.rmRep.sorted.bam --output-stats=deduplicatd -S deduplicated.PTBP1-Input1.align.rmRep.sorted.bam

umi_tools dedup -I PTBP1-Input2.align.rmRep.sorted.bam --output-stats=deduplicatd -S deduplicated.PTBP1-Input2.align.rmRep.sorted.bam



# sortSam
/cluster/facility/cxiao/pd/gatk-4.1.3.0/gatk --java-options "-Xmx2048m -XX:+UseParallelOldGC -XX:ParallelGCThreads=4 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" SortSam --INPUT=deduplicated.PTBP1-IP1.align.rmRep.sorted.bam --OUTPUT=deduplicated.PTBP1-IP1.align.rmRep.final.bam --VALIDATION_STRINGENCY=SILENT -SO=coordinate --CREATE_INDEX=true

/cluster/facility/cxiao/pd/gatk-4.1.3.0/gatk --java-options "-Xmx2048m -XX:+UseParallelOldGC -XX:ParallelGCThreads=4 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" SortSam --INPUT=deduplicated.PTBP1-IP2.align.rmRep.sorted.bam --OUTPUT=deduplicated.PTBP1-IP2.align.rmRep.final.bam --VALIDATION_STRINGENCY=SILENT -SO=coordinate --CREATE_INDEX=true

/cluster/facility/cxiao/pd/gatk-4.1.3.0/gatk --java-options "-Xmx2048m -XX:+UseParallelOldGC -XX:ParallelGCThreads=4 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" SortSam --INPUT=deduplicated.PTBP1-Input1.align.rmRep.sorted.bam --OUTPUT=deduplicated.PTBP1-Input1.align.rmRep.final.bam --VALIDATION_STRINGENCY=SILENT -SO=coordinate --CREATE_INDEX=true

/cluster/facility/cxiao/pd/gatk-4.1.3.0/gatk --java-options "-Xmx2048m -XX:+UseParallelOldGC -XX:ParallelGCThreads=4 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" SortSam --INPUT=deduplicated.PTBP1-Input2.align.rmRep.sorted.bam --OUTPUT=deduplicated.PTBP1-Input2.align.rmRep.final.bam --VALIDATION_STRINGENCY=SILENT -SO=coordinate --CREATE_INDEX=true


#创建用户自己的环境 clipper 
conda create -n clipper
#激活tf 环境并安装软件包
source activate clipper
conda install numpy（包名）

##加载anaconda 需要高级3以上的版本，适用python3 
module load anaconda/3.5.2
#download clipper-master.zip from github
cd clipper-master/

conda env create -f environment3.yml
conda activate clipper3
pip install clipper 



#on server environment:clipper3 and load bedtools
#source activate clipper3 
nohup clipper -b deduplicated.PTBP1-IP2.align.rmRep.final.bam -s mm10 -o deduplicated.PTBP1-IP2.rmRep.sorted.peaks.bed --processors 10

clipper -b deduplicated.PTBP1-IP1.align.rmRep.final.bam -s mm10 -o deduplicated.PTBP1-IP1.rmRep.sorted.peaks.bed --processors 10


#merge peak

wf_get_reproducible_eclip_peaks.cwl is in your $PATH.
Type: ./merge_peaks_2inputs.yaml

export PATH=${PWD}/bin:$PATH
export PATH=${PWD}/bin/perl:$PATH
export PATH=${PWD}/cwl:$PATH
export PATH=${PWD}/wf:$PATH;

species: mm10
samples:
  - 
    - name: "PTBP1_rep1"
      ip_bam: 
        class: File
        path: /cluster/facility/cxiao/KunyuLiao/20240717PTBP1/Align/sorted/deduplicated.PTBP1-IP1.align.rmRep.final.bam
      input_bam:
        class: File
        path: /cluster/facility/cxiao/KunyuLiao/20240717PTBP1/Align/sorted/deduplicated.PTBP1-Input1.align.rmRep.final.bam
      peak_clusters:
        class: File
        path: /cluster/facility/cxiao/KunyuLiao/20240717PTBP1/Align/sorted/deduplicated.PTBP1-IP1.rmRep.sorted.peaks.bed
    - name: "PTBP1_rep2"
      ip_bam: 
        class: File
        path: /cluster/facility/cxiao/KunyuLiao/20240717PTBP1/Align/sorted/deduplicated.PTBP1-IP2.align.rmRep.final.bam
      input_bam:
        class: File
        path: /cluster/facility/cxiao/KunyuLiao/20240717PTBP1/Align/sorted/deduplicated.PTBP1-Input2.align.rmRep.final.bam
      peak_clusters:
        class: File
        path: /cluster/facility/cxiao/KunyuLiao/20240717PTBP1/Align/sorted/deduplicated.PTBP1-IP2.rmRep.sorted.peaks.bed
chrom_sizes:
  class: File
  path: /cluster/database/Genomes/chrom.sizes/mm10.chrom.sizes



