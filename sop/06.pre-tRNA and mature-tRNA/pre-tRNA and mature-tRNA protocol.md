pre-tRNA and mature-tRNA protocol

1.对数据进行初步质控，两种方式

```
##bbduk软件
bbduk.sh in=PrectR0923CTC-LGJ24582_L1_1.fq.gz in2=PrectR0923CTC-LGJ24582_L1_2.fq.gz ref=adapter.fa out=CTC_trim_1.fq.gz out2=CTC_trim_2.fq.gz trimq=30 minlength=50

##trim软件
trim_galore -q 30 --length 50  --illumina -paired ${n}_1.fq.gz ${n}_2.fq.gz
```

2.对质控完的数据，进行mature-tRNA比对

```
#! /bin/sh

# index="/home/jj2024/00.db/01.hg38/05.hg38_trna_illumina/01.mature_tRNA/hg38-mature-tRNAs-CCA-cluster"
index="/home/jj2024/00.db/01.hg38/05.hg38_trna_illumina/01.mature_tRNA/03.addrRNA_snRNA_ivt3_dm_miala_ugc/hg38-mature-tRNAs-CCA-cluster"

for i in *.1_val_1.fq.gz;do
    j=${i/.1_val_1.fq.gz/}
        bowtie2 -p 8 --end-to-end --very-sensitive --dovetail  --no-mixed --no-discordant --fr -k 20 -X 1000  -x ${index} -1 ${j}.1_val_1.fq.gz -2 ${j}.2_val_2.fq.gz  | \
        samtools view -h -bF 4 - | \
        samtools sort -@ 4 - > ${j}.bam
done
```

3.对mature-tRNA进行reads计数

```
salmon quant -t /home/jj2024/00.db/01.hg38/05.hg38_trna_illumina/01.mature_tRNA/03.addrRNA_snRNA_ivt3_dm_miala_ugc/hg38-mature-tRNAs-CCA-cluster.fa -l A -a PrectR0923PIC-LGJ24583_L1.bam --output 01.pic
```

4.先取未必对上mature-tRNA的数据，再进行下一步的pre-tRNA比对

```
#! /bin/sh

index="/home/jj2024/00.db/01.hg38/05.hg38_trna_illumina/01.mature_tRNA/03.addrRNA_snRNA_ivt3_dm_miala_ugc/hg38-mature-tRNAs-CCA-cluster"

for i in *_1_val_1.fq.gz;do
    j=${i/_1_val_1.fq.gz/}
        bowtie2 -p 8 --end-to-end --very-sensitive --dovetail  --no-mixed --no-discordant --fr -k 20 -X 1000  -x ${index} -1 ${j}_1_val_1.fq.gz -2 ${j}_2_val_2.fq.gz  | \
        samtools view -h -bf 4 - | \
        samtools sort -@ 4 - > ${j}.bam
done

#samtools fastq -1 PrectR0923CTC-LGJ24582_L1_1.fq -2 PrectR0923CTC-LGJ24582_L1_2.fq  -n PrectR0923CTC-LGJ24582_L1.bam
```

5.pre-tRNA比对

```
salmon quant -i /home/jj2024/00.db/01.hg38/05.hg38_trna_illumina/02.Pre-tRNA/08.add_50bp_rRNA_snRNA_ivt3_dm_miala_ugc/hg38-tRNAs-50bp -l A -1 <fq> -2 <fq> --output <name>
```

