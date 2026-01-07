rm(list = ls())
# BiocManager::install("DiffBind")
library(readxl)
library(DiffBind)
library(tidyverse)
library(magrittr)
library(profileplyr)
setwd("E:\\temp\\03.atacseq")
data <- read.csv("diffbind.csv")
ta <- dba(sampleSheet = data)
# ta <- dba.blacklist(ta, blacklist=DBA_BLACKLIST_GRCH38)
#看样本间的相关性
# plot(ta,main="PCA_heatmap")
p <- dba.plotHeatmap(ta,margin=12)
ggsave("PCA_heatmap.png", plot=p, width=6, height=4)

#亲和结合矩阵，minOverlap最少重叠数，summits区间
ta_count <- dba.count(ta,minOverlap = 2)
dba.plotPCA(ta_count, attributes=DBA_TREATMENT, label=DBA_ID,ylim=c(-100,100),xlim=c(-300,300))

ta <- dba.normalize(ta_count)
# Establishing a contrast
#分组1，后面使用contrast=1单独查看
#name可以更改样品名
test <- dba.contrast(ta, contrast = c("Treatment","treat","control"), name1 = "polyic", name2 = "control")
#按照分组分别进行差异分析，默认使用DESeq2进行计算，可以选择method = DBA_EDGER(edgR)，或者两个都要method = DBA_ALL_METHODS
ta_analyze <- dba.analyze(test, method=DBA_DESEQ2)

# summary of results
dba.show(ta_analyze, bContrasts=T)
#保存文件
ta.DB_1 <- dba.report(ta_analyze,method=DBA_DESEQ2)
write.csv(ta.DB_1, file = "diffbind_polyic_control.csv")

#Rds文件内具体参数保存可以利用@
saveRDS(ta_analyze, file = "ta_analyze.Rds")

#bLoess中间的线,bFlip改变y轴(正负对调)，bxy改变中间线(水平或者斜线)
dba.plotMA(ta_analyze, contrast=1,method=DBA_DESEQ2, dotSize=.90, cex.axis=1, cex.lab=1.3, 
           bNormalized = FALSE, bLoess = 0, bFlip = FALSE,bXY = FALSE)

# overlapping peaks identified by the two different tools (DESeq2 and edgeR)
ta_analyze_1 <- dba.analyze(test, method=DBA_ALL_METHODS)
dba.plotVenn(ta_analyze_1,contrast=1,method=DBA_ALL_METHODS)
#查看差异分析的结果与导出为csv文件
ta.DB_1 <- dba.report(ta_analyze,contrast = 1)
write.csv(ta.DB_1, file = "diffbind_last.csv")
#以文件形式保存
comp1.edgeR <- dba.report(ta_analyze, method=DBA_EDGER, contrast = 1, th=1)
out <- as.data.frame(comp1.edgeR)
edge.bed <- out[ which(out$FDR < 0.05),]
write.csv(edge.bed, file="WT_vs_KD_edgeR.csv")
# Create bed files for each keeping only significant peaks (p < 0.05)
# EdgeR
out <- as.data.frame(comp1.edgeR)
edge.bed <- out[ which(out$FDR < 0.05), 
                 c("seqnames", "start", "end", "strand", "Fold")]
write.table(edge.bed, file="WT_vs_KD_edgeR.bed", sep="\t", quote=F, row.names=F, col.names=F)

##传统的上下调在找差异peak中称为“Gain”或“Loss”
#查看第1组的差异peak数量“Fold>0或<0”控制是“Gain”或“Loss”
sum(edge.bed$Fold>0,contrast=1)
sum(edge.bed$Fold<0,contrast=1)
#PCA
dba.plotPCA(ta_analyze) #这里可以使用所有样本进行PCA
dba.plotPCA(ta_analyze, contrast=1, label=DBA_FACTOR)#单独对分组1进行PCA

#MA plots
dba.plotMA(ta_analyze, contrast=1)

#Volcano plots
dba.plotVolcano(ta_analyze, contrast=1)

#Box plots
dba.plotBox(ta_analyze, contrast=1)

#Heatmap
hmap <- colorRampPalette(c("blue1", "grey", "red1"))(n = 10)
dba.plotHeatmap(ta_analyze,correlations=FALSE,scale="row",colScheme = hmap,margin=12)
dba.plotProfile(ta_analyze)

data <- data.frame(ta_analyze[["contrasts"]][[1]][["DESeq2"]][["de"]])
write.csv(data, file = "diffbind_all.csv")
