rm(list=ls())
setwd("E:\\02.employment\\03.tRNA\\11.BMDM_control_IL4_trnaanno")
library(ChIPseeker)
library(org.Mm.eg.db)
library(edgeR)
library(openxlsx)
filename <- dir()
data<-read.table(filename[5],header = F,row.names = 1)
colnames(data) <- c("CTR","LPS")
group <- factor(c(rep("CTR",1),rep("LPS",1)))
y <- DGEList(counts = data, genes = rownames(data), group = group)
#过滤
keep <- rowSums(cpm(y)>1) >= 1
y <- y[keep,,keep.lib.sizes=FALSE]
##TMM 标准化
y <- calcNormFactors(y)
y$samples
###推测离散度，若样本是人，设置bcv = 0.4，模式生物设置0.1(此处是根据有相关教程进行设置，也可以你根据你的结果自行设置，最终得到你的理想值)
bcv <- 0.1
et <- exactTest(y, dispersion=bcv^2)
topTags(et)

#矫正P值
DE <- et$table
res <- DE
head(DE)
res$FDRP <- p.adjust(res$PValue,method = "fdr",n=length(res$PValue))
head(res)

## 筛选出差异基因
## 筛选标准 FDR < 0.05，|logFC| > 1
diffsig <- res[(res$FDRP < 0.05 & abs(res$logFC) > 1),]
dim(diffsig)
## 新增一列，标显著性
res$Group <- "Not"
######
res$Group[which((res$FDRP < 0.05) & (res$logFC > 1))] = "Up"
res$Group[which((res$FDRP < 0.05) & (res$logFC < -1))] = "Down"
#### 查看DE数目
table(res$Group)
write.csv(res, file = "DE_sorted.csv")

###图形展示检验结果

library(ggplot2)

pdf('volcano.pdf',width = 7,height = 6.5)  ## 输出文件
numberup=res[res$Group=="Up",]
numberdn=res[res$Group=="Down",]
text <- paste("DE genesup = ",nrow(numberup),"\n","DE genesdn = ",nrow(numberdn),"\n","|logFC|>=0,p.adjust<0.05")
y_max <- (ceiling(max(-log10(res$PValue))/10)+ 1) *10
x_max <- ceiling(max(abs(res$logFC))) + 2
ggplot(res, aes(logFC, -log10(PValue)))+
  geom_point(aes(col=Group))+  #alpha=0.1, size=0.1
  scale_color_manual(values=c("#5069dd","#4f4f4f","#cc0001"))+
  labs(title = " ")+
  geom_vline(xintercept=c(-1,1), colour="grey", linetype="dashed")+
  geom_hline(yintercept = -log10(0.05),colour="grey", linetype="dashed")+
  labs(x="log2(FoldChange)",y="-log10(PValue)")+
  annotate("text", x=-x_max+8, y=y_max-15, label=text) +
  theme_classic()+
  theme(plot.title = element_text(size=15,hjust = 0.5),
        text = element_text(size = 14))+
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 


dev.off()

# library(org.Hs.eg.db)
# library(TxDb.Mmusculus.UCSC.mm10.knownGene)
# BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
library(TxDb.Mmusculus.UCSC.mm39.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm39.knownGene
data <- readPeakFile("gene_dn.txt")

data_anno <- annotatePeak(data, tssRegion=c(-1000, 1000), TxDb=txdb,annoDb="org.Mm.eg.db")

write.table(data_anno, file = "promoter_gene_dn.txt",sep = '\t', quote = FALSE, row.names = FALSE)
