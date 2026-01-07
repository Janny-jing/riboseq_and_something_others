rm(list = ls())
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
require(plyr)
library(edgeR)
setwd("Z:\\16.up_down_go_reactome\\01.WT_IFNR")
data<-read.table("00.raw_count.txt",header = T,row.names = 1)
expr <- na.omit(data)
expr <- expr[apply(expr,1,max)>=10, ]
dup_gene<-data.frame(gene=expr$gene_name,mean=apply(expr[,-grep("gene_name",colnames(expr))],1,mean))
expr <-expr[order(dup_gene$gene,dup_gene$mean,decreasing = T),]
expr<- expr[!duplicated(expr$gene_name),]
rownames(expr) <- expr$gene_name
expr <- expr[,-grep("gene_name",colnames(expr))]
# colnames(expr) <- c("WT_R1","WT_R2","WT_R3","WT_R4","mdx_R1","mdx_R2","mdx_R3","mdx_R4","mdx_R5")
colnames(expr) <- c("WT_R1","WT_R2","WT_R3","WT_R4","WT_R5","IFN_R1","IFN_R2","IFN_R3","IFN_R4","IFN_R5","IFN_R6","IFN_R7")
col.data = data.frame(Sample = rep(c("WT","IFN"), c(5,7)))
dds <- DESeqDataSetFromMatrix(countData = expr, colData = col.data, design = ~ Sample)
## Get normalized counts
dir.create("./correlation")
setwd("./correlation")
dds <- estimateSizeFactors(dds)
dds.cnt.nrm = counts(dds, normalized = TRUE)
write.table(dds.cnt.nrm, file="merged_count.normalized.txt", sep="\t", quote=FALSE, row.names = TRUE, col.names = TRUE)
## Correlation_heatmap
dds.cnt.nrm.l10 = log10(dds.cnt.nrm+1)

####### DE analysis  ########
dir.create("../deseq2")
setwd("./heatmap_volcano")
design(dds) = ~Sample
dds = DESeq(dds)
res  = results(dds, contrast = c("Sample", "IFN","WT"))
write.csv(res,file = "DEanalysis.csv")

# Get DEGs
up.lst <- row.names(res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange > 1 & res$padj < 0.05 & res$baseMean >100,])
dn.lst <- row.names(res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange < -1 & res$padj < 0.05 & res$baseMean >100,])
length(up.lst)
length(dn.lst)
write.table(sub(":.*$","",up.lst), file="DEGs_Up.ENSEMBL_ID.txt", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
write.table(sub(":.*$","",dn.lst), file="DEGs_Down.ENSEMBL_ID.txt", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
# Sort by p-value
res.srt = res[order(res$pvalue),]
write.csv(res.srt, file = "DE_sorted.csv")
write.csv(res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange > 1 & res$padj < 0.05,], file = "DEGs_Up_sorted.csv")
write.csv(res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange < -1 & res$padj < 0.05,], file = "DEGs_down_sorted.csv")
de.lst = unique(c(up.lst,dn.lst))

# DEG Heatmap
de.cnt.nrm = dds.cnt.nrm[row.names(dds.cnt.nrm) %in% de.lst, ]
de.cnt.ann = data.frame(
  row.names = de.lst,
  Expression = factor(
    rep("Unchanged", length(de.lst)),
    levels = c("Up-regulated", "Down-regulated", "Unchanged")
  )
)
de.cnt.ann$Expression[row.names(de.cnt.ann) %in% up.lst] = "Up-regulated"
de.cnt.ann$Expression[row.names(de.cnt.ann) %in% dn.lst] = "Down-regulated"
de.cnt.ann$Expression = factor(de.cnt.ann$Expression, levels = c("Up-regulated", "Down-regulated", "Unchanged"))

##改变上下调位置顺序
de.cnt.nrm_log <- log10(de.cnt.nrm+1)
matrix2 <- round(t(apply(de.cnt.nrm_log,1, scale)),2)
colnames(matrix2) <- colnames(de.cnt.nrm_log)
exprTable <- as.data.frame(matrix2)#若是镜像列名需要转置matrix,行不需要
row_dist = dist(exprTable,method = "euclidean")
hclust_1 <- hclust(row_dist)
manual_order = rev(rownames(de.cnt.nrm_log))
dend = reorder(as.dendrogram(hclust_1), wts=order(match(manual_order, rownames(exprTable))))
row_cluster <- as.hclust(dend)


pdf(file="allheatmap.pdf",width=7,height = 5)
pheatmap(
  de.cnt.nrm_log, show_rownames = FALSE, annotation_row = de.cnt.ann,scale = "row",
  border_color = NA,cluster_cols = F,,cluster_rows = row_cluster,
  annotation_colors = list(
    Expression=c("Up-regulated"="Red", "Down-regulated"="Blue", "Unchanged"="White")))
dev.off()
tiff(file="all_heatmap.tiff",width=20,height = 20,units = "cm", compression = "lzw",bg="white",res=600)
pheatmap(
  de.cnt.nrm_log, show_rownames = FALSE, annotation_row = de.cnt.ann,scale = "row",
  border_color = NA,cluster_cols = F,cluster_rows = row_cluster,
  annotation_colors = list(
    Expression=c("Up-regulated"="Red", "Down-regulated"="Blue", "Unchanged"="White")))
dev.off()


