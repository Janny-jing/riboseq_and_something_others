rm(list = ls())
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
require(plyr)
library(edgeR)
setwd("G:\\temp\\deseq2")
filename <- dir()
data<-read.table(filename[1],header = T,row.names = 1)
expr <- na.omit(data)
expr <- expr[apply(expr,1,max)>=10, ]
dup_gene<-data.frame(gene=expr$gene_name,mean=apply(expr[,-grep("gene_name",colnames(expr))],1,mean))
expr <-expr[order(dup_gene$gene,dup_gene$mean,decreasing = T),]
expr<- expr[!duplicated(expr$gene_name),]
rownames(expr) <- expr$gene_name
expr <- expr[,-grep("gene_name",colnames(expr))]
colnames(expr) <- c("WT_R1","WT_R2","WT_R3","WT_R4","MDX_R1","MDX_R2","MDX_R3","MDX_R4","MDX_R5")
col.data = data.frame(Sample = rep(c("WT","mdx"), c(4,5)))
dds <- DESeqDataSetFromMatrix(countData = expr, colData = col.data, design = ~ Sample)
## Get normalized counts
dir.create("./correlation")
setwd("./correlation")
dds <- estimateSizeFactors(dds)
dds.cnt.nrm = counts(dds, normalized = TRUE)
write.table(dds.cnt.nrm, file="merged_count.normalized.txt", sep="\t", quote=FALSE, row.names = TRUE, col.names = TRUE)
## Correlation_heatmap
dds.cnt.nrm.l10 = log10(dds.cnt.nrm+1)
tiff(file="Correlation_heatmap.tiff",width=20,height = 20,units = "cm", compression = "lzw",bg="white",res=600)
pheatmap(cor(dds.cnt.nrm.l10, method = "p"), 
         color = colorRampPalette(c("light pink","white","purple"))(50),
         cluster_rows = FALSE, cluster_cols = FALSE,
         display_numbers = TRUE, number_format = "%.3f")
dev.off()
pdf(file="Correlationheatmap.pdf",width=5,height = 5)
pheatmap(cor(dds.cnt.nrm.l10, method = "p"), 
         color = colorRampPalette(c("light pink","white","purple"))(50),
         cluster_rows = FALSE, cluster_cols = FALSE,
         display_numbers = TRUE, number_format = "%.3f")

dev.off()



## PCA analysis
rld <- rlog(dds, blind = TRUE)
pca.rld.data <- plotPCA(rld, intgroup = c("Sample"),returnData = TRUE, ntop = 1000)
percent.rld.var = round(100*attr(pca.rld.data, "percentVar"))
pdf(file="PCA.pdf",width=5,height = 5)
ggplot(pca.rld.data, aes(PC1, PC2, color=Sample)) +
  theme_bw() + 
  # xlim(-15,15) + ylim(-3, 3) +
  geom_point(size = 4) +
  geom_text(aes(label=name),hjust=-0.5, vjust=0.4) +
  xlab(paste0("PC1: ", percent.rld.var[1], "% variance")) +
  ylab(paste0("PC2: ", percent.rld.var[2], "% variance"))
dev.off()
tiff(file="PCA.tiff",width=20,height = 20,units = "cm", compression = "lzw",bg="white",res=600)
ggplot(pca.rld.data, aes(PC1, PC2, color=Sample)) +
  theme_bw() + 
  # xlim(-15,15) + ylim(-3, 3) +
  geom_point(size = 4) +
  geom_text(aes(label=name),hjust=-0.5, vjust=0.4) +
  xlab(paste0("PC1: ", percent.rld.var[1], "% variance")) +
  ylab(paste0("PC2: ", percent.rld.var[2], "% variance"))
dev.off()
####### DE analysis  ########
dir.create("../deseq2")
setwd("../deseq2")
design(dds) = ~Sample
dds = DESeq(dds)
res  = results(dds, contrast = c("Sample", "mdx","WT"))
write.csv(res,file = "DEanalysis.csv")

# Get DEGs
up.lst <- row.names(res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange > 1 & res$padj < 0.05,])
dn.lst <- row.names(res[!is.na(res$log2FoldChange) & !is.na(res$padj) & res$log2FoldChange < -1 & res$padj < 0.05,])
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

# MA-plot
tiff(file="MA.tiff",width=20,height = 20,units = "cm", compression = "lzw",bg="white",res=600)
par(mar=c(4,4,1,1), mfrow=c(1,1))
plot(
  log10(res$baseMean+1), res$log2FoldChange, xlab="Average read counts", ylab="log2fold(MDX/WT)",
  cex=.1, pch=19, col=alpha('black',0.5), xaxt = "n", xlim=c(0,6), ylim=c(-10,10)
)
axis(1, at = 0:6, labels = 10^(0:6))
abline(h=0, lwd=3, col=alpha('red',0.5))
points(log10(res$baseMean[row.names(res) %in% up.lst]+1), res$log2FoldChange[row.names(res) %in% up.lst], cex=0.3, pch=19, col='red')
points(log10(res$baseMean[row.names(res) %in% dn.lst]+1), res$log2FoldChange[row.names(res) %in% dn.lst], cex=0.3, pch=19, col='green')
dev.off()

pdf(file="MA.pdf",width=5,height = 5)
par(mar=c(4,4,1,1), mfrow=c(1,1))
plot(
  log10(res$baseMean+1), res$log2FoldChange, xlab="Average read counts", ylab="log2fold(MDX/WT)",
  cex=.1, pch=19, col=alpha('black',0.5), xaxt = "n", xlim=c(0,6), ylim=c(-10,10)
)
axis(1, at = 0:6, labels = 10^(0:6))
abline(h=0, lwd=3, col=alpha('red',0.5))
points(log10(res$baseMean[row.names(res) %in% up.lst]+1), res$log2FoldChange[row.names(res) %in% up.lst], cex=0.3, pch=19, col='red')
points(log10(res$baseMean[row.names(res) %in% dn.lst]+1), res$log2FoldChange[row.names(res) %in% dn.lst], cex=0.3, pch=19, col='green')
dev.off()

## difgene_pheatmap
expr_cpm <- log2(cpm(expr)+1)
res.srt <- na.omit(res.srt)
res.srt$regulated <- "normal"
loc_up <- intersect(which(res.srt$log2FoldChange > 1),which(res.srt$padj < 0.05))
loc_down <- intersect(which(res.srt$log2FoldChange < -1),which(res.srt$padj < 0.05))
res.srt$regulated[loc_up] <- "up"
res.srt$regulated[loc_down] <- "down"
DEGSeq2_sigGene <- rownames(res.srt)[res.srt$regulated!="normal"]
DEGSeq2_sig <- res.srt[res.srt$regulated!="normal",]
dat1 <- expr_cpm[match(DEGSeq2_sigGene,rownames(expr_cpm)),]
x=DEGSeq2_sig$log2FoldChange 
names(x)=rownames(DEGSeq2_sig) 
x=na.omit(x)
cg=c(names(head(sort(x),50)),
     names(tail(sort(x),50)))
dat1 <- as.data.frame(dat1)
dat_p <- dat1[cg,]
dat_p=na.omit(dat_p)
pheatmap(dat_p,
         show_rownames = T,
         show_colnames = T,
         cluster_cols = F,
         cluster_rows=F,
         filename='de_pheatmap.tiff',
         fontsize_row=15, 
         fontsize_col=15,
         height=40,  
         width =15,
         scale = "row",
         angle_col=45, 
         color =colorRampPalette(c("blue", "#ffffff","red"))(100),
         clustering_distance_rows = 'euclidean',
         clustering_method = 'single',
         cellwidth = 25,
         cellheight = 25
)
pheatmap(dat_p,
         show_rownames = T,
         show_colnames = T,
         cluster_cols = F,
         cluster_rows=F,
         filename='de_pheatmap.pdf',
         fontsize_row=15, 
         fontsize_col=15,
         height=40,  
         width =20,
         scale = "row",
         angle_col=45, 
         color =colorRampPalette(c("blue", "#ffffff","red"))(100),
         clustering_distance_rows = 'euclidean',
         clustering_method = 'single',
         cellwidth = 25,
         cellheight = 25
)
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
pdf(file="allheatmap.pdf",width=5,height = 5)
pheatmap(
  log10(de.cnt.nrm+1), scale = "row", show_rownames = FALSE, annotation_row = de.cnt.ann,
  border_color = NA,
  annotation_colors = list(
    Expression=c("Up-regulated"="Red", "Down-regulated"="Green", "Unchanged"="Grey")))
dev.off()
tiff(file="all_heatmap.tiff",width=20,height = 20,units = "cm", compression = "lzw",bg="white",res=600)
pheatmap(
  log10(de.cnt.nrm+1), scale = "row", show_rownames = FALSE, annotation_row = de.cnt.ann,
  border_color = NA,
  annotation_colors = list(
    Expression=c("Up-regulated"="Red", "Down-regulated"="Green", "Unchanged"="Grey")))
dev.off()