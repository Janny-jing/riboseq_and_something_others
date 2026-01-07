rm(list = ls())
setwd("Z:\\02.ATACseq_rnaseq\\01.A549_polyic_ATACseq_RNAseq\\01.A549_PIC_RNAseq")
filename <- dir()
data<-read.table(filename[5],header = T,row.names = 1)
data <- data[,6:8]
expr <- na.omit(data)
expr <- expr[apply(expr,1,max)>=2, ]
dup_gene<-data.frame(gene=expr$gene_name,mean=apply(expr[,-grep("gene_name",colnames(expr))],1,mean))
expr <-expr[order(dup_gene$gene,dup_gene$mean,decreasing = T),]
expr<- expr[!duplicated(expr$gene_name),]
rownames(expr) <- expr$gene_name
expr <- expr[,-grep("gene_name",colnames(expr))]
colnames(expr) <- c("A549_CTL_rnaseq","A549_PIC_rnaseq")
write.table(expr,file = "A549_RNASEQ_readcount.txt",sep = "\t",row.names = TRUE,col.names = TRUE)
