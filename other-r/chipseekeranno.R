rm(list=ls())
setwd("G:\\temp\\17.clipper\\hg19_amylee")
library(ChIPseeker)
# library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(clusterProfiler)
library(openxlsx)
# library(TxDb.Mmusculus.UCSC.mm10.knownGene)
# BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
columns(txdb)
data1 <- readPeakFile("elF3a_peaks.narrowPeak")
data2 <- readPeakFile("elF3g_peaks.narrowPeak")

peaks <- list(elF3a=data1, elF3g=data2)
peakAnnoList <- lapply(peaks, annotatePeak, TxDb=txdb,tssRegion=c(-100, 0), 
                       verbose=FALSE,addFlankGeneInfo=TRUE, flankDistance=5000,annoDb="org.Hs.eg.db")
# data_anno <- annotatePeak(data1, tssRegion=c(-500, 0), TxDb=txdb,annoDb="org.Mm.eg.db")
plotAnnoPie(peakAnnoList)

plotAnnoBar(peakAnnoList)


write.xlsx(peakAnnoList$elF3g, file = "elF3g_anno.xlsx")
data_anno_promoter <- subset(data_anno, annotation=="Promoter")
data_anno_promoter <- as.data.frame(data_anno_promoter)
data_anno_symbolid <- as.character(data_anno_promoter$SYMBOL)
geneid <- bitr(data_anno_symbolid, 
               fromType="SYMBOL",
               toType="ENTREZID", 
               OrgDb = "org.Hs.eg.db")
ego <- enrichGO(gene=geneid$ENTREZID, 
                OrgDb = org.Hs.eg.db, 
                ont = "ALL", 
                pAdjustMethod = "BH", 
                qvalueCutoff = 0.05, 
                readable = TRUE)
write.xlsx(ego@result,file ="go_all.xlsx")

ego_MF <- enrichGO(gene=geneid$ENTREZID, 
                   OrgDb = org.Hs.eg.db, 
                   ont = "MF", 
                   pAdjustMethod = "BH", 
                   qvalueCutoff = 0.05, 
                   readable = TRUE)
ego_CC <- enrichGO(gene=geneid$ENTREZID, 
                   OrgDb = org.Hs.eg.db, 
                   ont = "CC", 
                   pAdjustMethod = "BH", 
                   qvalueCutoff = 0.05, 
                   readable = TRUE)
ego_BP <- enrichGO(gene=geneid$ENTREZID, 
                   OrgDb = org.Hs.eg.db, 
                   ont = "BP", 
                   pAdjustMethod = "BH", 
                   qvalueCutoff = 0.05, 
                   readable = TRUE)


pdf(file="Go BP.pdf",width=10,height = 11)
tiff(file="Go BP.tiff",width=25,height = 20,units = "cm", compression = "lzw",bg="white",res=600)
barplot(ego_BP,showCategory=15,color="p.adjust",title="GO_BP Enrich",label_format=100)
dev.off()

png(file="ego_BP.png", width=7000, height=4000, res=600)
barplot(ego_BP, showCategory=15,title="EnrichmentGO_BP",label_format=100)
dev.off()
