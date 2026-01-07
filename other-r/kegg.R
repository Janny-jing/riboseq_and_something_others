rm(list=ls())
# library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(clusterProfiler)
library(topGO)
library(Rgraphviz)
library(pathview)
library(openxlsx)
library(R.utils)
library(DOSE)
setwd("F:\\02.employment\\04.zhangxu\\03.kegg_2\\proviral")
path="F:\\02.employment\\04.zhangxu\\03.kegg_2\\proviral"
filename=dir()
for (i in 1:length(filename)){
  datapath <- paste(path,filename[i],sep='/')
  setwd(datapath)
  dataname=dir()
  data_1 <- read.table(dataname[1])
  data_2 <- bitr(data_1$V1, fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
  R.utils::setOption("clusterProfiler.download.method",'auto') 
  data_kegg <- enrichKEGG(gene = data_2$ENTREZID,
                          organism = 'hsa',
                          pAdjustMethod = "BH",
                          qvalueCutoff = 1,
                          pvalueCutoff = 1)
  pdf(file="kegg.pdf",width=9,height = 11)
  tryCatch(
    {print(dotplot(data_kegg, title="Enrichment_kegg",showCategory=25,label_format=100))},
    warning = function(w) {message("warning")},
    error = function(e) {message("error")}
  ) 
  dev.off()
  tiff(file="kegg.tiff",width=25,height = 20,units = "cm", compression = "lzw",bg="white",res=600)
  tryCatch(
    {print(dotplot(data_kegg, title="Enrichment_kegg",showCategory=25,label_format=100))},
    warning = function(w) {message("warning")},
    error = function(e) {message("error")}
  ) 
  dev.off()
  data_kegg_genename <- setReadable(data_kegg,OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  write.xlsx(data_kegg_genename@result,file = "data_kegg_genename.xlsx")
  # data_kegg_genename@result$Description[1:25] 
  for (i in 1:25){
    select_pathway <- data_kegg_genename@result$ID[i] 
    tryCatch(
      {print(pathview(gene.data     = data_2$ENTREZID,
                      pathway.id    = select_pathway,
                      species       = 'hsa' ,      
                      kegg.native   = T,
                      new.signature = F, 
                      limit         = list(gene=1),
                      bins = list(gene = 20)))},
      warning = function(w) {message("warning")},
      error = function(e) {message("error")}
      )} 
}

