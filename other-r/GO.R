rm(list = ls())
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("pathview")
library(org.Hs.eg.db)
library(clusterProfiler)
library(topGO)
library(Rgraphviz)
library(pathview)
library(openxlsx)

setwd("E:\\02.employment\\04.zhangxu\\04.orf3a_go")
path="E:\\02.employment\\04.zhangxu\\04.orf3a_go"
filename=dir()
for (i in 1:length(filename)){
  datapath <- paste(path,filename[i],sep='/')
  setwd(datapath)
  dataname=dir()
  data_1 <- read.table(dataname[1])
  data_2 <- bitr(data_1$V1, fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
  data_ALL <- enrichGO(gene  = data_2$ENTREZID,
                           OrgDb         = org.Hs.eg.db,
                           ont           = "ALL" , 
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 0.05,
                           qvalueCutoff  = 0.05,
                           readable      = TRUE)
  
  write.xlsx(data_ALL@result,file ="GO_ALL.xlsx")
  png(file="GO_BP_1.png", width=7000, height=4000, res=600)
  tryCatch(
    {print(dotplot(data_ALL, showCategory=20,title="EnrichmentGO_BP",label_format = 100))},
    warning = function(w) {message("warning")},
    error = function(e) {message("error")}
  ) 
  dev.off()
  png(file="GO_BP_2.png", width=7000, height=4000, res=600)
  tryCatch(
    {print(barplot(data_ALL, showCategory=20,title="EnrichmentGO_BP",label_format = 100))},
    warning = function(w) {message("warning")},
    error = function(e) {message("error")}
  ) 
  dev.off()
  pdf(file="GO_1.pdf",width=9,height = 11)
  tryCatch(
    {print(dotplot(data_ALL, showCategory=20,title="EnrichmentGO_BP",label_format = 100))},
    warning = function(w) {message("warning")},
    error = function(e) {message("error")}
  ) 
  dev.off()
  pdf(file="GO_2.pdf",width=9,height = 11)
  tryCatch(
    {print(barplot(data_ALL, showCategory=20,title="EnrichmentGO_BP",label_format = 100))},
    warning = function(w) {message("warning")},
    error = function(e) {message("error")}
  ) 
  dev.off()
}


################
rm(list = ls())
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("pathview")
library(org.Mm.eg.db)
library(clusterProfiler)
library(topGO)
library(Rgraphviz)
library(pathview)
library(openxlsx)

setwd("E:\\temp\\GO")
path="E:\\temp\\GO"
filename=dir()
for (i in 1:length(filename)){
  datapath <- paste(path,filename[i],sep='\\')
  data_1 <- read.xlsx(datapath,colNames = T)
  data_2 <- bitr(data_1$Gene, fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Mm.eg.db)
  data_ALL <- enrichGO(gene  = data_2$ENTREZID,
                       OrgDb         = org.Mm.eg.db,
                       ont           = "ALL" , 
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05,
                       readable      = TRUE)
  
  write.xlsx(data_ALL@result,file =file.path(path, paste0(gsub("\\.xlsx$", "", basename(datapath)), "_GO.xlsx")))
  png(file=file.path(path, paste0(gsub("\\.xlsx$", "", basename(datapath)), "_GO.png")), width=7000, height=4000, res=600)
  tryCatch(
    {print(dotplot(data_ALL, showCategory=20,title="EnrichmentGO_BP",label_format = 100))},
    warning = function(w) {message("warning")},
    error = function(e) {message("error")}
  ) 
  dev.off()
  png(file=file.path(path, paste0(gsub("\\.xlsx$", "", basename(datapath)), "_GO2.png")), width=7000, height=4000, res=600)
  tryCatch(
    {print(barplot(data_ALL, showCategory=20,title="EnrichmentGO_BP",label_format = 100))},
    warning = function(w) {message("warning")},
    error = function(e) {message("error")}
  ) 
  dev.off()
  pdf(file=file.path(path, paste0(gsub("\\.xlsx$", "", basename(datapath)), "_GO.pdf")),width=9,height = 11)
  tryCatch(
    {print(dotplot(data_ALL, showCategory=20,title="EnrichmentGO_BP",label_format = 100))},
    warning = function(w) {message("warning")},
    error = function(e) {message("error")}
  ) 
  dev.off()
  pdf(file=file.path(path, paste0(gsub("\\.xlsx$", "", basename(datapath)), "_GO2.pdf")),width=9,height = 11)
  tryCatch(
    {print(barplot(data_ALL, showCategory=20,title="EnrichmentGO_BP",label_format = 100))},
    warning = function(w) {message("warning")},
    error = function(e) {message("error")}
  ) 
  dev.off()
}



