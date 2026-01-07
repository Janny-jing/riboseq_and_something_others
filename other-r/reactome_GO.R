# install.packages("D:\\pf\\R-4.4.1\\library\\reactome.db_1.88.0.tar.gz", repos = NULL, type = "source")
# BiocManager::install("ReactomePA")
# BiocManager::install("volcanoplot")

library(reactome.db)
library(ReactomePA)
library(openxlsx)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)

setwd("E:\\temp\\GO")
path="E:\\temp\\GO"
filename=dir()

for (i in 1:length(filename)){
  data_1 <- read.xlsx(filename[i],colNames = T)
  data_2 <- bitr(data_1$Gene, fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Mm.eg.db)
  data_ALL <- enrichGO(gene  = data_2$ENTREZID,
                       OrgDb         = org.Mm.eg.db,
                       ont           = "ALL" , 
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05,
                       readable      = TRUE)
  res <- enrichPathway(data_2$ENTREZID, 
                       organism = "mouse",
                       minGSSize = 10, maxGSSize = 200,
                       pvalueCutoff = 0.05, pAdjustMethod = "BH")
  res_genename <- setReadable(res,OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
  write.xlsx(data_ALL@result,file =file.path(path, paste0(gsub("\\.xlsx$", "", basename(filename[i])), "_GO.xlsx")))
  write.xlsx(res_genename@result,file =file.path(path, paste0(gsub("\\.xlsx$", "", basename(filename[i])), "_reactome.xlsx")))
  png(file=file.path(path, paste0(gsub("\\.xlsx$", "", basename(filename[i])), "_GO.png")), width=7000, height=4000, res=600)
  tryCatch(
    {print(dotplot(data_ALL, showCategory=20,title="EnrichmentGO_BP",label_format = 100))},
    warning = function(w) {message("warning")},
    error = function(e) {message("error")}
  ) 
  dev.off()
  png(file=file.path(path, paste0(gsub("\\.xlsx$", "", basename(filename[i])), "_GO2.png")), width=7000, height=4000, res=600)
  tryCatch(
    {print(barplot(data_ALL, showCategory=20,title="EnrichmentGO_BP",label_format = 100))},
    warning = function(w) {message("warning")},
    error = function(e) {message("error")}
  ) 
  dev.off()
  pdf(file=file.path(path, paste0(gsub("\\.xlsx$", "", basename(filename[i])), "_GO.pdf")),width=9,height = 11)
  tryCatch(
    {print(dotplot(data_ALL, showCategory=20,title="EnrichmentGO_BP",label_format = 100))},
    warning = function(w) {message("warning")},
    error = function(e) {message("error")}
  ) 
  dev.off()
  pdf(file=file.path(path, paste0(gsub("\\.xlsx$", "", basename(filename[i])), "_GO2.pdf")),width=9,height = 11)
  tryCatch(
    {print(barplot(data_ALL, showCategory=20,title="EnrichmentGO_BP",label_format = 100))},
    warning = function(w) {message("warning")},
    error = function(e) {message("error")}
  ) 
  dev.off()
  png(file=file.path(path, paste0(gsub("\\.xlsx$", "", basename(filename[i])), "_reactome.png")), width=7000, height=4000, res=600)
  tryCatch(
    {print(dotplot(res, showCategory=20,title="reactome pathway analysis",label_format = 100))},
    warning = function(w) {message("warning")},
    error = function(e) {message("error")}
  ) 
  dev.off()
  png(file=file.path(path, paste0(gsub("\\.xlsx$", "", basename(filename[i])), "_reactome2.png")), width=7000, height=4000, res=600)
  tryCatch(
    {print(barplot(res, showCategory=20,title="reactome pathway analysis",label_format = 100))},
    warning = function(w) {message("warning")},
    error = function(e) {message("error")}
  ) 
  dev.off()
  pdf(file=file.path(path, paste0(gsub("\\.xlsx$", "", basename(filename[i])), "_reactome.pdf")),width=9,height = 11)
  tryCatch(
    {print(dotplot(res, showCategory=20,title="reactome pathway analysis",label_format = 100))},
    warning = function(w) {message("warning")},
    error = function(e) {message("error")}
  ) 
  dev.off()
  pdf(file=file.path(path, paste0(gsub("\\.xlsx$", "", basename(filename[i])), "_reactome2.pdf")),width=9,height = 11)
  tryCatch(
    {print(barplot(res, showCategory=20,title="reactome pathway analysis",label_format = 100))},
    warning = function(w) {message("warning")},
    error = function(e) {message("error")}
  ) 
  dev.off()
}


