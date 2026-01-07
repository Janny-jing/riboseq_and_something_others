library(RcisTarget)
library(Seurat)
library(tidyverse)
library(openxlsx)
library(DT)
library(visNetwork)
library(reshape2)
library(webshot)
setwd("E:\\temp\\rcistarget")
filename =dir()
data <- read.csv(filename[1])
data_list <- as.list(data)


#这里有两个motif注释数据可以选  根据自己分析对象进行选择
# mouse:
data(motifAnnotations_mgi)
# human:
# data(motifAnnotations_hgnc)  #选择人的

#这个数据是比较难的，这个数据要从一个数据库上下载，
#很大，每一个文件大概1G，经常会下载不完整，会导致读入不行，会报错
#方法一，用R直接下载并读取
# featherURL <- "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-tss-centered-10kb-7species.mc9nr.feather" 
# download.file(featherURL, destfile=basename(featherURL)) #这样就会直接下载到当前工作目录下
#读入
# motifRankings <- importRankings("hg19-tss-centered-10kb-7species.mc9nr.feather")      #这里假如数据下载不全读取R会奔溃并restart，不是电脑的配置问题

#方法二，读取我自己提前下好的数据，
motifRankings <- importRankings("D:/pf/R-4.4.1/library/mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather")  #我的数据存放在cisTarget_databases文件夹下，跟网上下载有啥不同呢   hg19与hg38就是数据的版本不同，10kb与500bp就是想要研究基因上下调控的范围   10kb更大  所以读入hg19或hg38都可以。 

motifEnrichmentTable_wGenes <- cisTarget(data_list, 
                                         motifRankings,
                                         motifAnnot=motifAnnotations)#一步搞定分析

motifEnrichmentTable_wGenes_wLogo <- addLogo(motifEnrichmentTable_wGenes)

resultsSubset <- motifEnrichmentTable_wGenes_wLogo[1:10,]    #


# dtable <- datatable(resultsSubset[,-c("enrichedGenes", "TF_lowConf"), with=FALSE], 
#                     escape = FALSE, # To show the logo
#                     filter="top", options=list(pageLength=5))
dtable2 <- datatable(resultsSubset[, with=FALSE], 
                    escape = FALSE, # To show the logo
                    filter="top",extensions = 'Buttons', options = list(dom='Bfrtip',buttons=c('copy', 'csv', 'excel', 'print', 'pdf')))


html <- "dtable.html"
saveWidget(dtable2, html)
webshot(html, "dtable.pdf")

# 1. Calculate AUC
motifs_AUC <- calcAUC(data_list, motifRankings)

# 2. Select significant motifs, add TF annotation & format as table
motifEnrichmentTable <- addMotifAnnotation(motifs_AUC,
                                           nesThreshold=3,
                                           motifAnnot=motifAnnotations)

## 3. Identify significant genes for each motif
# (i.e. genes from the gene set in the top of the ranking)
# Note: Method 'iCisTarget' instead of 'aprox' is more accurate, but slower
motifEnrichmentTable_wGenes <- addSignificantGenes(motifEnrichmentTable, 
                                                   geneSets=data_list,
                                                   rankings=motifRankings, 
                                                   nCores=1,
                                                   method="aprox")

signifMotifNames <- motifEnrichmentTable$motif[1:10]

incidenceMatrix <- getSignificantGenes(data_list$Macrophage.C3, 
                                       motifRankings,
                                       signifRankingNames=signifMotifNames,
                                       plotCurve=TRUE, maxRank=5000, 
                                       genesFormat="incidMatrix",
                                       method="aprox")$incidMatrix

edges <- melt(incidenceMatrix)
edges <- edges[which(edges[,3]==1),1:2]
colnames(edges) <- c("from","to")
motifs <- unique(as.character(edges[,1]))
genes <- unique(as.character(edges[,2]))
nodes <- data.frame(id=c(motifs, genes),   
                    label=c(motifs, genes),    
                    title=c(motifs, genes), # tooltip 
                    shape=c(rep("diamond", length(motifs)), rep("elypse", length(genes))),
                    color=c(rep("pink", length(motifs)), rep("lightblue", length(genes))))
network <- visNetwork(nodes, edges,width = "100%",height = "1000px") %>% visOptions(highlightNearest = TRUE, 
                                         nodesIdSelection = TRUE) #这种网络图都是自定义的，也可以画TF与其互作基因的网络图
visSave(network, file = "network.html",background = "grey")
