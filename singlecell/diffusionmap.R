rm(list=ls())
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(slingshot)
library(SingleCellExperiment)

setwd("/storage/liuxiaodongLab/jiangjing/Projects/YutingFu/PD_YutingFu/slingshot")
data <- readRDS("new_umap_dim18_spread2.5_mindist0.7_celltypes_250416.rds")
data.sce <- as.SingleCellExperiment(data)
sce_slingshot1 <- slingshot(data.sce,      #输入单细胞对象
                            reducedDim = 'UMAP',  #降维方式
                            clusterLabels = data.sce$human_celltypes_250416,  #cell类型
                            start.clus = 'Rgl1-like',       #轨迹起点,也可以不定义
                            approx_points = 150)
SlingshotDataSet(sce_slingshot1) 
