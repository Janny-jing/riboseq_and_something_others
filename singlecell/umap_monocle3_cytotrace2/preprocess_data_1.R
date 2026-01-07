library(Seurat)
library(Matrix)
library(dplyr)

sobj <- readRDS("sampled_merge_data_rpca.rds")
sobj <- FindNeighbors(sobj, reduction = "integrated_rpca", dims = 1:30)
sobj <- FindClusters(sobj)
sobj <- RunUMAP(sobj, dims = 1:30, reduction = "integrated_rpca")
saveRDS(sobj, "sampled_merge_data_rpca_umap.rds")
