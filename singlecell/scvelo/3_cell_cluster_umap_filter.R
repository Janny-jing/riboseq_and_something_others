rm(list = ls())
library(Seurat)
library(magrittr)
library(dplyr)
setwd("/storage/liuxiaodongLab/jiangjing/Projects/YutingFu/PD_YutingFu/scvelo")
data1 <- readRDS("/storage/liuxiaodongLab/jiangjing/Projects/YutingFu/PD_YutingFu/scvelo/new_umap_dim18_spread2.5_mindist0.7_celltypes_250418_orig.cellnames.rds")
data2 <- readRDS("/storage/liuxiaodongLab/jiangjing/Projects/YutingFu/PD_YutingFu/scvelo/scanpy_UMAP_n_neighbors10_npcs23_mindist0.5_spread0.5_copy.rds")
data3 <- readRDS("/storage/liuxiaodongLab/jiangjing/Projects/YutingFu/PD_YutingFu/scvelo/p1ab.rds")
data2@meta.data$orig.cell.names <- data1@meta.data[rownames(data2@meta.data), "orig.cell.names"]
# 处理以-A结尾的情况
data3@meta.data$orig.cell.names <- gsub(
  pattern = "^(.*)-A$", 
  replacement = "P1A_\\1-1",
  x = data3@meta.data$orig.cell.names
)

# 处理以-B结尾的情况
data3@meta.data$orig.cell.names <- gsub(
  pattern = "^(.*)-B$", 
  replacement = "P1B_\\1-1",
  x = data3@meta.data$orig.cell.names
)

#write.csv(Embeddings(data, reduction = "umap"), file = "cell_embeddings.csv")
#write.csv(data@meta.data[, 'human_celltypes_250416', drop = FALSE], file = "cell_clusters.csv")
head(data@meta.data)
table(data$batch)

###提取P1A，P1B
umap_1 <- Embeddings(data, reduction = "umap") 
umap_2 <- umap_1[grepl("^P1[AB]", rownames(umap_1)), ]
write.csv(umap_2, file = "cell_embeddings_P1A_P1B.csv")
cluster_1 <- as.data.frame(data@meta.data[, 'human_celltypes_250416', drop = FALSE])
cluster_2 <- cluster_1[grepl("^P1[AB]", rownames(cluster_1)), ,drop = FALSE]
write.csv(cluster_2, file = "cell_clusters_P1A_P1B.csv")


##修改barcode
if (sum(is.na(data2@meta.data$orig.cell.names)) > 0) {
  warning("orig.cell.names 包含 NA 值，需修复！")
}

# 若存在 NA，用原始行名（如细胞 barcode）填充 NA 的 orig.cell.names
row_names <- ifelse(
  is.na(data2@meta.data$orig.cell.names),
  rownames(data2@meta.data),  # 原始行名（如细胞 barcode）
  data2@meta.data$orig.cell.names
)
# UMAP embeddings
umap_embeddings <- Embeddings(data2, reduction = "umap")
rownames(umap_embeddings) <- row_names
write.csv(umap_embeddings, "cell_embeddings.csv")

# meta.data
cell_clusters <- data2@meta.data[, 'human_celltypes_250418', drop = FALSE]
rownames(cell_clusters) <- row_names
write.csv(cell_clusters, "cell_clusters.csv")