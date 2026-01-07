library(Seurat)
library(tidyverse)
library(patchwork)
library(Matrix)
library(stringr)
setwd("/storage/liuxiaodongLab/jiangjing/Projects/YutingFu/PD_YutingFu/scvelo/data")
# Load the Seurat object
sobj <- readRDS("/storage/liuxiaodongLab/jiangjing/Projects/YutingFu/PD_YutingFu/scvelo/scanpy_UMAP_n_neighbors10_npcs23_mindist0.5_spread0.5_copy.rds")
table(sobj$batch)
# 1. 提取目标批次的细胞子集
target_batches <- c("P1A", "P1B")
sobj_filtered <- subset(sobj, subset = batch %in% target_batches)

# 2. 验证筛选结果
table(sobj_filtered$batch)  # 确认批次分布符合预期

# 3. 提取子集的计数矩阵
data_filtered <- LayerData(sobj, assay = "RNA", layer = "counts")

# 4. 获取特征和条形码信息
original_barcodes <- colnames(data_filtered)
orig_names <- sobj@meta.data$orig.cell.names
new_barcodes <- ifelse(is.na(orig_names), original_barcodes, as.character(orig_names))
colnames(data_filtered) <- new_barcodes
barcodes_filtered <- colnames(data_filtered)  # 子集细胞的条形码
features_filtered <- rownames(data_filtered)  # 基因名（与原始一致）

# 5. 创建输出目录
output_dir_filtered <- "./filtered_feature_bc_matrix/"
if (!dir.exists(output_dir_filtered)) {
  dir.create(output_dir_filtered, recursive = TRUE)
}

# 6. 导出barcodes.tsv.gz
write.table(
  barcodes_filtered,
  file = gzfile(file.path(output_dir_filtered, "barcodes.tsv.gz")),
  sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE
)

# 7. 导出features.tsv.gz (保持与原始一致)
features_df_filtered <- data.frame(
  feature_id = features_filtered,
  feature_name = features_filtered,
  feature_type = "Gene Expression"
)
write.table(
  features_df_filtered,
  file = gzfile(file.path(output_dir_filtered, "features.tsv.gz")),
  sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE
)

# 8. 导出稀疏矩阵matrix.mtx.gz
mtx_filtered <- as(data_filtered, "sparseMatrix")
temp_mtx_path <- tempfile()
Matrix::writeMM(mtx_filtered, temp_mtx_path)
system(sprintf("gzip -c %s > %s", 
               temp_mtx_path, 
               file.path(output_dir_filtered, "matrix.mtx.gz")))
