library(sceasy)
library(Seurat)
setwd("/storage/liuxiaodongLab/jiangjing/Projects/YutingFu/PD_YutingFu/umap_monocle3_cytoTRACE2")
# 转换 H5AD 为 Seurat 对象
expression_matrix <- read.csv("expression_matrix.csv", row.names=1)
cell_metadata <- read.csv("cell_metadata.csv", row.names=1)
umap_coords <- read.csv("umap_coordinates.csv", row.names=1)

# 创建 SingleCellExperiment 对象
library(SingleCellExperiment)
sce <- SingleCellExperiment(
  assays = list(counts = t(as.matrix(expression_matrix))),
  colData = cell_metadata,
  reducedDims = list(UMAP = as.matrix(umap_coords))
)

###CytoTRACE2
library(CytoTRACE2)

# 输入表达矩阵（需 raw counts）
counts_matrix <- as.matrix(expression_matrix)

# 运行 CytoTRACE2
results <- cytotrace2(counts_matrix)

# 整合 UMAP 坐标和细胞类型
results$metadata <- cbind(umap_coords, cell_metadata)
# extract annotation data
annotation <- cell_metadata$human_celltypes_250416

# generate prediction and phenotype association plots with plotData function
plots <- plotData(cytotrace2_result = results, 
                  annotation = annotation,
                  expression_data = counts_matrix)
plots$CytoTRACE2_UMAP
ggsave("CytoTRACE2_UMAP.pdf",  width = 9, height = 7, dpi = 300)
plots$CytoTRACE2_Potency_UMAP
ggsave("CytoTRACE2_Potency_UMAP.pdf",  width = 9, height = 7, dpi = 300)
plots$CytoTRACE2_Relative_UMAP
ggsave("CytoTRACE2_Relative_UMAP.pdf",  width = 9, height = 7, dpi = 300)
plots$Phenotype_UMAP
ggsave("Phenotype_UMAP.pdf",  width = 9, height = 7, dpi = 300)
plots$CytoTRACE2_Boxplot_byPheno
ggsave("CytoTRACE2_Boxplot_byPheno.pdf",  width = 9, height = 7, dpi = 300)

#######monocle3
library(monocle3)
# 从 Seurat/SingleCellExperiment 转换
cds <- as.cell_data_set(sce)  # 或 cds <- new_cell_data_set(...)

# 直接指定 UMAP 坐标（避免重新计算）
reducedDims(cds)$UMAP <- as.matrix(umap_coords)

# 设置分析参数（需与 Scanpy 一致）
cds <- cluster_cells(cds, 
                     reduction_method = "UMAP",
                     k = 10,  # 对应 n_neighbors
                     num_dim = 23,  # 对应 n_pcs
                     verbose = TRUE)

# 轨迹推断
cds <- learn_graph(cds, use_partition = FALSE)

