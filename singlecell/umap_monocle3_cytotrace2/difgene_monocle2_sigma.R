set.seed(123)

library(Seurat)
library(ggplot2)
library(dplyr)
library(monocle)

outpath <- "dofgene_monocle2_results"
if (!dir.exists(outpath)) dir.create(outpath)
print(outpath)
obj <- readRDS("scanpy_UMAP_n_neighbors10_npcs23_mindist0.5_spread0.5_copy.rds")
expr_matrix <- as(as.matrix(obj@assays$RNA@counts), 'sparseMatrix')
p_data <- obj@meta.data
f_data <- data.frame(gene_short_name = row.names(obj),row.names = row.names(obj))
pd <- new('AnnotatedDataFrame', data = p_data) 
fd <- new('AnnotatedDataFrame', data = f_data)
cds <- newCellDataSet(expr_matrix,phenoData = pd,
                       featureData = fd,
                       lowerDetectionLimit = 0.5,
                       expressionFamily = negbinomial.size())
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

Idents(obj) <- obj$human_celltypes_250418
all.markers <- FindAllMarkers(
  object = obj,
  only.pos = TRUE,
  min.diff.pct = 0.25,
  logfc.threshold = 0.25
)

# 每个 cluster 取 avg_log2FC 最高的前 300 个基因
markers_by_cluster <- all.markers %>%
  group_by(cluster) %>%
  slice_max(n = 300, order_by = avg_log2FC) %>%
  pull(gene) %>%
  unique()

# 取 top 400 上调基因
high_vs_low <- FindMarkers(
  object = obj,
  ident.1 = "DA progenitor-CHRNB3 high",
  ident.2 = "DA progenitor-CHRNB3 low",
  only.pos = TRUE,
  logfc.threshold = 0.25
)
low_vs_high <- FindMarkers(
  object = obj,
  ident.1 = "DA progenitor-CHRNB3 low",
  ident.2 = "DA progenitor-CHRNB3 high",
  only.pos = TRUE,
  logfc.threshold = 0.25
)

degs_low_vs_high <- head(arrange(low_vs_high, desc(avg_log2FC)), 400) %>%
  rownames() %>%
  as.vector()

# 取 top 400 上调基因
degs_high_vs_low <- head(arrange(high_vs_low, desc(avg_log2FC)), 400) %>%
  rownames() %>%
  as.vector()
final_gene_list <- unique(c(markers_by_cluster, degs_high_vs_low,degs_low_vs_high))

print(paste("Total unique genes:", length(final_gene_list)))

diff <- differentialGeneTest(cds[final_gene_list,],
              fullModelFormulaStr = "~human_celltypes_250418")

print(head(diff))
print(class(diff))
deg <- subset(diff, qval < 0.01)
deg <-deg[order(deg$qval,decreasing=F),]
ordergene <- rownames(deg) 
cds <- setOrderingFilter(cds, ordergene)

# 参数列表
sigmas <- seq(from = 0.1, to = 1, by = 0.1)

# 循环遍历所有参数组合
for (sigma in sigmas) {
    message(paste("sigma =", sigma))

    # 降维 + 排序
    cds_subset <- reduceDimension(cds,
                                  max_components = 2,
                                  auto_param_selection = FALSE,
				  sigma = sigma,
                                  reduction_method = 'DDRTree')

    cds_subset <- orderCells(cds_subset)
    root_cells <- row.names(subset(pData(cds_subset), human_celltypes_250418 == "Rgl1-like"))
    my_root_cell <- root_cells[1]
    cds_subset@auxOrderingData[[cds_subset@dim_reduce_type]]$root_cell <- my_root_cell
    cds_subset <- orderCells(cds_subset)
    # 构建文件名
    filename <- paste0(outpath, "/sigma_", sigma, ".png")

    # 绘图并保存
    png(filename, width = 9 * 100, height = 7 * 100, units = "px")
    plot1 <- plot_cell_trajectory(cds_subset, color_by = "human_celltypes_250418")
    print(plot1)
    dev.off()

    # 可选：带 facet_wrap 的版本
    plot2 <- plot_cell_trajectory(cds_subset, color_by = "human_celltypes_250418") +
      facet_wrap("~human_celltypes_250418", ncol = 3)

    filename2 <- paste0(outpath, "/sigma_", sigma, "_facet.png")
    png(filename2, width = 9 * 100, height = 7 * 100, units = "px")
    print(plot2)
    dev.off()
}


message("All combinations completed.")






