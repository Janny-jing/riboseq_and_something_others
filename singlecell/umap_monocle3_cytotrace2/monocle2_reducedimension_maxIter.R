set.seed(123)

library(Seurat)
library(ggplot2)
library(dplyr)
library(monocle)

outpath <- "monocle2_results"
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

disp_table <- dispersionTable(cds)
if (!"dispersion_norm" %in% names(disp_table)) {
  disp_table$dispersion_norm <- disp_table$dispersion_empirical / disp_table$dispersion_fit
}
expressed_genes <- disp_table %>% arrange(desc(dispersion_norm)) %>% slice(1:2000) %>% pull(gene_id)
diff <- differentialGeneTest(cds[expressed_genes,],
              fullModelFormulaStr = "~human_celltypes_250418")
print(head(diff))
print(class(diff))
deg <- subset(diff, qval < 0.01)
deg <-deg[order(deg$qval,decreasing=F),]
ordergene <- rownames(deg) 
cds <- setOrderingFilter(cds, ordergene)

# 参数列表
maxIters <- seq(from = 10, to = 100, by = 10)

# 循环遍历所有参数组合
for (maxIter in maxIters) {
    message(paste("maxIter =", maxIter))

    # 降维 + 排序
    cds_subset <- reduceDimension(cds,
                                  max_components = 2,
                                  auto_param_selection = FALSE,
				  maxIter = maxIter,
                                  reduction_method = 'DDRTree')

    cds_subset <- orderCells(cds_subset)
    root_cells <- row.names(subset(pData(cds_subset), human_celltypes_250418 == "Rgl1-like"))
    my_root_cell <- root_cells[1]
    cds_subset@auxOrderingData[[cds_subset@dim_reduce_type]]$root_cell <- my_root_cell
    cds_subset <- orderCells(cds_subset)
    # 构建文件名
    filename <- paste0(outpath, "/maxIter_", maxIter, ".pdf")

    # 绘图并保存
    png(filename, width = 9 * 100, height = 7 * 100, units = "px")
    plot1 <- plot_cell_trajectory(cds_subset, color_by = "human_celltypes_250418")
    print(plot1)
    dev.off()

    # 可选：带 facet_wrap 的版本
    plot2 <- plot_cell_trajectory(cds_subset, color_by = "human_celltypes_250418") +
      facet_wrap("~human_celltypes_250418", ncol = 3)

    filename2 <- paste0(outpath, "/maxIter_", maxIter, "_facet.pdf")
    png(filename2, width = 9 * 100, height = 7 * 100, units = "px")
    print(plot2)
    dev.off()
}


message("All combinations completed.")






