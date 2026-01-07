set.seed(123)

library(Seurat)
library(ggplot2)
library(dplyr)
library(monocle)

outpath <- "monocle2_dimension_ncenter_700"
if (!dir.exists(outpath)) dir.create(outpath)
print(outpath)
obj <- readRDS("subset_new_umap_dim18_spread2.5_mindist0.7_celltypes_250630_remove27_recluster_low_terminal_removeundefined_700.rds")
#expr_matrix <- as(as.matrix(obj@assays$RNA@counts), 'sparseMatrix')
expr_matrix <- GetAssayData(obj, assay = "RNA", slot = "counts")
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
              fullModelFormulaStr = "~human_celltypes_250702")
print(head(diff))
print(class(diff))
deg <- subset(diff, qval < 0.01)
deg <-deg[order(deg$qval,decreasing=F),]
ordergene <- rownames(deg) 
cds <- setOrderingFilter(cds, ordergene)

# 参数列表
ncenters=seq(from = 100, to = 1000, by = 100)
dimensions=seq(from = 2, to = 12, by = 1)

# 循环遍历所有参数组合
for (ncenter in ncenters) {
	for (dimension in dimensions) {

                 message(paste("Running with ncenter =", ncenter))

                 # 降维 + 排序
                 cds_subset <- reduceDimension(cds,
                                  max_components = dimension,
                                  auto_param_selection = FALSE,
                                  ncenter = ncenter,
                                  reduction_method = 'DDRTree')

                 cds_subset <- orderCells(cds_subset)
                 root_cells <- row.names(subset(pData(cds_subset), human_celltypes_250630 == "Rgl1-like"))
                 my_root_cell <- root_cells[1]
                 cds_subset@auxOrderingData[[cds_subset@dim_reduce_type]]$root_cell <- my_root_cell
                 cds_subset <- orderCells(cds_subset)
                 # 构建文件名
                 filename <- paste0(outpath, "/ncenter_", ncenter,"_dimension_",dimension,  ".png")

                 # 绘图并保存
                 png(filename, width = 9 * 100, height = 7 * 100)
                 plot1 <- plot_cell_trajectory(cds_subset, color_by = "human_celltypes_250702")
                 print(plot1)
                 dev.off()

                 # 可选：带 facet_wrap 的版本
                 plot2 <- plot_cell_trajectory(cds_subset, color_by = "human_celltypes_250702") +
                    facet_wrap("~human_celltypes_250702", ncol = 3)

                 filename2 <- paste0(outpath, "/ncenter_", ncenter,"_dimension_",dimension, "_facet.png")
                 png(filename2, width = 9 * 100, height = 7 * 100)
                 print(plot2)
                 dev.off()
         }
}


message("All combinations completed.")






