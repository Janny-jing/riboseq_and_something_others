set.seed(123)

library(Seurat)
library(ggplot2)
library(dplyr)
library(monocle)

outpath <- "monocle2_results"
if (!dir.exists(outpath)) dir.create(outpath)
print(outpath)
# 如果输出文件夹不存在，则进行创建
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

# 取 top 400 上调基因
degs_high_vs_low <- head(arrange(high_vs_low, desc(avg_log2FC)), 400) %>%
  rownames() %>%
  as.vector()
final_gene_list <- unique(c(markers_by_cluster, degs_high_vs_low))
print(paste("Total unique genes:", length(final_gene_list)))

diff <- differentialGeneTest(cds[final_gene_list,],
              fullModelFormulaStr = "~human_celltypes_250418")
print(head(diff))
print(class(diff))
deg <- subset(diff, qval < 0.01)
deg <-deg[order(deg$qval,decreasing=F),]
write.table(deg,file="celltype_deg_monocle.xls",col.names=T,row.names=F,sep="\t",quote=F)
ordergene <- rownames(deg) 
cds <- setOrderingFilter(cds, ordergene)
cds <- reduceDimension(cds, max_components = 2,
    method = 'DDRTree')

cds <- orderCells(cds)
unique(pData(cds)$human_celltypes_250418)
root_cells <- row.names(subset(pData(cds), human_celltypes_250418 == "Rgl1-like"))
my_root_cell <- root_cells[1]
cds@auxOrderingData[[cds@dim_reduce_type]]$root_cell <- my_root_cell
cds <- orderCells(cds)

save(cds,file=paste0(outpath,"/monocle2.cds"))

pdf("train_monocle2_Pseudotime.pdf",width=7,height=7)
plot_cell_trajectory(cds,color_by="Pseudotime") 
dev.off()

pdf("train_monocle2_stat1.pdf",width=7,height=7)
plot_cell_trajectory(cds,color_by="State")
dev.off()

pdf("train_monocle2_stat2.pdf",width=7,height=7)
plot_cell_trajectory(cds, color_by = "State") + facet_wrap("~State", ncol = 3)
dev.off()

pdf("train_monocle2_celltype1.pdf",width=7,height=7)
plot_cell_trajectory(cds,color_by="human_celltypes_250418",cell_size=1)+ geom_point(alpha = 0.1)
dev.off()

pdf("train_monocle2_celltype2.pdf",width=7,height=7)
plot_cell_trajectory(cds,color_by="human_celltypes_250418",cell_size=1,alpha=0.5)+ facet_wrap("~human_celltypes_250418", ncol = 3)
dev.off()

