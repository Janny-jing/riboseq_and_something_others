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
write.table(deg,file="celltype_deg_monocle.xls",col.names=T,row.names=F,sep="\t",quote=F)
ordergene <- rownames(deg) 
cds <- setOrderingFilter(cds, ordergene)
head(fData(cds))
p1=plot_ordering_genes(cds)
print(p1)
ggsave(filename = paste0(outpath,"/plot_ordering_genes.pdf"),plot = p1,width = 6,height = 6, units ="in",device = "pdf")
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
plot_cell_trajectory(cds,color_by="human_celltypes_250418")
dev.off()

pdf("train_monocle2_celltype2.pdf",width=7,height=7)
plot_cell_trajectory(cds,color_by="human_celltypes_250418")+ facet_wrap("~human_celltypes_250418", ncol = 3)
dev.off()

