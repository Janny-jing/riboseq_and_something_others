library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(patchwork)
library(tidyverse)
library(writexl)
library(Matrix)
library(loupeR)

data1 <- readRDS("/storage/liuxiaodongLab/jiangjing/Projects/YutingFu/PD_YutingFu/20250604_brain_mouse/allsamples_6replicates_mnn_rm_mixed_rm_genes_add_1.5m_umap_annotationed_250514.rds")

rat <- data1
head(rat@meta.data)
mouse <- readRDS("/storage/liuxiaodongLab/jiangjing/Projects/YutingFu/PD_YutingFu/20250604_brain_mouse/rat_converted_sc_data.rds")
head(mouse@meta.data)

mouse$celltype <- mouse$cell_type1


mouse$dataset <- "Ximerakis"
rat$dataset <- "Our_data"
rat$sample <- rat$orig.ident
mouse$sample <- mouse$donor_age
table(rat$sample)
head(rat@meta.data)
rat$celltype <- rat$celltypes_250514_new

rownames(rat)
rownames(mouse)

rat <- subset(rat, idents = c("Graft-AC","Graft-Neuron","Graft-OLG"), invert = TRUE)
# 获取当前所有基因名
all_genes <- rownames(rat)

# 步骤1: 删除以"GRCh38"开头的基因
# 找出不以"GRCh38"开头的基因
keep_genes <- !grepl("^GRCh38", all_genes)

# 步骤2: 去除"mRatBN7.2-"前缀
# 先创建基因名副本
new_gene_names <- all_genes
# 只处理大鼠基因（保留不以"GRCh38"开头的基因）
new_gene_names[keep_genes] <- sub("^mRatBN7\\.2-", "", new_gene_names[keep_genes])

# 创建基因过滤后的对象
rat_filtered <- rat[keep_genes, ]

# 更新基因名
rownames(rat_filtered@assays$RNA@counts) <- new_gene_names[keep_genes]
rownames(rat_filtered@assays$RNA@data) <- new_gene_names[keep_genes]
rownames(rat_filtered) <- new_gene_names[keep_genes]  # 更新主rownames

sobj <- merge(rat_filtered, mouse)
sobj[["RNA"]] <- JoinLayers(sobj[["RNA"]])
counts <- LayerData(sobj, assay = "RNA", layer = "counts")

# non_zero_values <- counts@x
# hist(non_zero_values, breaks = 50, main = "Value Distribution of Counts", xlab = "Counts", ylab = "Frequency", col = "blue")

sobj <- CreateSeuratObject(counts = counts, meta.data = sobj@meta.data)

head(sobj@meta.data)
table(sobj$sample)
Idents(sobj) <- "sample"
levels(sobj)
levels(sobj)
new.cluster.ids=c("12M","18M","1.5M_5.5M_3.5M","1.5M_5.5M_3.5M","8M","1.5M_5.5M_3.5M","P1A","P1B","week_7_8","week_9_10_11","week_9_10_11","week_9_10_11","week_6","week_7_8")
names(new.cluster.ids) <- levels(sobj)
sobj <- RenameIdents(sobj, new.cluster.ids)
table(Idents(sobj))
sobj$sample <- Idents(sobj)


sobj[["RNA"]] <- split(sobj[["RNA"]], f = sobj$sample)

sobj <- NormalizeData(sobj)
sobj <- FindVariableFeatures(sobj)
sobj <- ScaleData(sobj)
sobj <- RunPCA(sobj)



# Integration method and identifier
integration_method <- "HarmonyIntegration"
identifier <- "harmony"

# integration_method <- "FastMNNIntegration"
# identifier <- "mnn"

# integration_method <- "RPCAIntegration"
# identifier <- "rpca"

# Construct the name for the new reduction
reduction_name <- paste('integrated', identifier, sep = '.')




# Perform integration dynamically
sobj <- IntegrateLayers(
  object = sobj, method = integration_method,
  new.reduction = reduction_name,
  verbose = TRUE)

# Re-join layers after integration
sobj[["RNA"]] <- JoinLayers(sobj[["RNA"]])

# Further analysis steps
sobj <- FindNeighbors(sobj, reduction = reduction_name, dims = 1:18)
sobj <- FindClusters(sobj, resolution = seq(0.1, 1, 0.1))
sobj <- RunUMAP(sobj, dims = 1:18, reduction = reduction_name)
head(sobj@meta.data)
# Save plots
DimPlot(sobj, group.by = "celltype", label=TRUE,repel=T,split.by = "dataset",pt.size=0.6)
DimPlot(sobj, group.by = "dataset",pt.size=0.7)
saveRDS(sobj, "../output/subcluster_human_cells_La_Manno_harmony_overlap_genes.rds")



sobj <- RunHarmony(
  sobj,
  group.by.vars = "dataset",
  reduction = "pca",
  reduction.save = "harmony",
  assay.use = "RNA",
  project.dim = FALSE
)

# Re-join layers after integration
sobj[["RNA"]] <- JoinLayers(sobj[["RNA"]])

# Further analysis steps
sobj <- FindNeighbors(sobj, reduction = "harmony", dims = 1:18)
sobj <- FindClusters(sobj, resolution = seq(0.1, 1, 0.1))
sobj <- RunUMAP(sobj, dims = 1:18, reduction = "harmony")
head(sobj@meta.data)
# Save plots
DimPlot(sobj, group.by = "celltype", label=TRUE,repel=T,split.by = "dataset",pt.size=0.6)
DimPlot(sobj, group.by = "dataset",pt.size=0.7)
saveRDS(sobj, "./last_merge_data.rds")

png("umap_plot_1.png", width = 4000, height = 1500, res = 300)
DimPlot(sobj,
        group.by = "celltype",
        label = TRUE,
        repel = TRUE,
        split.by = "dataset",
        pt.size = 0.6)
dev.off()


png("umap_plot_2.png", width = 2000, height = 1500, res = 300)
DimPlot(sobj,
        group.by = "dataset",
        pt.size=0.7)
dev.off()

pdf("umap_plot_1.pdf",width=20,height=8)
DimPlot(sobj,
        group.by = "celltype",
        label = FALSE,
        split.by = "dataset",
        pt.size = 0.6)
dev.off()

