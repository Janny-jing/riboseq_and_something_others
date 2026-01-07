library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(patchwork)
library(tidyverse)
library(writexl)
library(Matrix)
library(loupeR)

all <- readRDS("/storage/liuxiaodongLab/yanxixi/Projects/YutingFu/PD_YutingFu/2025-03-13_102_Heatmap-of-AUROC-in-human-cells-with-public-data/data/our_data_human_overlap_genes.rds")

human <- all
head(human@meta.data)
chen <- readRDS("/storage/liuxiaodongLab/yanxixi/Projects/YutingFu/PD_YutingFu/2025-03-13_102_Heatmap-of-AUROC-in-human-cells-with-public-data/data/human_ventral_midbrain_2016_cells_overlap_genes.rds")
head(chen@meta.data)

chen$celltype <- chen$Cell_Type


chen$dataset <- "La_manno"
human$dataset <- "Our_data"
human$sample <- human$stages
chen$sample <- chen$Timepoint
table(human$sample)
head(human@meta.data)
human$celltype <- human$human_celltypes_250330

rownames(human)
rownames(chen)

Idents(human) <- "human_celltypes_250330"
human <- subset(human, idents=c("DAN-1","DAP-3-wanted","DAP-3-unwanted","DAP-1","DAP-2","DAN-2","DAP-4"))
Idents(chen) <- "Cell_Type"
chen <- subset(chen, idents=c("hDA0","hDA1","hDA2","hNProg","hProgBP","hProgFPL","hProgFPM","hProgM","hRgl1"))

sobj <- merge(human, chen)
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
# integration_method <- "HarmonyIntegration"
# identifier <- "harmony"

# integration_method <- "FastMNNIntegration"
# identifier <- "mnn"

integration_method <- "RPCAIntegration"
identifier <- "rpca"

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
saveRDS(sobj, "../output/subcluster_human_cells_La_Manno_rpca_overlap_genes.rds")
