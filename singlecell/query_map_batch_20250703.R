
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(Matrix)
library(readr)
library(reticulate)
library(dplyr)

ref <- readRDS("human_clustering_20250624.rds")
query1 <- readRDS("human_cs7_reanno.Rds")
query2 <- readRDS("human_cs8_reanno.Rds")

original_umap_model <- ref[["umap"]]
  
ref <- NormalizeData(ref) %>% FindVariableFeatures(nfeatures = 2000)  %>% ScaleData() %>% RunPCA()


query1 <- NormalizeData(query1) %>% FindVariableFeatures(nfeatures = 2000) %>% ScaleData() %>% RunPCA()
#p <- ElbowPlot(query1)
#ggsave("query1_elbowplot_pca.png", plot = p, width = 8, height = 6, dpi = 300)
query1 <- FindNeighbors(query1, dims = 1:30) %>% FindClusters()
#ref <- SetIdent(ref, value = "RNA_snn_res.1")
query1 <- RunUMAP(query1, dims = 1:30 ,seed.use = 42)


query2 <- NormalizeData(query2) %>% FindVariableFeatures(nfeatures = 2000) %>% ScaleData() %>% RunPCA()
p <- ElbowPlot(query2)
ggsave("query2_elbowplot_pca.png", plot = p, width = 8, height = 6, dpi = 300)
query2 <- FindNeighbors(query2, dims = 1:30) %>% FindClusters()
#ref <- SetIdent(ref, value = "RNA_snn_res.1")
query2 <- RunUMAP(query2, dims = 1:30 ,seed.use = 42)


# Step 6: 寻找锚点并迁移注释
anchors <- FindTransferAnchors(
  reference = ref,
  query = query1,
  dims = 1:30,
  reference.reduction = "pca"
)

predictions <- TransferData(anchorset = anchors, refdata = ref$reanno, dims = 1:30)
query1 <- AddMetaData(query1, metadata = predictions)

ref <- RunUMAP(ref, dims = 1:30, reduction = "pca", return.model = TRUE)
# 提取原始 UMAP embeddings（已处理过列名）
new_embeddings <- as.matrix(original_umap_model@cell.embeddings)
colnames(new_embeddings) <- c("umap_1", "umap_2")

# 替换 cell.embeddings（绘图用）
ref@reductions$umap@cell.embeddings <- new_embeddings

# 替换 misc$model$embedding（映射用）
ref@reductions$umap@misc$model$embedding <- new_embeddings

query1 <- MapQuery(anchorset = anchors, reference = ref, query = query1,
    refdata = list(celltype = "reanno"), reference.reduction = "pca", reduction.model = "umap")

p1 <- DimPlot(ref,reduction = "umap", group.by = "reanno", label = TRUE, label.size = 3,
    repel = TRUE) + NoLegend() + ggtitle("Reference annotations")
p2 <- DimPlot(query1,reduction = "ref.umap", group.by = "final_anno", label = TRUE,
    label.size = 3, repel = TRUE) + NoLegend() +  ggtitle("Query transferred labels")
p3 <- DimPlot(ref,reduction = "umap", group.by = "lineage", label = TRUE, label.size = 3,
    repel = TRUE) + NoLegend() + ggtitle("Reference annotations")
p4 <- DimPlot(ref,reduction = "umap", group.by = "orig.ident", label = TRUE, label.size = 3,
    repel = TRUE) + NoLegend() + ggtitle("Reference annotations")    
# 组合图形
combined_plot <- (p1 + p2 )/ (p3 + p4)

# 保存为 PNG 或 PDF 文件
ggsave("reference_and_query1_umap_integrate.png", plot = combined_plot, width = 16, height = 15, dpi = 300)




godsnot_102 <- c(
  "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
  "#FFDBE5", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF",
  "#997D87", "#5A0007", "#809693", "#6A3A4C", "#1B4400", "#4FC601", "#3B5DFF",
  "#4A3B53", "#FF2F80", "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92",
  "#FF90C9", "#B903AA", "#D16100", "#DDEFFF", "#000035", "#7B4F4B", "#A1C299",
  "#300018", "#0AA6D8", "#013349", "#00846F", "#372101", "#FFB500", "#C2FFED",
  "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09", "#00489C", "#6F0062",
  "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66", "#885578",
  "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
  "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F",
  "#938A81", "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757",
  "#C8A1A1", "#1E6E00", "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C",
  "#772600", "#D790FF", "#9B9700", "#549E79", "#FFF69F", "#201625", "#72418F",
  "#BC23FF", "#99ADC0", "#3A2465", "#922329", "#5B4534", "#FDE8DC", "#404E55",
  "#0089A3", "#CB7E98", "#A4E804", "#324E72"
)


# Define the ordered labels
final_anno_labels = c('TE', 'CTB_1','CTB_2', 'STB_1', 'STB_2', 'STB_3', 'EVT_1', 'EVT_2', 
                      'Epiblast_1', 'Epiblast_2', 'Epiblast_3','Ectodrm',
                      'Amniontic.epi','Amniontic.ectoderm',
                      'PGC',
                      'Primitive.streak',
                      'Neuromesodermal.progenitor',
                      'Neural.crest', 'Neural.ectoderm.forebrain', 'Neural.ectoderm.hindbrain', 'Neural.ectoderm.midbrain','Spinal.cord',
                      'Axial.mesoderm','Emergent.mesoderm','Pre-somatic.mesoderm','Somite', 'Myocyte.progenitor', 'Lateral.plate.mesoderm_1',
                      'Lateral.plate.mesoderm_2','Cardiac.mesoderm','Connecting.stalk','Amniotic.mesoderm','Exe.meso.progenitor','YS.mesoderm_1', 'YS.mesoderm_2',
                      'Hypoblast_1', 'Hypoblast_2', 'AVE', 'VE', 'YS.endoderm',
                      'DE','Gut',
                      'Notochord',
                      'Hemogenic.endothelial.progenitor','Endothelium','Erythroid','Primitive.megakaryocyte','Myeloid.progenitor'
)


query1$final_anno <- factor(query1$final_anno, levels = final_anno_labels, ordered = TRUE)

# Optionally, you can drop unused levels
query1$final_anno <- droplevels(query1$final_anno)

# Get the unique clusters from your dataset
unique_clusters <- unique(query1$final_anno)

# Ensure the number of colors matches the number of clusters
color_palette <- setNames(godsnot_102[1:length(final_anno_labels)], final_anno_labels)


# Generate the plot with the custom palette
p <- DimPlot(query1, reduction = "ref.umap", group.by = "final_anno") +
  scale_color_manual(values = color_palette) +  # Apply the custom color palette
  ggtitle("UMAP Plot with Custom Color Palette") 

pdf("human_CS7_RCTD_full_top_optimized.pdf", width = 13, height = 10)
p

dev.off()



# DimPlot 设置点的大小 + 自定义颜色 + 标题
p <- DimPlot(
  query1,
  reduction = "ref.umap",
  group.by = "final_anno",
  pt.size = 1.5,                # 调整点的大小（默认是 1）
  label = FALSE               # 是否显示标签（可选）
) +
  scale_color_manual(values = color_palette) +  # 应用自定义颜色
  ggtitle("UMAP Plot with Custom Color Palette") +
  theme(plot.title = element_text(hjust = 0.5))  # 标题居中

# 输出为 PDF
pdf("human_CS7_RCTD_full_top_optimized.pdf", width = 15, height = 8)
print(p)
dev.off()








# Step 6: 寻找锚点并迁移注释
anchors <- FindTransferAnchors(
  reference = ref,
  query = query2,
  dims = 1:30,
  reference.reduction = "pca"
)

predictions <- TransferData(anchorset = anchors, refdata = ref$reanno, dims = 1:30)
query2 <- AddMetaData(query2, metadata = predictions)

query2 <- MapQuery(anchorset = anchors, reference = ref, query = query2,
    refdata = list(celltype = "reanno"), reference.reduction = "pca", reduction.model = "umap")

p1 <- DimPlot(ref,reduction = "umap", group.by = "reanno", label = TRUE, label.size = 3,
    repel = TRUE) + NoLegend() + ggtitle("Reference annotations")
p2 <- DimPlot(query2,reduction = "ref.umap", group.by = "final_anno", label = TRUE,
    label.size = 3, repel = TRUE) + NoLegend() +  ggtitle("Query transferred labels")
p3 <- DimPlot(ref,reduction = "umap", group.by = "lineage", label = TRUE, label.size = 3,
    repel = TRUE) + NoLegend() + ggtitle("Reference annotations")
p4 <- DimPlot(ref,reduction = "umap", group.by = "orig.ident", label = TRUE, label.size = 3,
    repel = TRUE) + NoLegend() + ggtitle("Reference annotations")    
# 组合图形
combined_plot <- (p1 + p2 )/ (p3 + p4)

# 保存为 PNG 或 PDF 文件
ggsave("reference_and_query2_umap_integrate.png", plot = combined_plot, width = 16, height = 15, dpi = 300)


# Reorder the 'first_type' column according to the custom list and remove those that don't exist
query2$final_anno <- factor(query2$final_anno, levels = final_anno_labels, ordered = TRUE)

# Optionally, you can drop unused levels
query2$final_anno <- droplevels(query2$final_anno)

# Get the unique clusters from your dataset
unique_clusters <- unique(query2$final_anno)

# Ensure the number of colors matches the number of clusters
color_palette <- setNames(godsnot_102[1:length(final_anno_labels)], final_anno_labels)


library(ggplot2)
library(Seurat)
library(ggpubr)

# DimPlot 设置点的大小 + 自定义颜色 + 标题
p <- DimPlot(
  query2,
  reduction = "ref.umap",
  group.by = "final_anno",
  pt.size = 1.5,                # 调整点的大小（默认是 1）
  label = FALSE               # 是否显示标签（可选）
) +
  scale_color_manual(values = color_palette) +  # 应用自定义颜色
  ggtitle("UMAP Plot with Custom Color Palette") +
  theme(plot.title = element_text(hjust = 0.5))  # 标题居中

# 输出为 PDF
pdf("human_CS8_RCTD_full_top_optimized.pdf", width = 15, height = 8)
print(p)
dev.off()
