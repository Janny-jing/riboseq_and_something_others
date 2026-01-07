
# 加载必要的包
library(Seurat)
library(tidyverse)
library(patchwork)
library(Matrix)
library(future)
library(stringr)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(future)
plan(multisession)  # 设置一个执行计划
# 提高限制（例如允许传输 40GB）
options(future.globals.maxSize = 40 * 1024^3) 

read_custom_10x <- function(data.dir) {
  # 检查文件是否存在
  required_files <- c("matrix.mtx.gz", "barcodes.tsv.gz", "features.tsv.gz")
  if(!all(file.exists(file.path(data.dir, required_files)))) {
    stop(paste("Missing files in", data.dir))
  }
  
  # 读取matrix
  matrix_path <- file.path(data.dir, "matrix.mtx.gz")
  mat <- Matrix::readMM(gzfile(matrix_path))
  
  # 读取barcodes
  barcode_path <- file.path(data.dir, "barcodes.tsv.gz")
  barcodes <- readLines(gzfile(barcode_path))
  
  # 读取features - 保留所有列
  feature_path <- file.path(data.dir, "features.tsv.gz")
  features_df <- read.delim(gzfile(feature_path), header = FALSE, stringsAsFactors = FALSE)
  
  # 检查列数
  if(ncol(features_df) == 1) {
    warning(paste("Features file in", data.dir, "has only 1 column. Using as gene names."))
    features_df$V2 <- features_df$V1  # 复制第一列为第二列
  } else if(ncol(features_df) < 2) {
    stop(paste("Invalid features file in", data.dir, "with", ncol(features_df), "columns"))
  }

  # 返回矩阵，确保行名为基因名，列名为细胞barcode
  colnames(mat) <- barcodes
  rownames(mat) <- features_df$V2  # 使用第二列为基因名

  return(mat)
}
# 1. 批量读取并处理h9样本 - 保持原始名称
ips_samples <- c("iPS3M_A", "iPS3M_B", "iPS3M_C", "iPS3M_D","iPS3M_E","iPS3M_F","iPS8M_A", "iPS8M_B", "iPS8M_C", "iPS8M_D","iPS8M_E","iPS8M_F")
ips_list <- list()
qc_list <- list()

# 创建保存图片的主目录
output_dir <- "qc_violin_plots"
dir.create(output_dir, showWarnings = FALSE)

for (sample_name in ips_samples) {
  cat("\nProcessing sample:", sample_name, "\n")
  
  # 使用自定义函数读取
  data <- read_custom_10x(data.dir = sample_name)
  
  # 检查数据是否为空
  if(ncol(data) == 0 || nrow(data) == 0) {
    warning(paste("Empty data for sample:", sample_name))
    next
  }
  
  # 创建Seurat对象 - 使用原始名称
  seu <- CreateSeuratObject(
    counts = data,
    project = sample_name,
    min.cells = 3,
    min.features = 200
  )
  
  # 记录原始细胞数
  orig_cells <- ncol(seu)
  cat("Original cells:", orig_cells, "\n")
  
  # 添加元数据
  seu$batch <- sample_name
  seu$dataset <- "New_ips"
  
  # 自动检测MT基因模式
  mt_patterns <- c("^MT-", "^Mt-", "^mt-", "-MT$", "-mt$", "^MT[0-9]", "^mt[0-9]")
  mt_found <- sapply(mt_patterns, function(p) any(grepl(p, rownames(seu))))
  
  if(any(mt_found)) {
    mt_pattern <- mt_patterns[which(mt_found)[1]]
    seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = mt_pattern)
    cat("Using mitochondrial pattern:", mt_pattern, "\n")
  } else {
    seu[["percent.mt"]] <- 0
    warning("No mitochondrial genes detected in", sample_name)
  }

  # +++ 质控前 小提琴图 +++
  p_before <- VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0)
  ggsave(filename = file.path(output_dir, paste0(sample_name, "_before_qc.png")),
         plot = p_before, width = 12, height = 4, dpi = 300)

  # 质控过滤 - 严格标准
  keep_cells <- seu$nFeature_RNA > 300 & 
                seu$nFeature_RNA < 8000 &
                seu$nCount_RNA > 2000 &
                seu$nCount_RNA < 25000 & 
                seu$percent.mt < 5
  
  # 处理可能的NA值
  keep_cells[is.na(keep_cells)] <- FALSE
  
  seu <- seu[, keep_cells]
  
  # 记录过滤后细胞数
  filtered_cells <- ncol(seu)
  cat("Filtered cells:", filtered_cells, 
      "(", round(filtered_cells/orig_cells*100, 1), "% retained)\n")
  
  # +++ 质控后 小提琴图 +++
  p_after <- VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0)
  ggsave(filename = file.path(output_dir, paste0(sample_name, "_after_qc.png")),
         plot = p_after, width = 12, height = 4, dpi = 300)

  # 创建QC指标数据框用于记录
  qc_stats <- data.frame(
    Sample = sample_name,
    Original_Cells = orig_cells,
    Filtered_Cells = filtered_cells,
    Genes_Median = median(seu$nFeature_RNA),
    UMI_Median = median(seu$nCount_RNA),
    MT_Median = median(seu$percent.mt)
  )
  
  # 添加到结果列表
  ips_list[[sample_name]] <- seu
  qc_list[[sample_name]] <- qc_stats
}

sample_names <- names(ips_list) 
ips_list <- merge(x = ips_list[[1]], 
                       y = ips_list[-1], 
                       add.cell.ids = sample_names)

ips_list <- NormalizeData(ips_list) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA()

ips_list <- IntegrateLayers(
  object = ips_list, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca",
  verbose = FALSE
)

ips_list[["RNA"]] <- JoinLayers(ips_list[["RNA"]])

ips_list <- FindNeighbors(ips_list, reduction = "integrated.rpca")
ips_list <- FindClusters(ips_list, resolution = 0.01)
ips_list <- RunUMAP(ips_list, dims = 1:30, reduction = "integrated.rpca", seed.use = 42)
saveRDS(ips_list, "ips_3m_8m_integrated_rpca_raw.rds")


ips_list <- ips
all_clusters <- as.character(unique(Idents(ips_list)))

# 创建一个空列表，保存每个 cluster 处理后的子集
filtered_clusters <- list()

# 遍历每个 cluster
for (clust in all_clusters) {
  # 提取当前 cluster
  sub <- subset(ips_list, subset = RNA_snn_res.0.01 == clust)
  
  # 根据 cluster 编号筛选基因
  if (clust %in% c("0", "1", "2")) {
    # cluster 0~2：只保留鼠源基因（以 mRatBN7.2- 开头）
    species_genes <- grep("^mRatBN7\\.2\\-", rownames(sub), value = TRUE)
  } else if (clust == "3") {
    # cluster 3：只保留人源基因（以 GRCh38.109- 开头）
    species_genes <- grep("^GRCh38\\.109\\-", rownames(sub), value = TRUE)
  } else {
    # 其他 cluster（如果有的话）保留所有基因
    species_genes <- rownames(sub)
  }
  
  # 筛选表达矩阵
  sub_species <- sub[species_genes, ]
  
  # 保存到列表中
  filtered_clusters[[clust]] <- sub_species
}

# 合并所有 cluster
ips_filtered <- merge(x = filtered_clusters[[1]], y = filtered_clusters[-1])
ips_filtered[["RNA"]] <- JoinLayers(ips_filtered[["RNA"]])


existing_list <- SplitObject(ips_filtered, split.by = "orig.ident")

sample_names <- names(existing_list) 
existing_list <- merge(x = existing_list[[1]], 
                       y = existing_list[-1], 
                       add.cell.ids = sample_names)

existing_list <- NormalizeData(existing_list) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA()

existing_list <- IntegrateLayers(
  object = existing_list, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca",
  verbose = FALSE
)

existing_list[["RNA"]] <- JoinLayers(existing_list[["RNA"]])

existing_list <- FindNeighbors(existing_list, reduction = "integrated.rpca")
existing_list <- FindClusters(existing_list)
existing_list <- RunUMAP(existing_list, dims = 1:30, reduction = "integrated.rpca", seed.use = 42)
saveRDS(existing_list, "umap_ips_3m_8m_integrated_rpca.rds")
############################################################################################################################
# 检查并清理 RNA assay 中无用的 scale.data.* layers
assay <- existing_list[["RNA"]]
# 获取所有 layer 名称
layer_names <- names(assay@layers)
# 使用正则表达式匹配以 "scale.data." 开头且后跟数字的 layer，但排除 "scale.data" 本身
to_remove <- grep("^scale\\.data\\.[0-9]+$", layer_names, value = TRUE)
if (length(to_remove) > 0) {
  # 批量删除
  assay@layers[to_remove] <- NULL
  # 更新回 Seurat 对象
  existing_list[["RNA"]] <- assay
  message("✅ 已删除以下无用 layers: ", paste(to_remove, collapse = ", "))
} else {
  message("🔍 未发现多余的 scale.data.* layers。")
}
# 验证结果
cat("\n📋 当前 RNA assay 中的 layers:\n")
print(names(existing_list[["RNA"]]@layers))

############################################################################################################################
# 3. 处理现有数据集
ref <- readRDS("allsamples_6replicates_mnn_rm_mixed_rm_genes_add_1.5m_umap_annotationed_250514.rds")
existing_list <- readRDS("umap_ips_3m_8m_integrated_rpca.rds")
anchors <- FindTransferAnchors(
  reference = ref,
  query = existing_list,
  dims = 1:30,
  reference.reduction = "pca"
)

predictions <- TransferData(anchorset = anchors, refdata = ref$celltypes_250514_new, dims = 1:30)
existing_list <- AddMetaData(existing_list, metadata = predictions)

ref <- RunUMAP(ref, dims = 1:30, reduction = "pca", return.model = TRUE)

existing_list <- MapQuery(anchorset = anchors, reference = ref, query = existing_list,
    refdata = list(celltype = "celltypes_250514_new"), reference.reduction = "mnn", reduction.model = "umap")

p1 <- DimPlot(ref,reduction = "umap", group.by = "celltypes_250514_new", label = TRUE, label.size = 3,
    repel = TRUE)  + ggtitle("Reference annotations")
p2 <- DimPlot(existing_list,reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE,
    label.size = 3, repel = TRUE)  +  ggtitle("Query transferred labels")
p3 <- DimPlot(ref,reduction = "umap", group.by = "orig.ident", label = TRUE, label.size = 3,
    repel = TRUE) + ggtitle("Reference annotations")
p4 <- DimPlot(existing_list,reduction = "ref.umap", group.by = "orig.ident", label = TRUE, label.size = 3,
    repel = TRUE) + ggtitle("Reference annotations")    
# 组合图形
combined_plot <- (p1 + p2 )/ (p3 + p4)

# 保存为 PNG 或 PDF 文件
ggsave("reference_and_query1_umap_integrate.png", plot = combined_plot, width = 20, height = 15, dpi = 300)


existing_list <- FindNeighbors(existing_list, dims = 1:30) %>% FindClusters(resolution = 1.1)
#existing_list <- SetIdent(existing_list, value = "RNA_snn_res.1")
existing_list <- RunUMAP(existing_list, dims = 1:30 seed.use = 42)



p <- DimPlot(exiting_list,reduction = "umap", group.by = "RNA_snn_res.1.1", label = TRUE, label.size = 3,
    repel = TRUE) 

# 保存为 PNG 或 PDF 文件
ggsave("umap_1.1.png", plot = p, width = 20, height = 15, dpi = 300)
p <- DimPlot(exiting_list,reduction = "umap", group.by = "RNA_snn_res.1.3", label = TRUE, label.size = 3,
    repel = TRUE) 

# 保存为 PNG 或 PDF 文件
ggsave("umap_1.3.png", plot = p, width = 20, height = 15, dpi = 300)

p <- DimPlot(exiting_list,reduction = "umap", group.by = "RNA_snn_res.1.4", label = TRUE, label.size = 3,
    repel = TRUE) 

# 保存为 PNG 或 PDF 文件
ggsave("umap_1.4.png", plot = p, width = 20, height = 15, dpi = 300)

p <- DimPlot(exiting_list,reduction = "umap", group.by = "predicted.celltype", label = TRUE, label.size = 3,
    repel = TRUE) 

# 保存为 PNG 或 PDF 文件
ggsave("umap_predicted_celltype.png", plot = p, width = 20, height = 15, dpi = 300)









marker_genes <- list(
  # 移植来源细胞类型
  "Graft-Neuron" = c("GAD2", "INA", "NEFL", "TMEM130"),
  "Graft-AC" = c("GJA1", "AQP4", "GFAP", "SLC1A3"),
  "Graft-OLG" = c("APOD", "OLIG1", "OLIG2", "SOX10", "MBP"),
  # 纹状体中型多棘神经元
  "MSN-Drd1" = c("Drd1", "Tac1", "Rgs9", "Ppp1r1b"),
  "MSN-Drd2" = c("Drd2", "Penk", "Gpr88", "Rgs9"),
  # 兴奋性谷氨酸能神经元 (TEGLU)
  "TEGLU-Slc30a3" = c("Slc30a3", "Pdzrn3", "Cdh12", "Col6a1", "Slc17a7", "Mef2c", "Nptxr"),
  "TEGLU-Rxfp1" = c("Rxfp1", "Hs3st2", "Col6a1", "Slc30a3", "Slc17a7", "Mef2c", "Cemip2"),
  "TEGLU-Nxph3" = c("Nxph3", "Zfpm2", "Slc17a7", "Hs3st4", "Foxp2", "Sv2b"),
  "TEGLU-Abi3bp" = c("Abi3bp", "Slc17a7", "Cpne4", "Mef2c"),
  "TEGLU-ETV1" = c("Etv1", "Bmpr1b", "Sv2b", "Slc17a7", "Cadps2"),
  "TEGLU-Nr4a2" = c("Nr4a2", "Tfap2d", "Hs3st4", "Slc17a7", "Sv2b"),
  "TEGLU-Grp" = c("Grp", "Hs3st4", "Sv2b", "Slc17a7", "Plcxd3"),
  # 抑制性 GABA 能神经元 (TEINI)
  "TEINI-Gpc3" = c("Gpc3", "Gad1", "Gad2", "St8sia4"),
  "TEINI-Meis2" = c("Meis2", "Gad1", "Gad2", "Tac1", "Cntnap5c"),
  "TEINI-Sst" = c("Sst", "Reln", "Cntnap4", "Lypd6b", "Gad1", "Gad2"),
  "TEINI-Hpse" = c("Hpse", "Pld5", "Rerg", "Tac1", "Gad1", "Gad2"),
  "TEINI-Vip" = c("Vip", "Dlx1", "Npy", "Gad1", "Cck", "Slc32a1", "Gad2", "Cnr1"),
  "TEINI-Il1rapl2" = c("Il1rapl2", "Lypd6b", "Gad1", "Gad2"),
  # 胶质细胞
  "AC" = c("Gfap", "Aqp4", "Phkg1", "Slco1c1"),
  "OLG" = c("Mog", "Mobp", "Mbp", "Mal", "Gjc3"),
  "OPC" = c("Pdgfra", "Stk32a", "Vcan", "Olig2", "Sox10", "Olig1", "Ptprz1"),
  # 非神经元脑间质细胞
  "EPEN" = c("Ccdc153", "Aqp4", "Gfap", "Sox2"),
  "VLMC" = c("Col1a1", "Igf2", "Ptgds", "Pdgfrb", "Prdm6", "Eya2", "Foxd1"),
  # 血管相关细胞
  "Vascular cells" = c("Mecom", "Flt1", "Abcc9", "Vtn", "Prom1", "Myl9")
)

names(marker_genes)
############################################################################################################33
####这里缺少了一个重复基因处理
# 获取当前基因名（直接从 Seurat 对象提取）
old_genes <- rownames(existing_list)
# 去除前缀（以第一个 "-" 之后的内容为准）
new_genes <- str_replace(old_genes, "^[^-]+-", "")
# 查看前后对比
head(data.frame(Old = old_genes, New = new_genes))
# 获取默认 assay 名称
assay_name <- DefaultAssay(existing_list)
# 更新基因名（适用于 Assay5 / SpectraAssay）
rownames(existing_list) <- new_genes
# 更新 Seurat 对象（确保内部一致性）
existing_list <- UpdateSeuratObject(existing_list)
# 确认是否成功
head(rownames(existing_list))


# 设置输出目录
output_dir <- "marker_umap_plots"
dir.create(output_dir, showWarnings = FALSE)

# 遍历每个细胞类型
for (celltype in names(marker_genes)) {
  
  # 获取当前细胞类型的 marker 基因
  markers <- marker_genes[[celltype]]
  
  # 检查这些基因是否存在于 Seurat 对象中
  valid_markers <- intersect(markers, rownames(existing_list))
  
  if (length(valid_markers) == 0) {
    cat("No valid markers for", celltype, "\n")
    next
  }
  
  # 创建每个 marker 的 UMAP 图
  plots <- list()
  for (gene in valid_markers) {
    p <- FeaturePlot(existing_list, features = gene, reduction = "umap", raster = FALSE) +
      ggtitle(gene)
    plots[[gene]] <- p
  }
  
  # 将图组合成每行最多4个图
  combined_plot <- wrap_plots(plots, ncol = 4)
  
  # 静默保存图像（不显示）
  png_file <- file.path(output_dir, paste0(celltype, "_markers_umap.png"))
  grDevices::png(
    filename = png_file,
    width = 1600,
    height = 400 * ceiling(length(plots) / 4),
    res = 150
  )
  print(combined_plot)
  dev.off()
  
  cat("Saved plot for", celltype, "\n")
}


########################################################################################


# 1. 定义注释列表
cluster_annotation <- list(
  # 移植来源细胞类型
  "Graft-Neuron" = c(7,22,24,36,37),
  "Graft-AC" = c(14),
  "Graft-OLG" = c(28),
  # 纹状体中型多棘神经元
  "MSN-Drd1" = c(0,6,16),
  "MSN-Drd2" = c(1,3,26),
  # 兴奋性谷氨酸能神经元 (TEGLU)
  "TEGLU-Slc30a3" = c(2,9,25,38,40,46),
  "TEGLU-Rxfp1" = c(4),
  "TEGLU-Nxph3" = c(12),
  "TEGLU-Abi3bp" = c(21),
  "TEGLU-ETV1" = c(33),
  "TEGLU-Nr4a2" = c(31),
  "TEGLU-Grp" = c(32),
  # 抑制性 GABA 能神经元 (TEINI)
  "TEINI-Gpc3" = c(13,20),
  "TEINI-Meis2" = c(10),
  "TEINI-Sst" = c(15),
  "TEINI-Hpse" = c(11),
  "TEINI-Vip" = c(18,41),
  "TEINI-Il1rapl2" = c(29),
  "TEINI-ChaT" = c(34),
  # 胶质细胞
  "AC" = c(5,27,48),
  "OLG" = c(8,39),
  "OPC" = c(17,43),
  # 非神经元脑间质细胞
  "EPEN" = c(35),
  "VLMC" = c(30),
  # 血管相关细胞
  "Vascular cells" = c(19,44,45,47),
  "Microglia" = c(23,42)
)

# 1. 设置当前 ident 为 RNA_snn_res.1.1
Idents(existing_list) <- "RNA_snn_res.1.1"

# 2. 提取 cluster ID（注意：它是字符型，需转为 numeric）
clusters <- as.numeric(as.character(Idents(existing_list)))

# 3. 创建 cell_type 向量，默认为 "Unknown"
cell_types <- rep("Unknown", length(clusters))

# 4. 根据你的 annotation list 进行映射
for (cell_type_name in names(cluster_annotation)) {
  target_clusters <- cluster_annotation[[cell_type_name]]
  # 匹配 cluster IDs 并赋值
  cell_types[clusters %in% target_clusters] <- cell_type_name
}

# 5. 将注释写入 meta.data
existing_list$cell_type <- factor(cell_types,
                                  levels = c(names(cluster_annotation), "Unknown"))

# 6. （可选）将注释后的 cell_type 设为当前 ident
Idents(existing_list) <- "cell_type"

# 7. 查看注释分布
table(Idents(existing_list))


p <- DimPlot(existing_list,reduction = "umap", group.by = c("RNA_snn_res.1.1","cell_type"), label = TRUE, label.size = 3,
    repel = TRUE) 

# 保存为 PNG 或 PDF 文件
ggsave("umap_1.1_celltype.png", plot = p, width = 30, height = 15, dpi = 300)

saveRDS(existing_list, "umap_ips_3m_8m_integrated_rpca_annocelltype_250722.rds")