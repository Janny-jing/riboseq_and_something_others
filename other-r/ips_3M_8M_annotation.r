
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
ips_samples <- c("ips1.5M_A","ips1.5M_B","ips1.5M_B","iPS3M_A", "iPS3M_B", "iPS3M_C", "iPS3M_D","iPS3M_E","iPS3M_F","iPS8M_A", "iPS8M_B", "iPS8M_C", "iPS8M_D","iPS8M_E","iPS8M_F")
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
saveRDS(ips_list, "ips_1.5m_3m_8m_integrated_rpca_raw.rds")

########################################################################################3
##筛选人和鼠的cluster，提取人的cluster，重新降维聚类，并与已注释数据进行印射
# 提取 meta.data
ips_list <- readRDS("ips_1.5m_3m_8m_integrated_rpca_raw.rds")
all_genes <- rownames(ips_list)
is_human <- grepl("^GRCh38\\.109-", all_genes)
is_rat   <- grepl("^mRatBN7\\.2-", all_genes)
# 获取人类和大鼠基因列表
human_genes <- all_genes[is_human]
rat_genes   <- all_genes[is_rat]
ips_list[["percent_human"]] <- PercentageFeatureSet(ips_list, features = human_genes)
ips_list[["percent_rat"]]   <- PercentageFeatureSet(ips_list, features = rat_genes)

metadata <- ips_list@meta.data
metadata$is_human <- ifelse(metadata$percent_human > 50, "Human", "Rat")
# 查看每 cluster 中 Human/Rat 细胞数量分布
table(metadata$seurat_clusters, metadata$is_human)
# 绘制小提琴图
png("violin_1.png", width = 4000, height = 1500, res = 300)
ggplot(metadata, aes(x = factor(RNA_snn_res.0.01), y = percent_human)) +
  geom_violin(trim = FALSE) +
  labs(title = "Human Gene Expression Proportion by Cluster",
       x = "Cluster",
       y = "Percent Human") +
  theme_classic()
dev.off()

p <- VlnPlot(ips_list, features = "percent_human", pt.size = 0)
ggsave(filename =  "violin_2.png",
         plot = p, width = 12, height = 4, dpi = 300)


# 绘制umap图
# 定义一个从蓝色（低）到红色（高）的颜色渐变
my_cols <- brewer.pal(n = 11, name = "RdBu")  # 11阶红蓝渐变
my_cols <- colorRampPalette(my_cols)(100)     # 插值到100种颜色
# 绘图并保存
png("feature_1.png", width = 2000, height = 1500, res = 300)
FeaturePlot(ips_list, features = "percent_human")  # 使用自定义颜色
dev.off()

p <- DimPlot(ips_list,group.by = "RNA_snn_res.0.01", label = TRUE, label.size = 3,
    repel = TRUE)  
ggsave("umap_filter_human.png", plot = p, width = 12, height = 10, dpi = 300)

#提取cluster3，并剔除其中的鼠源基因
table(Idents(ips_list))
cluster2 <- subset(ips_list, subset = RNA_snn_res.0.01 == 2)
dim(cluster2)
# 找出人源基因
human_genes <- grep("^GRCh38\\.109-", rownames(cluster2), value = TRUE)

# 筛选表达矩阵
cluster2_human <- cluster2[human_genes, ]

# 查看维度变化
dim(cluster2)         # 原始数据
dim(cluster2_human)   # 只保留 human 基因后


existing_list <- SplitObject(cluster2_human, split.by = "orig.ident")
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
p <- ElbowPlot(existing_list)
ggsave("elbowplot_pca.png", plot = p, width = 8, height = 6, dpi = 300)
existing_list <- FindNeighbors(existing_list, reduction = "integrated.rpca", dims = 1:30)
existing_list <- FindClusters(existing_list)
existing_list <- RunUMAP(existing_list, dims = 1:30, reduction = "integrated.rpca", seed.use = 42)
# 查看原始基因名
head(rownames(existing_list))
# 去除前缀 GRCh38.109-
new_gene_names <- sub("^GRCh38\\.109-", "", rownames(existing_list))
# 设置新的基因名
rownames(existing_list) <- new_gene_names
# 检查结果
head(rownames(existing_list))
saveRDS(existing_list, "umap_ips_3m_8m_integrated_rpca.rds")
############################################################################################################################
##与h9注释rds印射
# 3. 处理现有数据集
ips <- readRDS("umap_ips_3m_8m_integrated_rpca.rds")
h9 <- readRDS("new_umap_dim18_spread2.5_mindist0.7_celltypes_250418_orig.cellnames.rds")
##去掉移植前的细胞
h9_filtered <- subset(h9, subset = stages != "P1A" & stages != "P1B")
original_umap_model <- h9_filtered[["umap"]]
anchors <- FindTransferAnchors(
  reference = h9_filtered,
  query = ips,
  dims = 1:30,
  reference.reduction = "pca"
)

predictions <- TransferData(anchorset = anchors, refdata = h9_filtered$human_celltypes_250418, dims = 1:30)
ips <- AddMetaData(ips, metadata = predictions)

h9_filtered <- RunUMAP(h9_filtered, dims = 1:30, reduction = "pca", return.model = TRUE)
# 提取原始 UMAP embeddings（已处理过列名）
new_embeddings <- as.matrix(original_umap_model@cell.embeddings)
colnames(new_embeddings) <- c("umap_1", "umap_2")

# 替换 cell.embeddings（绘图用）
h9_filtered@reductions$umap@cell.embeddings <- new_embeddings

# 替换 misc$model$embedding（映射用）
h9_filtered@reductions$umap@misc$model$embedding <- new_embeddings
ips <- MapQuery(anchorset = anchors, reference = h9_filtered, query = ips,
    refdata = list(celltype = "human_celltypes_250418"), reference.reduction = "pca", reduction.model = "umap")

p1 <- DimPlot(h9_filtered,reduction = "umap", group.by = "human_celltypes_250418", label = TRUE, label.size = 3,
    repel = TRUE)  + ggtitle("Reference annotations")
p2 <- DimPlot(ips,reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE,
    label.size = 3, repel = TRUE)  +  ggtitle("Query transferred labels")
p3 <- DimPlot(h9_filtered,reduction = "umap", group.by = "orig.ident", label = TRUE, label.size = 3,
    repel = TRUE) + ggtitle("Reference annotations")
p4 <- DimPlot(ips,reduction = "ref.umap", group.by = "orig.ident", label = TRUE, label.size = 3,
    repel = TRUE) + ggtitle("Reference annotations")    
# 组合图形
combined_plot <- (p1 + p2 )/ (p3 + p4)

# 保存为 PNG 或 PDF 文件
ggsave("reference_and_query1_umap_integrate.png", plot = combined_plot, width = 16, height = 15, dpi = 300)
############################################################################################################################
# 检查并清理 RNA assay 中无用的 scale.data.* layers
assay <- ips[["RNA"]]
# 获取所有 layer 名称
layer_names <- names(assay@layers)
# 使用正则表达式匹配以 "scale.data." 开头且后跟数字的 layer，但排除 "scale.data" 本身
to_remove <- grep("^scale\\.data\\.[0-9]+$", layer_names, value = TRUE)
if (length(to_remove) > 0) {
  # 批量删除
  assay@layers[to_remove] <- NULL
  # 更新回 Seurat 对象
  ips[["RNA"]] <- assay
  message("✅ 已删除以下无用 layers: ", paste(to_remove, collapse = ", "))
} else {
  message("🔍 未发现多余的 scale.data.* layers。")
}
# 验证结果
cat("\n📋 当前 RNA assay 中的 layers:\n")
print(names(ips[["RNA"]]@layers))
############################################################################################################################
# 确保原始对象不被修改
sobj <- ips

# 创建参数组合网格
param_grid <- expand.grid(
  dims = seq(20, 31, 5),         
  k = seq(20, 31, 5),            
  #resolution = seq(0.5, 2.1, 1), 
  min.dist = seq(0.5, 1, 0.2),  
  spread = seq(0.5, 3, 0.5)      
)

# 循环测试所有参数组合
for(i in 1:nrow(param_grid)) {

  # 获取当前参数
  curr <- param_grid[i, ]
  d <- curr$dims
  k_val <- curr$k
  #res <- curr$resolution
  md <- curr$min.dist
  sp <- curr$spread
  
  tryCatch({
    # 使用当前参数处理数据
    sobj <- FindNeighbors(sobj, 
                          reduction = "integrated.rpca", 
                          dims = 1:d, 
                          k = k_val)
    
    sobj <- FindClusters(sobj, resolution = 0.5)
    
    sobj <- RunUMAP(sobj, 
                    dims = 1:d,
                    reduction = "integrated.rpca",
                    min.dist = md,
                    spread = sp,
                    n.neighbors = k_val,  # 与FindNeighbors的k保持一致
                    seed.use = 42)
    
    # 生成唯一文件名（包含所有参数）
    fname <- sprintf(
      "umap_d%02d_k%02d_md%.1f_sp%.1f.png",
      d, k_val, md, sp
    )
    param_title <- paste(
       sprintf("d:%d k:%d  min.dist:%.1f spread:%.1f",
                d, k_val, md, sp)
    )
    
    # 保存高质量UMAP图
    png(fname, width = 4000, height = 1800, res = 300)  # 增加高度容纳标题
    
    print(
      DimPlot(sobj, 
              group.by = c("seurat_clusters","predicted.celltype"), 
              label = TRUE,
              repel = TRUE,
              pt.size = 0.6) +
        ggtitle(param_title) +  # 添加参数标题
        theme(plot.title = element_text(size = 22, face = "bold", 
                                        hjust = 0.5, margin = margin(b = 15)))
    )
    dev.off()
    
    # 可选：保存中间结果（谨慎使用，文件较大）
    # saveRDS(sobj, paste0("result_", fname, ".rds"))
    
    message(sprintf("已完成: %s", fname))
    
  }, error = function(e) {
    message(sprintf("参数组合失败 (d=%d,k=%d,md=%.1f,sp=%.1f): %s",
                    d, k_val, md, sp, e$message))
  })
}



ips <- FindNeighbors(ips, dims = 1:20, reduction = "integrated.rpca") %>% FindClusters(resolution = 0.5)
#ips <- SetIdent(ips, value = "RNA_snn_res.0.5")
ips <- RunUMAP(ips, dims = 1:20 ,seed.use = 42, reduction = "integrated.rpca")
p <- DimPlot(ips,reduction = "umap", group.by = "cell_type", label = TRUE, label.size = 3,repel = TRUE)
# 保存为 PNG 或 PDF 文件
ggsave("umap_all.png", plot = p, width = 20, height = 10, dpi = 300)

##################################################################################################################
markers <- list(
  "Neural stem cells" = c("TNFRSF12A","MKI67","TOP2A","PAX3","SOX1","IGFBPL1"),
  "Rgl1-like" = c("top2a", "mki67", "Lmx1a", "Apcdd1", "OTX2", "TTR"),
  "Glial progenitors" = c("nos2", "top2a", "mki67", "col4a6", "col4a5", "NES"),
  "Neural progenitors" = c("TLX3", "HES6", "NEUROD4", "IGFBPL1", "PAX3", "CXCR4"),
  "Astrocyte progenitor cells" = c("AQP4", "GFAP", "SPARC", "SPARCL1"),
  "OPC" = c("PDGFRA", "OLIG1", "OLIG2", "SOX10"),
  "OLG" = c("MAG", "MOG", "MYRF", "SOX10"),
  "DA progenitor-FABP7" = c("FABP7", "TTR", "Lmx1a", "Apcdd1", "OTX2"),
  "DA progenitor-CHRNB3 high" = c("CHRNA6", "CHRNB3", "LMX1A"),
  "DA progenitor-CHRNB3 low" = c("DDC", "ASCL1", "LMX1A"),
  "DA progenitor-terminal" = c("LMX1A", "EN1", "NR4A2","DDC"),
  "Immature DA" = c("KCNJ6", "TH", "SLC6A3", "SLC18A2", "PITX3", "NR4A2"),
  "Mature DA" = c("TH", "ALDH1A1", "NR4A2", "SLC18A2", "KCNJ6", "SLC6A3"),
  "VLMC" = c("FN1", "S100A11", "CDH6", "FBLN1"),
  "Astrocyte" = c("ALDH1L1", "SLC1A3", "AQP4", "GFAP"),
  "GLUL-1" = c("LHX9", "GRM4", "GRIN1", "SLC17A6","LHX2", "GLS2"),
  "GLUL-2" = c("GRM4", "GRIN1", "GRIN2B", "SLC17A6", "GLS2"),
  "Unknown-2" = c("TAC1", "C1QL2", "C1QL4","NKX6-1", "NKX6-2","SLC17A6"),
  "Unknown-1" = c("PROX1", "ONECUT1", "ONECUT2", "ONECUT3"),
  "GABAL-1" = c("SLC32A1", "GAD2", "GAD1"),
  "GABAL-2" = c("SLC32A1", "GAD2", "GAD1","SST","DLX1"),
  "NE neuron" = c("DBH", "PHOX2B", "ASCL1"),
  "Motor neuron" = c("SLC18A3", "SLC5A7", "CHAT", "ISL1"),
  "5-HT" = c("SLC6A4", "TPH2", "GATA3")
)



##########################################################################################################
##选择合适的umap进行注释
markers <- list(
  "Rgl1-like" = c("top2a", "mki67", "Lmx1a", "Apcdd1", "OTX2", "TTR"),
  "DA progenitor-FABP7" = c("FABP7", "TTR", "Lmx1a", "Apcdd1", "OTX2"),
  "DA progenitor-CHRNB3 high" = c("CHRNA6", "CHRNB3", "LMX1A","SCN1A-AS1", "SCN9A", "SCN5A", "SCNN1B", "SLC24A3"),
  "DA progenitor-CHRNB3 low" = c("DDC", "ASCL1", "LMX1A"),
  "DA progenitor-terminal" = c("LMX1A", "EN1", "NR4A2", "DDC"),
  "Immature DA" = c("KCNJ6", "TH", "SLC6A3", "SLC18A2", "PITX3", "NR4A2"),
  "Mature DA" = c("TH", "ALDH1A1", "NR4A2", "SLC18A2", "KCNJ6", "SLC6A3")
)
get_valid_genes <- function(sobj, gene_list) {
  assay_name <- DefaultAssay(sobj)
  all_genes <- rownames(GetAssayData(sobj, assay = assay_name, slot = "data"))
  
  valid_genes <- intersect(tolower(gene_list), tolower(all_genes))
  matched_genes <- all_genes[match(valid_genes, tolower(all_genes))]
  return(matched_genes)
}
output_dir <- "marker_gene_umap_plots_subset_10"
dir.create(output_dir, showWarnings = FALSE)

for (cell_type in names(markers)) {
  cat(paste0("Processing ", cell_type, "... \n"))
  
  genes <- markers[[cell_type]]
  valid_genes <- get_valid_genes(ips, genes)
  
  if (length(valid_genes) == 0) {
    cat("  No valid genes found.\n")
    next
  }
  
  plots <- lapply(valid_genes, function(gene) {
    p <- FeaturePlot(ips, features = gene, reduction = "umap", pt.size = 0.5) +
      ggtitle(gene) +
      theme(legend.position = "right")
    return(p)
  })
  
  # 拼图（横向）
  combined_plot <- wrap_plots(plots, nrow = 1)
  
  # 保存图像
  file_name <- paste0(cell_type, "_marker_genes.png")
  save_path <- file.path(output_dir, file_name)
  ggsave(save_path, plot = combined_plot, width = 4 * length(valid_genes), height = 4, dpi = 300)
}

###########################################################################################################################3
##另一个marker gene 画图方法

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


############################################################################################################################
####提出DAP相关簇，重新聚类分群注释
# 2. 查看原始 cluster 分布（确认）
table(ips@meta.data$seurat_clusters)
# 3. 提取目标 cluster：2, 4, 5, 6, 14, 19
target_clusters <- c(2, 4, 5, 6, 14, 19)
sub_ips <- subset(ips, subset = seurat_clusters %in% target_clusters)
new_rna_assay <- CreateAssayObject(counts = GetAssayData(sub_ips, assay = "RNA", layer = "counts"))
# 2. 将新创建的、索引正确的 assay 赋值回原对象。
#    这会替换掉那个损坏的 assay。
sub_ips[["RNA"]] <- new_rna_assay

# 3. 现在，尝试进行标准化等操作。
#    由于 assay 的内部结构已修复，这些操作应该可以正常进行。
sub_ips <- NormalizeData(sub_ips)
sub_ips <- FindVariableFeatures(sub_ips, selection.method = "vst", nfeatures = 2000)
sub_ips <- ScaleData(sub_ips)
sub_ips <- RunPCA(sub_ips)
p <- ElbowPlot(sub_ips)
ggsave("elbowplot_pca.png", plot = p, width = 8, height = 6, dpi = 300)
sub_ips <- FindNeighbors(sub_ips, dims = 1:10)
sub_ips <- FindClusters(sub_ips)  # 可尝试 0.4~1.0
# 8. 重新运行 UMAP（使用子集数据）
# 使用整合后的 PCA 的前20个维度来计算 UMAP
sub_ips <- RunUMAP(sub_ips, reduction = "integrated.rpca", dims = 1:10)

p <- DimPlot(sub_ips,reduction = "umap", group.by = "new_group", label = TRUE, label.size = 3,repel = TRUE)
# 保存为 PNG 或 PDF 文件
ggsave("umap_orig.ident.png", plot = p, width = 20, height = 10, dpi = 300)




##################################################################################################################################
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





############################################################################################################################
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

############################################################################################################################
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

###########################
table(Idents(ips_list))
cluster3 <- subset(ips_list, subset = RNA_snn_res.0.01 == 3)
dim(cluster3)
# 找出人源基因
human_genes <- grep("^GRCh38\\.109-", rownames(cluster3), value = TRUE)

# 筛选表达矩阵
cluster3_human <- cluster3[human_genes, ]
######################################################################################
# 提取下划线 "_" 之前的部分作为新列
ips@meta.data$sample_group <- sapply(strsplit(ips@meta.data$orig.ident, "_"), function(x) x[1])

library(Seurat)
library(readr)
library(loupeR)
all_barcodes <- read_tsv("barcodes.tsv", col_names = FALSE)
all_barcodes_vec <- as.character(all_barcodes[[1]])

set.seed(123)
selected_barcodes <- sample(all_barcodes_vec, size = ncol(ips), replace = FALSE)

sample_name <- ips@meta.data$sample_group

new_cell_name <- paste0(sample_name,"_",selected_barcodes)
names(new_cell_name) <- Cells(ips)
ips <- RenameCells(ips,new.names=new_cell_name)

create_loupe_from_seurat(ips, output_name = "umap_ips_3m_8m_integrated_rpca")

######################################################################################
####取特定细胞，查看其在umap上的分布
ref <- readRDS("/storage2/liuxiaodongLab/jiangjing/Projects/YutingFu/PD_YutingFu/20250718_ips_3_8M_all_umap/umap_ips_3m_8m_integrated_rpca_annocelltype_250722.rds")
ips <- readRDS("umap_ips_3m_8m_integrated_rpca.rds")
# 1. 提取 ips 中指定 cluster 的细胞 ID
selected_clusters <- c("0", "7", "9", "11", "17")
cells_of_interest <- rownames(ips@meta.data)[ips@meta.data$seurat_clusters %in% selected_clusters]

# 2. 检查 ref 是否为 Seurat 对象
if ("Seurat" %in% class(ref)) {
  # 处理细胞 ID 可能不匹配的情况（常见于整合后对象）
  # 检查是否有 '-1' 或来源前缀
  common_cells <- intersect(cells_of_interest, colnames(ref))
  
  # 如果没找到，尝试添加来源前缀（如 "ips_"）
  if (length(common_cells) == 0) {
    cells_of_interest_with_prefix <- paste0("ips_", cells_of_interest)
    common_cells <- intersect(cells_of_interest_with_prefix, colnames(ref))
  }
  
  if (length(common_cells) == 0) {
    stop("在 ref 中找不到 ips 的细胞，请检查细胞 ID 命名规则是否一致")
  }
  
  # 添加高亮列
  ref$highlight <- "Others"
  ref$highlight[common_cells] <- "Selected (ips clusters 0,7,9,11,17)"
  ref$highlight <- factor(ref$highlight, levels = c("Others", "Selected (ips clusters 0,7,9,11,17)"))
  
  # 3. 绘图并保存到文件
  library(Seurat)
  library(ggplot2)

  p <- DimPlot(ref, group.by = "highlight", cols = c("lightgray", "red"), label = FALSE, reduction = "umap", pt.size = 0.5) +
    ggtitle("UMAP: ips clusters 0,7,9,11,17 in ref") +
    theme(plot.title = element_text(size = 12, hjust = 0.5))

  # 保存图像
  ggsave("ips_selected_clusters_on_ref_umap.png", plot = p, width = 10, height = 8, dpi = 300, bg = "white")
  print("UMAP 图像已保存为 'ips_selected_clusters_on_ref_umap.png'")

} else {
  stop("ref 不是 Seurat 对象，请确认其类型。当前 class(ref): ", paste(class(ref), collapse = ", "))
}

######################################################################################
###注释
# 1. 定义注释列表
cluster_annotation <- list(
  "OPC" = c(10,23),
  "OLG" = c(20,22),
  "GLUL" = c(0,7,9,11,17),
  "Mature DA" = c(19),
  "5-HT" = c(15,24),
  "GABAL" = c(1,21),
  "DA progenitor-CHRNB3 high" = c(2),
  "DA progenitor-CHRNB3 low" = c(4,14),
  "Immature DA" = c(5,6),
  "Astrocyte progenitor cells" = c(12,13,25),
  "VLMC" = c(18),
  "Astrocyte" = c(3,8,16)
)

# 1. 设置当前 ident 为 RNA_snn_res.0.5
Idents(ips) <- "RNA_snn_res.0.5"

# 2. 提取 cluster ID（注意：它是字符型，需转为 numeric）
clusters <- as.numeric(as.character(Idents(ips)))

# 3. 创建 cell_type 向量，默认为 "Unknown"
cell_types <- rep("Unknown", length(clusters))

# 4. 根据你的 annotation list 进行映射
for (cell_type_name in names(cluster_annotation)) {
  target_clusters <- cluster_annotation[[cell_type_name]]
  # 匹配 cluster IDs 并赋值
  cell_types[clusters %in% target_clusters] <- cell_type_name
}

# 5. 将注释写入 meta.data
ips$cell_type <- factor(cell_types,
                                  levels = c(names(cluster_annotation), "Unknown"))

# 6. （可选）将注释后的 cell_type 设为当前 ident
Idents(ips) <- "cell_type"

# 7. 查看注释分布
table(Idents(ips))
ips$cell_type <- factor(ips$cell_type, levels = names(cluster_annotation))

p <- DimPlot(ips,reduction = "umap", group.by = c("RNA_snn_res.0.5","cell_type"), label = TRUE, label.size = 3,
    repel = TRUE) 

# 保存为 PNG 或 PDF 文件
ggsave("umap_0.5_celltype.png", plot = p, width = 30, height = 15, dpi = 300)

saveRDS(ips, "umap_ips_3m_8m_human_integrated_rpca_annocelltype_250728.rds")

