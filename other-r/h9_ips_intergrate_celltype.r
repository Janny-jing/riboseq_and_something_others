library(Seurat)
library(dplyr)
library(harmony)
library(SeuratWrappers)
library(tidyverse)
h9 <- readRDS("scanpy_UMAP_n_neighbors15_npcs22_mindist1.0_spread0.5.rds")
ips <- readRDS("ips_DAP_after_qc.rds")
h9_p1 <- subset(h9, subset = sample %in% c("P1A", "P1B"))

h9_p1$orig.ident <- h9_p1$sample


# 4. 准备全局整合数据
# 拆分现有数据集为多个样本对象
existing_list <- SplitObject(h9_p1, split.by = "sample")
ips <- SplitObject(ips, split.by = "orig.ident")

# 合并所有样本
all_samples <- c(ips, existing_list)

merged_seurat <- merge(x = all_samples[[1]], 
                       y = all_samples[-1])

merged_seurat <- NormalizeData(merged_seurat) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA(npcs = 50, verbose = FALSE)  # 关键：必须运行PCA

merged_seurat[["RNA"]] <- JoinLayers(merged_seurat[["RNA"]])
print(names(merged_seurat@reductions)) # 应显示"pca"

# 1. 创建新的分组列 dataset_new
merged_seurat$dataset_new <- as.character(merged_seurat$orig.ident)


merged_seurat[["RNA"]] <- split(merged_seurat[["RNA"]], f = merged_seurat$dataset_new)

obj <- IntegrateLayers(
  object = merged_seurat, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca",
  verbose = FALSE
)

sobj <- obj


sobj <- RunHarmony(
  merged_seurat,
  group.by.vars = "dataset_new",
  reduction = "pca",
  reduction.save = "harmony",
  assay.use = "RNA",
  project.dim = FALSE
)

sobj[["RNA"]] <- JoinLayers(sobj[["RNA"]])

# Further analysis steps
sobj <- FindNeighbors(sobj, reduction = "integrated.rpca", dims = 1:50,k=20)
sobj <- FindClusters(sobj, resolution = 0.5)
sobj <- RunUMAP(sobj, dims = 1:50, 
                reduction = "integrated.rpca", 
                min.dist = 0.1, 
                spread=1,
                n.neighbors =20,
                seed.use = 42)
head(sobj@meta.data)
saveRDS(sobj, "integrated_rpca.rds")


png("umap_plot_human_celltypes_250418_dims50_mindist0.5_spread0.5_neighbors20_rpca_1.png", width = 8000, height = 5000, res = 300)
DimPlot(sobj, 
        group.by = "human_celltypes_250630", 
        label = TRUE,
        repel = TRUE,
        split.by = "orig.ident",
        pt.size = 0.6)
dev.off()

######################################################################################################################


# 确保原始对象不被修改
sobj_backup <- sobj

# 创建参数组合网格
param_grid <- expand.grid(
  dims = seq(20, 31, 5),         
  k = seq(20, 31, 5),            
  resolution = seq(0.5, 2.1, 1), 
  min.dist = seq(0.5, 1, 0.2),  
  spread = seq(0.5, 3, 0.5)      
)

# 循环测试所有参数组合
for(i in 1:nrow(param_grid)) {
  # 从备份恢复对象
  sobj <- sobj_backup
  
  # 获取当前参数
  curr <- param_grid[i, ]
  d <- curr$dims
  k_val <- curr$k
  res <- curr$resolution
  md <- curr$min.dist
  sp <- curr$spread
  
  tryCatch({
    # 使用当前参数处理数据
    sobj <- FindNeighbors(sobj, 
                          reduction = "integrated.rpca", 
                          dims = 1:d, 
                          k = k_val)
    
    sobj <- FindClusters(sobj, resolution = res)
    
    sobj <- RunUMAP(sobj, 
                    dims = 1:d,
                    reduction = "integrated.rpca",
                    min.dist = md,
                    spread = sp,
                    n.neighbors = k_val,  # 与FindNeighbors的k保持一致
                    seed.use = 42)
    
    # 生成唯一文件名（包含所有参数）
    fname <- sprintf(
      "umap_d%02d_k%02d_res%.1f_md%.1f_sp%.1f.png",
      d, k_val, res, md, sp
    )
    param_title <- paste(
       sprintf("d:%d k:%d res:%.1f min.dist:%.1f spread:%.1f",
                d, k_val, res, md, sp)
    )
    
    # 保存高质量UMAP图
    png(fname, width = 4000, height = 1800, res = 300)  # 增加高度容纳标题
    
    print(
      DimPlot(sobj, 
              group.by = "human_celltypes_250630", 
              label = TRUE,
              repel = TRUE,
              split.by = "orig.ident",
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
    message(sprintf("参数组合失败 (d=%d,k=%d,res=%.1f,md=%.1f,sp=%.1f): %s",
                    d, k_val, res, md, sp, e$message))
  })
}

# 恢复原始对象
sobj <- sobj_backup

####################################################################################################################

library(Seurat)
library(readr)

all_barcodes <- read_tsv("barcodes.tsv", col_names = FALSE)
all_barcodes_vec <- as.character(all_barcodes[[1]])

set.seed(123)
selected_barcodes <- sample(all_barcodes_vec, size = ncol(sobj), replace = FALSE)

sample_name <- sobj@meta.data$dataset_new

new_cell_name <- paste0(sample_name,"_",selected_barcodes)
names(new_cell_name) <- Cells(sobj)
sobj <- RenameCells(sobj,new.names=new_cell_name)

create_loupe_from_seurat(sobj, output_name = "integrated_rpca")




##########################################################################################################
##两数剧间印射

library(Seurat)
library(SeuratDisk)  # 如果需要读写 .h5ad 文件
library(dplyr)
library(future)
plan(multisession)  # 设置一个执行计划

# 提高限制（例如允许传输 20GB）
options(future.globals.maxSize = 40 * 1024^3) 

# 加载数据
h9 <- readRDS("new_count_umap_dim18_spread2.5_mindist0.7_celltypes_250630_updated_original.rds")
ips <- readRDS("ips_DAP_after_qc.rds")

h9_p1 <- subset(h9, subset = sample %in% c("P1A", "P1B"))

h9_p1 <- NormalizeData(h9_p1) %>%
  FindVariableFeatures(nfeatures = 2000) %>%
  ScaleData() %>%
  RunPCA(npcs = 50, verbose = FALSE)


h9_p1[["RNA"]] <- split(h9_p1[["RNA"]], f = h9_p1$sample)
# Step 4: 使用 integrated.rpca 方法进行整合（如果需要）
h9_p1 <- IntegrateLayers(
  object = h9_p1,
  method = RPCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.rpca",
  verbose = FALSE
)
h9_p1[["RNA"]] <- JoinLayers(h9_p1[["RNA"]])

# Step 4: 构建邻居图并聚类
h9_p1 <- FindNeighbors(h9_p1, reduction = "integrated.rpca", dims = 1:50)
h9_p1 <- FindClusters(h9_p1, resolution = 1)
h9_p1 <- RunUMAP(h9_p1, reduction = "integrated.rpca", dims = 1:50, seed.use = 42)

# Step 5: query 数据预处理
ips <- NormalizeData(ips) %>%
  FindVariableFeatures(nfeatures = 2000) %>%
  ScaleData() %>%
  RunPCA(npcs = 50, verbose = FALSE)

# Elbow plot
p <- ElbowPlot(ips)
ggsave("elbowplot_pca.png", plot = p, width = 8, height = 6, dpi = 300)

# 聚类和 UMAP
ips <- FindNeighbors(ips, dims = 1:20) %>% FindClusters(resolution = 0.5)
#ips <- SetIdent(ips, value = "RNA_snn_res.1")
ips <- RunUMAP(ips, dims = 1:10 ,min.dist = 0.3, spread=1.5, n.neighbors = 30, seed.use = 42)

# Step 6: 寻找锚点并迁移注释
anchors <- FindTransferAnchors(
  reference = h9_p1,
  query = ips,
  dims = 1:20,
  reference.reduction = "integrated.rpca"
)

predictions <- TransferData(anchorset = anchors, refdata = h9_p1$human_celltypes_250630, dims = 1:20)
ips <- AddMetaData(ips, metadata = predictions)

#h9_p1 <- RunUMAP(h9_p1, dims = 1:30, reduction = "integrated.rpca", return.model = TRUE)
#ips <- MapQuery(anchorset = anchors, reference = h9_p1, query = ips,
#    refdata = list(celltype = "human_celltypes_250630"), reference.reduction = "integrated.rpca", reduction.model = "umap")



p1 <- DimPlot(h9_p1,reduction = "umap", group.by = "human_celltypes_250630", label = TRUE, label.size = 3,
    repel = TRUE) + NoLegend() + ggtitle("Reference annotations")
p2 <- DimPlot(ips,reduction = "ref.umap", group.by = "predicted.id", label = TRUE,
    label.size = 3, repel = TRUE) +  ggtitle("Query transferred labels")
# 组合图形
combined_plot <- p1 + p2

# 保存为 PNG 或 PDF 文件
ggsave("reference_and_query_umap.png", plot = combined_plot, width = 12, height = 6, dpi = 300)


####查看基因表达

marker_genes <- list(
  "Rgl1-like" = c("top2a", "mki67", "Lmx1a", "Apcdd1", "OTX2", "TTR"),
  "DA progenitor-FABP7" = c("FABP7", "TTR", "Lmx1a", "Apcdd1", "OTX2"),
  "DA progenitor-CHRNB3 high" = c("CHRNA6", "CHRNB3", "LMX1A"),
  "DA progenitor-CHRNB3 low" = c("DDC", "ASCL1", "LMX1A"),
  "DA progenitor-terminal" = c("LMX1A", "EN1", "NR4A2", "DDC")
)
all_marker_genes <- unique(unlist(marker_genes))

get_valid_genes <- function(sobj, gene_list) {
  assay_name <- DefaultAssay(sobj)
  
  # 获取当前 assay 的基因名
  all_genes <- rownames(GetAssayData(sobj, assay = assay_name, slot = "data"))
  
  # 忽略大小写匹配
  valid_genes <- intersect(tolower(gene_list), tolower(all_genes))
  
  # 返回原始大小写格式的匹配结果
  matched_genes <- all_genes[match(valid_genes, tolower(all_genes))]
  return(matched_genes)
}

# Step 1: 获取 reference 和 query 中真实存在的基因
ref_genes <- get_valid_genes(h9_p1, all_marker_genes)
query_genes <- get_valid_genes(ips, all_marker_genes)

# Step 2: 取交集，确保两个对象都有的基因才进行对比
genes_to_plot <- intersect(ref_genes, query_genes)
genes_to_plot <- sort(genes_to_plot)

# Step 3: 打印结果
print(paste0("Valid genes for plotting (", length(genes_to_plot), "): "))
print(genes_to_plot)


# Step 6: 绘图函数（安全版本）
plot_gene_expression <- function(sobj, reduction = "umap", genes, title_prefix = "") {
  plots <- lapply(genes, function(gene) {
    p <- FeaturePlot(
      object = sobj,
      features = gene,
      reduction = reduction,
      pt.size = 0.5,
      label = FALSE
    ) +
      ggtitle(paste0(title_prefix, ": ", gene)) +
      theme(legend.position = "right")
    
    return(p)
  })
  
  names(plots) <- genes
  return(plots)
}

# Step 7: 绘图
ref_plots <- plot_gene_expression(h9_p1, reduction = "umap", 
                                 genes = genes_to_plot, 
                                 title_prefix = "Reference")

query_plots <- plot_gene_expression(ips, reduction = "ref.umap", 
                                    genes = genes_to_plot, 
                                    title_prefix = "Query")

# Step 8: 组合对比图
combined_plots <- mapply(function(ref_p, query_p) {
  if (is.null(ref_p) || is.null(query_p)) return(NULL)
  ref_p + query_p
}, ref_plots, query_plots, SIMPLIFY = FALSE)

# 删除 NULL 元素
combined_plots <- Filter(Negate(is.null), combined_plots)

# Step 9: 保存图像
output_dir <- "gene_expression_comparison"
dir.create(output_dir, showWarnings = FALSE)

for (i in seq_along(combined_plots)) {
  gene_name <- names(combined_plots)[i]
  save_path <- file.path(output_dir, paste0(gene_name, "_comparison.png"))
  ggsave(save_path, plot = combined_plots[[i]], width = 12, height = 18, dpi = 300)
}



################################################################################


sobj <- SetIdent(sobj, value = RNA_snn_res.1)
markers <- FindAllMarkers(sobj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# 每个 cluster 取 top 5 marker 基因用于辅助人工判断
top_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 5) %>%
  ungroup() %>%
  split(.$cluster)

# 打印出来供你查看
print("Top marker genes for each cluster:")
print(top_markers)


# 定义标记基因字典
markers <- list(
  "Rgl1-like" = c("top2a", "mki67", "Lmx1a", "Apcdd1", "OTX2", "TTR"),
  "DA progenitor-FABP7" = c("FABP7", "TTR", "Lmx1a", "Apcdd1", "OTX2"),
  "DA progenitor-CHRNB3 high" = c("CHRNA6", "CHRNB3", "LMX1A","SCN1A-AS1", "SCN9A", "SCN5A", "SCNN1B", "SLC24A3"),
  "DA progenitor-CHRNB3 low" = c("DDC", "ASCL1", "LMX1A"),
  "DA progenitor-terminal" = c("LMX1A", "EN1", "NR4A2", "DDC"),
  "Immature DA" = c("KCNJ6", "TH", "SLC6A3", "SLC18A2", "PITX3", "NR4A2"),
  "Mature DA" = c("TH", "ALDH1A1", "NR4A2", "SLC18A2", "KCNJ6", "SLC6A3")
)





sobj <- RunUMAP(sobj, dims = 1:20, 
                n.neighbors =20,
                seed.use = 42)
p <- DimPlot(ips,reduction="umap")
ggsave("umap_by_new_sample.png", plot = p, width = 8, height = 6, dpi = 300)


##############################################################################################
library(Seurat)
library(SeuratDisk)  # 如果需要读写 .h5ad 文件
library(dplyr)
library(future)
plan(multisession)  # 设置一个执行计划

# 提高限制（例如允许传输 20GB）
options(future.globals.maxSize = 20 * 1024^3) 

ips <- readRDS("ips_DAP_after_qc.rds")

# Step 5: query 数据预处理
ips <- NormalizeData(ips) %>%
  FindVariableFeatures(nfeatures = 2000) %>%
  ScaleData() %>%
  RunPCA(npcs = 50, verbose = FALSE)


# 聚类和 UMAP
ips <- FindNeighbors(ips, dims = 1:30) %>% FindClusters(resolution = 1)
#ips <- SetIdent(ips, value = "RNA_snn_res.1")
ips <- RunUMAP(ips, dims = 1:30 ,seed.use = 42)

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
output_dir <- "marker_gene_umap_plots"
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