# 单细胞数据整合分析脚本 - 修复样本重命名时机
# 版本：v3.2 - 调整重命名时机
# 日期：2024-06-17

# 加载必要的包
library(Seurat)
library(tidyverse)
library(patchwork)
library(Matrix)
library(future)
library(stringr)


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
  } 
  # 保留原列数检查（确保至少有1列）
  else if(ncol(features_df) < 2) {
    stop(paste("Invalid features file in", data.dir, "with", ncol(features_df), "columns"))
  }
  
  # 后续处理保持不变（使用V2列）...
  rat_gene_indices <- grep("^GRCh38\\.109_", features_df$V2)
  
  if(length(rat_gene_indices) == 0) {
    stop(paste("No rat genes found in", data.dir))
  }
  
  # 筛选矩阵和特征
  mat <- mat[rat_gene_indices, ]
  features_df <- features_df[rat_gene_indices, ]
  
  # 处理基因名 - 使用第二列
  features_df$V2 <- str_remove(features_df$V2, "^GRCh38\\.109_")
  
  # 处理基因名中的下划线
  features_df$V2 <- str_replace_all(features_df$V2, "_", "-")
  
  # 处理重复基因名
  if(any(duplicated(features_df$V2))) {
    features_df$V2 <- make.unique(features_df$V2, sep = "-")
  }
  
  # 设置行名 - 使用处理后的基因名
  rownames(mat) <- features_df$V2
  colnames(mat) <- barcodes
  
  # 记录信息
  cat("Rat genes found:", nrow(mat), "\n")
  
  return(mat)
}

# 1. 批量读取并处理h9样本 - 保持原始名称
h9_samples <- c("h9-1.5M_A", "h9-1.5M_B", "h9-1.5M_C", "h9-DAP")
h9_list <- list()
qc_list <- list()

output_dir <- "qc_violin_plots"
dir.create(output_dir, showWarnings = FALSE)
for (sample_name in h9_samples) {
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
  seu$dataset <- "New_H9"
  
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
  h9_list[[sample_name]] <- seu
  qc_list[[sample_name]] <- qc_stats
}

# 3. 处理现有数据集
existing_seu <- readRDS("new_umap_dim18_spread2.5_mindist0.7_celltypes_250418_orig.cellnames.rds")

# 删除P1A和P1B样本
if ("orig.ident" %in% colnames(existing_seu@meta.data)) {
  cells_to_keep <- !(existing_seu$orig.ident %in% c("P1A", "P1B"))
  existing_seu <- existing_seu[, cells_to_keep]
}

# 重命名元数据列
existing_seu$orig.batch <- existing_seu$batch 
existing_seu$orig.dataset <- existing_seu$dataset 

# 添加新元数据
existing_seu$dataset <- "Existing"

# 确保基因名格式一致
rownames(existing_seu) <- str_replace_all(rownames(existing_seu), "_", "-")

# 4. 准备全局整合数据
# 拆分现有数据集为多个样本对象
existing_ids <- unique(existing_seu$orig.ident)
existing_list <- SplitObject(existing_seu, split.by = "orig.ident")

# 合并所有样本
all_samples <- c(h9_list, existing_list)

merged_seurat <- merge(x = all_samples[[1]], 
                       y = all_samples[-1], 
                       add.cell.ids = sample_names)

merged_seurat <- NormalizeData(merged_seurat) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA(npcs = 50, verbose = FALSE)  # 关键：必须运行PCA

merged_seurat[["RNA"]] <- JoinLayers(merged_seurat[["RNA"]])
print(names(merged_seurat@reductions)) # 应显示"pca"

# 1. 创建新的分组列 dataset_new
merged_seurat$dataset_new <- as.character(merged_seurat$orig.ident)

# 2. 定义需要合并的组
rat1.5M_group <- c("h9-1.5M_A", "h9-1.5M_B", "h9-1.5M_C", 
                   "rat1.5M_A", "rat1.5M_B", "rat1.5M_C", 
                   "rat1.5M_D", "rat1.5M_E", "rat1.5M_F")

rat12M_BC_group <- c("rat12M_B", "rat12M_C")

# 3. 更新分组标识
merged_seurat$dataset_new[merged_seurat$orig.ident %in% rat1.5M_group] <- "rat1.5M"
merged_seurat$dataset_new[merged_seurat$orig.ident %in% rat12M_BC_group] <- "rat12M_BC"

saveRDS(merged_seurat, "global_before_integrated_20250618.rds")

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
sobj <- FindNeighbors(sobj, reduction = "harmony", dims = 1:50,k=20)
sobj <- FindClusters(sobj, resolution = 0.5)
sobj <- RunUMAP(sobj, dims = 1:50, 
                reduction = "harmony", 
                min.dist = 0.1, 
                spread=1,
                n.neighbors =20,
                seed.use = 42)
head(sobj@meta.data)
saveRDS(sobj, "global_integrated_harmony_20250618.rds")


png("umap_plot_human_celltypes_250418_dims50_mindist0.1_spread1_neighbors20.png", width = 4000, height = 1500, res = 300)
DimPlot(sobj, 
        group.by = "human_celltypes_250418", 
        label = TRUE,
        repel = TRUE,
        #split.by = "dataset",
        pt.size = 0.6)
dev.off()

######################################################################################################################


# 确保原始对象不被修改
sobj_backup <- sobj

# 创建参数组合网格
param_grid <- expand.grid(
  dims = seq(20, 31, 5),         
  k = seq(20, 31, 5),            
  resolution = seq(1, 2.1, 1), 
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
                          reduction = "harmony", 
                          dims = 1:d, 
                          k = k_val)
    
    sobj <- FindClusters(sobj, resolution = res)
    
    sobj <- RunUMAP(sobj, 
                    dims = 1:d,
                    reduction = "harmony",
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
              group.by = "human_celltypes_250418", 
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
    message(sprintf("参数组合失败 (d=%d,k=%d,res=%.1f,md=%.1f,sp=%.1f): %s",
                    d, k_val, res, md, sp, e$message))
  })
}

# 恢复原始对象
sobj <- sobj_backup








