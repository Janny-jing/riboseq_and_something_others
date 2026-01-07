rm(list = ls())
setwd("Y:\\temp\\20250731_singlecell_twotwo")

library(ggplot2)
library(dplyr)
library(ggrepel)

de_genes$pvalue <- -log10(de_genes$p_val)
# 分别选出上调和下调前15个pvalue最小的基因
top_upregulated_by_pvalue <- head(de_genes[order(-de_genes$avg_log2FC), ], 10) # 上调基因按log2FC排序，然后按pvalue排序
top_downregulated_by_pvalue <- head(de_genes[order(de_genes$avg_log2FC), ], 10) # 下调基因按log2FC排序，然后按pvalue排序
# 合并两个数据框，用于在图中标注
to_label <- rbind(top_upregulated_by_pvalue, top_downregulated_by_pvalue)
# 标记需要标注的基因点
de_genes$highlight <- ifelse(rownames(de_genes) %in% rownames(to_label), "Highlighted", "Normal")

pdf(file="volcano.pdf",width=5.5,height = 6)
# 绘制火山图
ggplot(de_genes, aes(x = avg_log2FC, y = pvalue)) +
  # 上调基因用红色点表示
  geom_point(data = subset(de_genes, avg_log2FC > 1 & highlight == "Normal" & p_val < 0.05), aes(color = "Upregulated"), size = 2.5) +
  # 下调基因用蓝色点表示
  geom_point(data = subset(de_genes, avg_log2FC < -1 & highlight == "Normal"& p_val < 0.05), aes(color = "Downregulated"), size = 2.5) +
  # log2FC在-1到1之间的基因用灰色点表示
  geom_point(data = subset(de_genes, avg_log2FC >= -1 & avg_log2FC <= 1 | p_val > 0.05), aes(color = "Non-significant"), size = 2.5) +
  # 高亮前15个基因点为橙色
  geom_point(data = subset(de_genes, highlight == "Highlighted"), aes(color = "Highlighted"), size = 2.5) +
  # 标注上调和下调的前15个基因，使用橙色
  geom_text_repel(
    data = to_label,
    aes(label = rownames(to_label)),
    color = "black",
    size = 5,
    show.legend = FALSE,
    max.overlaps = Inf
  ) +
  xlim(-11,11)+
  scale_color_manual(values = c("Upregulated" = "#ED9090", "Downregulated" = "#055394", "Highlighted" = "#B6802C", "Non-significant" = "grey")) + # 设置颜色
  theme_minimal() +
  labs(title = "Top genes based on log2Fold change",
       x = "log2FC",
       y = "-log10(pvalue)",
       color = "Regulation") +
  theme_classic() +
  geom_vline(xintercept = 1, linetype = "dashed", color = "#999999", size = 0.5) +
  geom_vline(xintercept = -1, linetype = "dashed", color = "#999999", size = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#999999", size = 0.5) +
  # 调整主题以使标题居中并增大字体
  theme(
    plot.title = element_text(hjust = 0.5, size = 15),  # 居中并增大字体
    axis.text = element_text(size = 15),
    axis.title = element_text(size=15),
    legend.position = "none",
    axis.line = element_line(color = "black",size = 0.5),
    axis.ticks.length=unit(0.4, "cm"),
    axis.ticks=element_line(size = 0.5)
  )
dev.off()
#####################################################################


####################################################################
rm(list = ls())
setwd("Y:\\temp\\20250731_singlecell_twotwo")
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(cowplot)
library(readxl)
library(xlsx)
library(dplyr)
library(tibble)
crispr.data <- Read10X(data.dir = "Y:\\temp\\20250731_singlecell_twotwo\\01.cellrangercount\\outs\\filtered_feature_bc_matrix")
head(crispr.data[[1]],n=5)[,1:5]
head(crispr.data[[2]],n=5)[,1:5]
gex_matrix <- crispr.data[[1]]  # 基因表达
crispr_matrix <- crispr.data[[2]]  # sgRNA 计数

# 3. 创建 Seurat 对象（使用基因表达为主数据）
seu <- CreateSeuratObject(
  counts = gex_matrix,
  project = "CRISPR_screen",
  assay = "RNA"
)

# 4. 添加 CRISPR 数据作为第二个 Assay
crispr_assay <- CreateAssayObject(counts = crispr_matrix)
seu[["CRISPR"]] <- crispr_assay

# 5. QC（基于基因表达数据）
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
seu <- subset(seu,
              subset = nFeature_RNA > 200 &
                nFeature_RNA < 7500 &
                percent.mt < 10)

# 6. 标准化和降维（基于基因表达）
seu <- NormalizeData(seu) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA(features = VariableFeatures(.))

# 7. 构建邻居图 & UMAP（基于基因表达）
seu <- FindNeighbors(seu, dims = 1:30) %>%
  FindClusters() %>%
  RunUMAP(dims = 1:30)

crispr_call_file <- "Y:\\temp\\20250731_singlecell_twotwo\\01.cellrangercount\\outs\\crispr_analysis/protospacer_calls_per_cell.csv"
crispr_call <- read.csv(crispr_call_file, stringsAsFactors = FALSE)

# 设置行名为 cell barcode
rownames(crispr_call) <- crispr_call$cell_barcode
crispr_call$cell_barcode <- NULL  # 移除原列

# 找到与 Seurat 对象匹配的细胞
common_cells <- intersect(colnames(seu), rownames(crispr_call))

feature_call <- crispr_call[match(colnames(seu), rownames(crispr_call)), "feature_call"]
num_features <- crispr_call[match(colnames(seu), rownames(crispr_call)), "num_features"]

# 4. 添加到 Seurat
seu[["feature_call"]] <- feature_call
seu[["num_features"]] <- num_features

# 查看 metadata
head(seu@meta.data[, c("feature_call", "num_features")])

# 查看分类分布
table(seu[["feature_call"]])
DimPlot(seu, reduction = "umap", group.by = "feature_call", repel = TRUE, label = TRUE)
saveRDS(seu,"20250824_twotwocombine_crispr.rds")

library(reticulate)
seu <- readRDS("20250824_twotwocombine_crispr.rds")
seu[["RNA"]] <- as(seu[["RNA"]],"Assay")  #这步运行不成功有可能是seurat包不是v5，这个针对v5包
sceasy::convertFormat(seu,from="seurat",to="anndata",outFile="20250824_twotwocombine_crispr.h5ad")

fc <- seu@meta.data[["feature_call"]]

# 查看结构
str(fc)  # 看是 factor 还是 character

# 强制转为 character，保留 NA
fc <- as.character(fc)

# 现在 fc 是真正的字符向量，NA 是真正的 NA
head(fc, 20)

# 定义 untransduced 类别
untransduced_cats <- c("", "No confident call", "No guide molecules", NA, "None")
fc[fc %in% untransduced_cats] <- "untransduced"

# 更新到 Seurat
seu[["feature_call"]] <- fc
table(seu[["feature_call"]], useNA = "no")
# 3. 排除 barcode2|barcode3|barcode4 的 8 个细胞
seu_filtered <- subset(seu, subset = feature_call != "barcode2|barcode3|barcode4")

table(seu_filtered[["feature_call"]], useNA = "no")

# 4. 创建分组变量
# 提取 feature_call
fc <- seu_filtered@meta.data[["feature_call"]]

# 判断是否为 multiple
is_multiple <- grepl("\\|", fc) & (fc != "untransduced")

# 添加新列
seu_filtered[["perturbation_type"]] <- ifelse(fc == "untransduced", "untransduced",
                                              ifelse(is_multiple, "multiple", "single"))
table(seu_filtered[["perturbation_type"]], useNA = "no")


# 5. 创建输出目录
dir.create("crispr_de_results", showWarnings = FALSE)

# 6. 函数：保存差异表达结果 + 火山图
save_de_results <- function(deg, group1, group2, name) {
  # 保存表格
  file <- file.path("crispr_de_results", paste0(name, "_deg.csv"))
  write.csv(deg, file, row.names = TRUE)
  
  # 火山图
  deg_plot <- deg %>%
    rownames_to_column("gene") %>%
    mutate(group = ifelse(p_val_adj < 0.05 & abs(avg_log2FC) > 0.25,
                          ifelse(avg_log2FC > 0, "Up", "Down"),
                          "NS"))
  
  p <- ggplot(deg_plot, aes(x = avg_log2FC, y = -log10(p_val_adj), color = group)) +
    geom_point(alpha = 0.8, size = 1) +
    scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "gray")) +
    theme_bw() +
    labs(title = paste("DE:", group1, "vs", group2),
         x = "log2 Fold Change",
         y = "-log10(p-value)",
         color = "Significance") +
    theme(legend.position = "right")
  
  ggsave(file.path("crispr_de_results", paste0(name, "_volcano.png")), 
         plot = p, width = 6, height = 5, dpi = 300)
  
  return(p)
}

# ===================================================
# 比较 1: 单个 vs untransduced（分开做）
# ===================================================
for (bc in c("barcode2", "barcode3", "barcode4")) {
  cat("Running DE:", bc, "vs untransduced\n")
  
  sub <- subset(seu_filtered, subset = feature_call %in% c(bc, "untransduced"))
    # 关键：设置 Idents 为 feature_call（确保只有两个 group）
  Idents(sub) <- "feature_call"
  deg <- FindMarkers(sub, ident.1 = bc, ident.2 = "untransduced",
                     test.use = "wilcox", logfc.threshold = 0.25, min.pct = 0.1, 
                     return.thresh = 0.05)
  
  save_de_results(deg, bc, "untransduced", paste0(bc, "_vs_untransduced"))
}

# ===================================================
# 比较 2: 单个 vs 单个（两两比较）
# ===================================================
single_pairs <- combn(c("barcode2", "barcode3", "barcode4"), 2, simplify = FALSE)
for (pair in single_pairs) {
  bc1 <- pair[1]
  bc2 <- pair[2]
  cat("Running DE:", bc1, "vs", bc2, "\n")
  
  sub <- subset(seu_filtered, subset = feature_call %in% c(bc1, bc2))
  Idents(sub) <- "feature_call"
  
  deg <- FindMarkers(sub, ident.1 = bc2, ident.2 = bc1,
                     test.use = "wilcox", logfc.threshold = 0.25, min.pct = 0.1, 
                     return.thresh = 0.05)
  
  save_de_results(deg, bc1, bc2, paste0(bc1, "_vs_", bc2))
}

# ===================================================
# 比较 3: 单个 vs 多个（按主扰动分组）
# ===================================================
# 定义每种单扰动对应的多重组合
multi_map <- list(
  barcode2 = c("barcode2|barcode3", "barcode2|barcode4"),
  barcode3 = c("barcode2|barcode3", "barcode3|barcode4"),
  barcode4 = c("barcode2|barcode4", "barcode3|barcode4")
)
# 1. 修复 feature_call 结构
if (is.data.frame(seu_filtered@meta.data[["feature_call"]])) {
  cat("⚠️ 检测到 feature_call 为 data.frame，正在修复...\n")
  fc_vec <- seu_filtered[["feature_call"]][["feature_call"]]  # 提取内部向量
  seu_filtered[["feature_call"]] <- fc_vec
}

# 3. 获取真实存在的 barcodes
existing_barcodes <- unique(seu_filtered@meta.data[["feature_call"]])

# 4. 重新运行你的循环
for (bc in c("barcode2", "barcode3", "barcode4")) {
  multi_groups <- multi_map[[bc]]
  
  if (!any(multi_groups %in% existing_barcodes)) {
    cat("Warning:", bc, "has no corresponding multiple groups. Skipping.\n")
    next
  }
  
  cat("Running DE:", bc, "single vs multiple\n")
  
  sub <- subset(seu_filtered, subset = feature_call %in% c(bc, multi_groups))
  
  group <- ifelse(sub[["feature_call"]] == bc, "single", "multiple")
  Idents(sub) <- factor(group, levels = c("single", "multiple"))
  
  deg <- FindMarkers(sub, ident.1 = "multiple", ident.2 = "single",
                     test.use = "wilcox", logfc.threshold = 0.25, min.pct = 0.1, 
                     return.thresh = 0.05)
  
  save_de_results(deg, "single", "multiple", paste0(bc, "_single_vs_multiple"))
}
# ===================================================
# 比较 4: 多个 vs untransduced（修复版）
# ===================================================
cat("Running DE: multiple vs untransduced\n")

# 获取所有多重感染组
multiple_groups <- unique(subset(seu_filtered@meta.data, perturbation_type == "multiple")$feature_call)

# 子集
sub <- subset(seu_filtered, subset = feature_call %in% c(multiple_groups, "untransduced"))

# 临时分组
group <- ifelse(sub[["feature_call"]] == "untransduced", "untransduced", "multiple")
Idents(sub) <- factor(group, levels = c("untransduced", "multiple"))

deg <- FindMarkers(sub, ident.1 = "multiple", ident.2 = "untransduced",
                   test.use = "wilcox", logfc.threshold = 0.25, min.pct = 0.1, 
                   return.thresh = 0.05)

save_de_results(deg, "multiple", "untransduced", "multiple_vs_untransduced")

# ===================================================
# 比较 5: 多个 vs untransduced（分别）
# ===================================================
cat("Running DE: each multiple group vs untransduced\n")

# 获取所有多重感染组
multiple_barcodes <- unique(subset(seu_filtered@meta.data, perturbation_type == "multiple")$feature_call)

# 逐个比较
for (multi_bc in multiple_barcodes) {
  cat("Running DE:", multi_bc, "vs untransduced\n")
  
  # 子集：只保留该多重组合 + untransduced
  sub <- subset(seu_filtered, subset = feature_call %in% c(multi_bc, "untransduced"))
  
  # 设置分组
  Idents(sub) <- "feature_call"
  
  # 差异表达
  deg <- FindMarkers(
    sub,
    ident.1 = multi_bc,
    ident.2 = "untransduced",
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    return.thresh = 0.05
  )
  
  # 保存结果
  save_de_results(deg, multi_bc, "untransduced", paste0(gsub("\\|", "_", multi_bc), "_vs_untransduced"))
}
# ===================================================
# 比较 6: 多重 vs 单重（每个多重组合 vs 其组分单重）
# ===================================================
cat("Running DE: multiple vs single components\n")

# 定义每个多重组合的“组分”
multi_to_single <- list(
  "barcode2|barcode3" = c("barcode2", "barcode3"),
  "barcode2|barcode4" = c("barcode2", "barcode4"),
  "barcode3|barcode4" = c("barcode3", "barcode4")
)

# 遍历每个多重组合
for (multi_bc in names(multi_to_single)) {
  singles <- multi_to_single[[multi_bc]]
  
  # 检查该多重组合是否存在
  if (!(multi_bc %in% seu_filtered@meta.data[["feature_call"]])) {
    cat("Warning: ", multi_bc, " not found. Skipping.\n")
    next
  }
  
  # 检查两个单重是否都存在
  available_singles <- singles[singles %in% seu_filtered@meta.data[["feature_call"]]]
  if (length(available_singles) == 0) {
    cat("Warning: no single components for ", multi_bc, " found. Skipping.\n")
    next
  }
  
  # 对每个单重组分进行比较
  for (single_bc in available_singles) {
    cat("Running DE:", multi_bc, "vs", single_bc, "\n")
    
    # 子集：多重 + 单重
    sub <- subset(seu_filtered, subset = feature_call %in% c(multi_bc, single_bc))
    
    # 设置分组
    Idents(sub) <- "feature_call"
    
    # 差异表达：多重 vs 单重
    deg <- FindMarkers(
      sub,
      ident.1 = multi_bc,
      ident.2 = single_bc,
      test.use = "wilcox",
      logfc.threshold = 0.25,
      min.pct = 0.1,
      return.thresh = 0.05
    )
    
    # 保存结果
    name <- paste0(gsub("\\|", "_", multi_bc), "_vs_", single_bc)
    save_de_results(deg, multi_bc, single_bc, name)
  }
}
# ===================================================
# 比较 7: 多重 vs 多个（每个多重组合 vs 多个）
# ===================================================
dir.create("crispr_de_results/multiple_vs_multiple", showWarnings = FALSE, recursive = TRUE)
dir.create("crispr_de_results/supplement", showWarnings = FALSE, recursive = TRUE)

# 定义比较组合
multi_pairs <- list(
  "barcode2_barcode3_vs_barcode2_barcode4" = c("barcode2|barcode3", "barcode2|barcode4"),
  "barcode2_barcode3_vs_barcode3_barcode4" = c("barcode2|barcode3", "barcode3|barcode4"),
  "barcode2_barcode4_vs_barcode3_barcode4" = c("barcode2|barcode4", "barcode3|barcode4")
)

# 循环比较
for (name in names(multi_pairs)) {
  pair <- multi_pairs[[name]]
  bc1 <- pair[1]
  bc2 <- pair[2]
  
  cat("Running DE:", bc1, "vs", bc2, "\n")
  
  # 子集
  sub <- subset(seu_filtered, subset = feature_call %in% pair)
  
  # ✅ 修复 1: 确保 feature_call 是字符
  sub@meta.data[["feature_call"]] <- as.character(sub@meta.data[["feature_call"]])
  
  # ✅ 修复 2: 创建临时分组（最安全）
  group <- ifelse(sub[["feature_call"]] == bc1, bc1, bc2)
  Idents(sub) <- factor(group, levels = c(bc2, bc1))  # ident.2 = bc2, ident.1 = bc1
  
  # 差异表达
  deg <- FindMarkers(
    sub,
    ident.1 = bc1,
    ident.2 = bc2,
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    return.thresh = 0.05
  )
  
  # 保存 DEG 结果
  deg_file <- file.path("crispr_de_results/multiple_vs_multiple", paste0(name, "_deg.csv"))
  write.csv(deg, deg_file, row.names = TRUE)
  
  # ================================
  # GO/KEGG 富集分析
  # ================================
  run_enrichment(deg_file, name = name, output_dir = "crispr_de_results/enrichment")
}



# ===================================================
# 7. 综合热图（top DEGs from barcode2/3 vs untransduced）
# ===================================================
top_genes <- c()
for (bc in c("barcode2", "barcode3","barcode4")) { 
  deg_file <- file.path("crispr_de_results", paste0(bc, "_vs_untransduced_deg.csv"))
  if (file.exists(deg_file)) {
    deg <- read.csv(deg_file, row.names = 1)
    top_up <- rownames(deg[deg$p_val_adj < 0.05 & deg$avg_log2FC > 0.5, ])[1:5]
    top_down <- rownames(deg[deg$p_val_adj < 0.05 & deg$avg_log2FC < -0.5, ])[1:5]
    top_genes <- c(top_genes, top_up, top_down)
  }
}
top_genes <- unique(top_genes)
top_genes <- top_genes[top_genes %in% rownames(seu_filtered)]

if (length(top_genes) > 0) {
  pheatmap(
    GetAssayData(seu_filtered, assay = "RNA", slot = "data")[top_genes, ],
    annotation_col = seu_filtered@meta.data[colnames(seu_filtered), c("feature_call", "perturbation_type")],
    show_rownames = TRUE,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    filename = "crispr_de_results/top_genes_heatmap.png",
    width = 12,
    height = 8
  )
}

# ===================================================
# 8. 保存处理后的对象
# ===================================================
saveRDS(seu_filtered, "crispr_de_results/seu_filtered.rds")

cat("✅ 所有分析完成！结果已保存至 'crispr_de_results' 目录\n")

# ===================================================
# 函数：对 DEG 表做 GO/KEGG 富集分析
# 输入：deg_csv 文件路径
# 输出：GO_BP, GO_MF, KEGG 结果 + 气泡图
# ===================================================
library(clusterProfiler)
library(org.Hs.eg.db)  # 如果是小鼠，用 org.Mm.eg.db
library(enrichplot)
run_enrichment <- function(deg_file, name, output_dir = "crispr_de_results/enrichment") {
  # 创建输出目录
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # 1. 检查文件是否存在
  if (!file.exists(deg_file)) {
    cat("⚠️ 文件不存在:", deg_file, "\n")
    return(NULL)
  }
  
  # 2. 读取 DEG 表
  deg <- tryCatch({
    read.csv(deg_file, row.names = 1)
  }, error = function(e) {
    cat("⚠️ 读取失败:", deg_file, "\n")
    NULL
  })
  
  if (is.null(deg)) return(NULL)
  
  # 3. 筛选显著上调基因（ident.1 中高表达）
  sig_up <- deg %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::filter(p_val_adj < 0.05 & avg_log2FC > 0) %>%
    dplyr::arrange(p_val_adj)
  
  if (nrow(sig_up) == 0) {
    cat("⚠️", name, ": 无显著上调基因（p_val_adj < 0.05 & avg_log2FC > 0），跳过富集分析\n")
    return(NULL)
  }
  
  cat("✅", name, ": 发现", nrow(sig_up), "个显著上调基因，开始富集分析\n")
  
  # 4. 基因名转换：SYMBOL → ENTREZID
  # 假设你的基因是 SYMBOL（如 TP53, CDKN1A）
  gene_vector <- tryCatch({
    bitr(
      geneID = sig_up$gene,
      fromType = "SYMBOL",        # ✅ 根据你的数据调整：SYMBOL / ENSEMBL / ENTREZID
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db        # 如果是小鼠，改为 org.Mm.eg.db
    )
  }, error = function(e) {
    cat("⚠️ bitr 错误:", e$message, "\n")
    NULL
  })
  
  if (is.null(gene_vector) || nrow(gene_vector) == 0) {
    cat("⚠️", name, ": 基因转换失败（可能基因名类型错误），跳过\n")
    return(NULL)
  }
  
  # 5. GO 富集分析（生物过程 BP）
  go_bp <- tryCatch({
    enrichGO(
      gene = gene_vector$ENTREZID,
      OrgDb = org.Hs.eg.db,
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.05,
      minGSSize = 5,      # 允许较小通路
      maxGSSize = 500
    )
  }, error = function(e) {
    cat("⚠️ GO 富集失败:", e$message, "\n")
    NULL
  })
  
  # 6. KEGG 富集分析
  kegg <- tryCatch({
    enrichKEGG(
      gene = gene_vector$ENTREZID,
      organism = "human",  # 如果是小鼠，改为 "mouse"
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05
    )
  }, error = function(e) {
    cat("⚠️ KEGG 富集失败:", e$message, "\n")
    NULL
  })
  
  # 7. 保存结果
  if (!is.null(go_bp) && length(go_bp) > 0) {
    go_df <- as.data.frame(go_bp)
    write.csv(go_df, file.path(output_dir, paste0(name, "_GO_BP.csv")), row.names = FALSE)
    
    # 画图
    tryCatch({
      p1 <- dotplot(go_bp, showCategory = 10) + ggtitle(paste("GO BP Enrichment\n", name))
      ggsave(file.path(output_dir, paste0(name, "_GO_BP.png")), plot = p1, width = 8, height = 6, dpi = 300)
    }, error = function(e) {
      cat("⚠️ GO 画图失败:", e$message, "\n")
    })
  }
  
  if (!is.null(kegg) && length(kegg) > 0) {
    kegg_df <- as.data.frame(kegg)
    write.csv(kegg_df, file.path(output_dir, paste0(name, "_KEGG.csv")), row.names = FALSE)
    
    # 画图
    tryCatch({
      p2 <- dotplot(kegg, showCategory = 10) + ggtitle(paste("KEGG Pathway Enrichment\n", name))
      ggsave(file.path(output_dir, paste0(name, "_KEGG.png")), plot = p2, width = 8, height = 6, dpi = 300)
    }, error = function(e) {
      cat("⚠️ KEGG 画图失败:", e$message, "\n")
    })
  }
  
  # 8. 返回结果
  cat("✅", name, "富集分析完成\n")
  return(list(
    go_bp = go_bp,
    kegg = kegg,
    n_up_genes = nrow(sig_up),
    n_mapped = nrow(gene_vector)
  ))
}



# ===================================================
# 批量对所有差异分析结果做富集
# ===================================================
enrichment_results <- list()

# 获取所有 _deg.csv 文件
deg_files <- list.files("crispr_de_results/", pattern = "_deg.csv$", full.names = TRUE)

for (file in deg_files) {
  # 提取名字（如 barcode2_vs_untransduced）
  name <- tools::file_path_sans_ext(basename(file))
  
  # 跳过非 DEG 文件
  if (!grepl("_vs_", name)) next
  
  # 运行富集
  result <- run_enrichment(file, name)
  if (!is.null(result)) {
    enrichment_results[[name]] <- result
  }
}

seu <- readRDS("crispr_de_results/seu_filtered.rds")
# 正常着色
DimPlot(seu, group.by = "perturbation_type", label = TRUE)
DimPlot(seu, group.by = "feature_call", label = TRUE, repel = TRUE)
# 画 UMAP 并标注 cluster
DimPlot(seu,  label = TRUE)

# 检查每个 cluster 中的 perturbation_type 分布
table(seu@meta.data[["perturbation_type"]], seu@meta.data[["seurat_clusters"]])


# 只保留转染细胞
seu_transfected <- subset(seu, subset = perturbation_type %in% "multiple")

# 重新计算 HVGs 和 UMAP（仅用于探索）
seu_transfected <- FindVariableFeatures(seu_transfected, nfeatures = 500)
seu_transfected <- RunPCA(seu_transfected, features = VariableFeatures(seu_transfected))
seu_transfected <- RunUMAP(seu_transfected, dims = 1:10)

# 看是否能分开
DimPlot(seu_transfected, group.by = "feature_call", label = FALSE)

#################################################################3
# ===================================================
# 分析: barcode2 vs barcode3 在 multiple 细胞中的 UMAP 分布
# ===================================================

# 1. 读取 protospacer_calls_per_cell.csv
crispr_call_file <- "Y:/temp/20250731_singlecell_twotwo/01.cellrangercount/outs/crispr_analysis/protospacer_calls_per_cell.csv"
crispr_call <- read.csv(crispr_call_file, stringsAsFactors = FALSE)

# 设置行名为 cell barcode
rownames(crispr_call) <- crispr_call$cell_barcode
crispr_call$cell_barcode <- NULL

# 1. 提取 barcode2|barcode3 的细胞名
cells_multi <- rownames(crispr_call)[crispr_call$feature_call == "barcode2|barcode3"]

# 2. 检查这些细胞是否在 seu 中
cells_in_seu <- cells_multi[cells_multi %in% colnames(seu)]
cat("Found", length(cells_in_seu), "barcode2|barcode3 cells in Seurat object\n")

# 3. 提取 num_umis 并拆分
umis <- crispr_call[cells_in_seu, "num_umis"]
umis_split <- strsplit(umis, "\\|")

# 4. 转为数值
umis_b2 <- sapply(umis_split, function(x) as.numeric(x[1]))
umis_b3 <- sapply(umis_split, function(x) as.numeric(x[2]))

# 5. 创建分组
group <- ifelse(umis_b2 > umis_b3, "b2>b3", 
                ifelse(umis_b3 > umis_b2, "b3>b2", "equal"))

# 6. 创建子集 Seurat 对象
seu_multi <- subset(seu, cells = cells_in_seu)

# 7. 添加分组信息
seu_multi[["b2_vs_b3"]] <- factor(group, levels = c("b2>b3", "b3>b2", "equal"))

# 8. 重新计算 UMAP（可选，但推荐）
# 因为你只关心这 251 个细胞的内部结构
seu_multi <- RunPCA(seu_multi, features = VariableFeatures(seu))  # 使用原 HVGs
seu_multi <- RunUMAP(seu_multi, dims = 1:20)

# 9. 画图
DimPlot(seu_multi, group.by = "b2_vs_b3", label = TRUE, repel = TRUE) +
  ggtitle("UMAP: barcode2 vs barcode3 UMI in barcode2|barcode3 Cells")


seu_multi[["group"]] <- factor(group, levels = c("b2>b3", "b3>b2", "equal"))
seu_multi[["umi_b2"]] <- umis_b2
seu_multi[["umi_b3"]] <- umis_b3

# 6. 重新计算 UMAP（可选）
seu_multi <- RunPCA(seu_multi, features = VariableFeatures(seu))
seu_multi <- RunUMAP(seu_multi, dims = 1:20)

# 7. 画图：分组 + UMI 强度
library(cowplot)

# 主图：分组
p1 <- DimPlot(seu_multi, group.by = "group", label = TRUE, repel = TRUE) +
  ggtitle("Group: b2>b3 vs b3>b2")

# 分组内 barcode2 UMI
p2 <- FeaturePlot(seu_multi, features = "umi_b2", 
                  label = TRUE, 
                  cols = c("lightgrey", "red")) +
  ggtitle("UMI of barcode2 in each group")

# 分组内 barcode3 UMI
p3 <- FeaturePlot(seu_multi, features = "umi_b3", 
                  label = TRUE, 
                  cols = c("lightgrey", "blue")) +
  ggtitle("UMI of barcode3 in each group")

# 拼图
plot_grid(p1, p2, p3, nrow = 1)





