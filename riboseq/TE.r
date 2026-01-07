rm(list = ls())
setwd("Y:\\temp\\20251122-riboseq_new_3sample\\TE")
library(dplyr)
library(ggplot2)
library(ggrepel)
library(DESeq2)
library(reshape2)
library(xtail)
library(sva)
library(org.Hs.eg.db)  # 添加基因名转换包

### 0. 工具函数 ##########
mypca <- function(data, sampleTable){
  tmm <- data[apply(t(data), 2, sd) != 0, ]
  pca <- prcomp(t(tmm), center = TRUE, scale. = TRUE)
  df <- as.data.frame(pca$x)
  summ1 <- base::summary(pca)
  xlab1 <- paste0("PC1(", round(summ1$importance[2,1]*100, 2), "%)")
  ylab1 <- paste0("PC2(", round(summ1$importance[2,2]*100, 2), "%)")
  
  df$condition <- sampleTable$condition[sampleTable$sample %in% row.names(df)]
  
  p.pca1 <- ggplot(data = df, aes(x = .data$PC1, y = .data$PC2, shape = .data$condition)) +
    ggforce::geom_mark_ellipse(aes(fill = .data$condition), color = NA) +
    geom_point(size = 3, color = "gray23") +
    labs(x = xlab1, y = ylab1, color = "", shape = "", title = "") +
    scale_shape_manual(values = c(18:0)) +
    theme_bw() +
    theme(plot.title = element_text(hjust = .5, vjust = 0, size = 15),
          axis.text = element_text(size = 15), axis.title = element_text(size = 15),
          legend.text = element_text(size = 10), legend.title = element_text(size = 12),
          plot.margin = unit(c(1,1,1,1), "in"), panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    ggrepel::geom_text_repel(data = df, aes(x = .data$PC1, y = .data$PC2,
                                            label = rownames(df)), size = 5, max.overlaps = 100)
  
  return(p.pca1)
}

# 添加基因名转换函数
convert_ensembl_to_symbol <- function(ensembl_ids) {
  gene_symbols <- mapIds(org.Hs.eg.db,
                         keys = ensembl_ids,
                         column = "SYMBOL",
                         keytype = "ENSEMBL",
                         multiVals = "first")
  
  # 处理NA值，保留原始Ensembl ID
  result <- ifelse(is.na(gene_symbols), ensembl_ids, gene_symbols)
  
  # 统计转换成功率
  converted_count <- sum(!is.na(gene_symbols))
  cat("基因名转换统计:\n")
  cat("总基因数:", length(ensembl_ids), "\n")
  cat("成功转换:", converted_count, "\n")
  cat("转换成功率:", round(converted_count/length(ensembl_ids)*100, 2), "%\n")
  
  return(result)
}

# 添加带基因名的数据保存函数
save_data_with_gene_names <- function(data, filename, add_gene_symbols = TRUE) {
  if (add_gene_symbols) {
    gene_symbols <- convert_ensembl_to_symbol(rownames(data))
    data_with_names <- data.frame(
      gene_id = rownames(data),
      gene_symbol = gene_symbols,
      data,
      row.names = NULL,
      check.names = FALSE
    )
  } else {
    data_with_names <- data.frame(
      gene_id = rownames(data),
      data,
      row.names = NULL,
      check.names = FALSE
    )
  }
  write.csv(data_with_names, filename, row.names = FALSE)
  cat("已保存:", filename, "\n")
}

### 1. 导入数据 ####
# 读取RNA-seq数据
rna <- read.table("rnaseq_count.txt", header = T, row.names = 1)
# 修改列名以匹配脚本要求
colnames(rna) <- c("EV_1", "EV_2", "EV_3", "orf3a_1", "orf3a_2", "orf3a_3")

# 读取Ribo-seq数据
ribo <- read.table("riboseq_counts.txt", header = T, row.names = 1)
# 修改列名以匹配脚本要求
colnames(ribo) <- c("EV_1", "EV_2", "EV_3", "orf3a_1", "orf3a_2", "orf3a_3")

# 检查数据
cat("RNA-seq数据维度:", dim(rna), "\n")
cat("Ribo-seq数据维度:", dim(ribo), "\n")

# 找到共有的基因
common_genes <- intersect(rownames(rna), rownames(ribo))
cat("共有基因数量:", length(common_genes), "\n")

# 筛选共有的基因
rna <- rna[common_genes, ]
ribo <- ribo[common_genes, ]

### 2. 去除批次效应 ####
sva_remove <- function(data, pca.save = NULL){
  data <- round(data, 0)
  sampleTable <- data.frame(
    sample = colnames(data),
    condition = unlist(lapply(strsplit(colnames(data), "_"), function(x) x[[1]][1])),
    batch = unlist(lapply(strsplit(colnames(data), "_"), function(x) x[[2]][1]))
  )
  removed <- ComBat_seq(as.matrix(data),
                        group = sampleTable$condition,
                        batch = sampleTable$batch)
  if (!is.null(pca.save)) {
    pca_org <- mypca(data, sampleTable)
    pca <- mypca(removed, sampleTable)
    ggsave(paste0(pca.save, "_pca_original.jpg"), pca_org, height = 6, width = 6)
    ggsave(paste0(pca.save, "_pca_remove_batch.jpg"), pca, height = 6, width = 6)
  }
  return(removed)
}

# 创建输出目录
dir.create("./plot", showWarnings = FALSE, recursive = TRUE)

# 去除批次效应
cat("去除Ribo-seq批次效应...\n")
ribo <- sva_remove(ribo, "./plot/ribo")
cat("去除RNA-seq批次效应...\n")
rna <- sva_remove(rna, "./plot/rna")

### 3. 运行xTail分析 #####
xtail_pipeline <- function(ribo, rna, bins = 10000){
  ribo <- ribo[rowSums(ribo) != 0, ] %>% na.omit()
  ribo <- round(ribo, 0)
  
  rna <- rna[rowSums(rna) != 0, ] %>% na.omit()
  rna <- round(rna, 0)
  
  condition <- unlist(lapply(strsplit(colnames(ribo), "_"), function(x) x[[1]][1]))
  diffTE <- xtail(rna, ribo, condition, bins = bins)
  xtail_result <- resultsTable(diffTE, log2FC = T, log2R = T)
  plot_table <- resultsTable(diffTE)
  
  R_plot <- xtail::plotRs(object = diffTE)
  FC_plot <- xtail::plotFCs(object = diffTE)
  return(list(data = xtail_result, R = R_plot, FC = FC_plot))
}

# 运行xTail分析
cat("运行xTail分析...\n")
xtail_result <- xtail_pipeline(ribo, rna)
saveRDS(xtail_result, "xtail_result.rds")

### 4. TE结果分析和可视化 #########
# 提取TE数据
R_data <- xtail_result$R$data

# 使用xTail最终结果定义差异基因
R_data$Category <- ifelse(R_data$log2FC_TE_final >= 1 & R_data$pvalue_final < 0.05, "TE_Up",
                          ifelse(R_data$log2FC_TE_final <= -1 & R_data$pvalue_final < 0.05, "TE_Down", 
                                 "Unchanged"))

R_data$Category <- factor(R_data$Category, levels = c("TE_Up", "TE_Down", "Unchanged"))

# 添加基因符号列
R_data$gene_id <- rownames(R_data)
R_data$gene_symbol <- convert_ensembl_to_symbol(rownames(R_data))

# 重新排列列顺序，让基因信息在前
column_order <- c("gene_id", "gene_symbol", "Category", 
                  "EV_log2TE", "orf3a_log2TE", "log2FC_TE_final", "pvalue_final",
                  setdiff(colnames(R_data), c("gene_id", "gene_symbol", "Category", 
                                              "EV_log2TE", "orf3a_log2TE", "log2FC_TE_final", "pvalue_final")))
R_data <- R_data[, column_order]

# 统计各分类基因数量
cat("\n翻译效率变化基因统计:\n")
cat("TE_Up (翻译效率上调):", sum(R_data$Category == "TE_Up"), "\n")
cat("TE_Down (翻译效率下调):", sum(R_data$Category == "TE_Down"), "\n")
cat("Unchanged (无显著变化):", sum(R_data$Category == "Unchanged"), "\n")

# 设置颜色
te_colors <- c("TE_Up" = "#BF1D2D", "TE_Down" = "#283890", "Unchanged" = "grey")


# 绘制TE图 - 手动控制图层和图例
R_plot <- ggplot() +
  # 先画所有点（用于图例）
  geom_point(data = R_data,
             aes(x = EV_log2TE, y = orf3a_log2TE, color = Category), 
             size = 2, alpha = 1) +
  # 再覆盖一层灰色点（底层）
  geom_point(data = subset(R_data, Category == "Unchanged"),
             aes(x = EV_log2TE, y = orf3a_log2TE), 
             color = "grey", size = 2, alpha = 1) +
  # 最后覆盖差异基因点（最上层）
  geom_point(data = subset(R_data, Category != "Unchanged"),
             aes(x = EV_log2TE, y = orf3a_log2TE, color = Category), 
             size = 2, alpha = 1) +
  scale_color_manual(values = te_colors) +
  geom_hline(yintercept = c(-1, 1), linetype = 2, alpha = 0.5) +
  geom_vline(xintercept = c(-1, 1), linetype = 2, alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, color = "black", alpha = 0.3) +
  xlab("log2 EV (RFP/mRNA)") +
  ylab("log2 orf3a (RFP/mRNA)") +
  labs(title = "Translational Efficiency Changes",
       subtitle = paste0("TE_Up: ", sum(R_data$Category == "TE_Up"), 
                         " genes | TE_Down: ", sum(R_data$Category == "TE_Down"), " genes")) +
  theme_classic() +
  theme(panel.border = element_rect(size = 1, fill = NA),
        axis.text = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(size = 12, hjust = 0.5))
### 5. 保存结果 ####
# 保存图片
ggsave("TE_plot.pdf", R_plot, height = 7, width = 8)
ggsave("TE_plot.jpg", R_plot, height = 7, width = 8, dpi = 300)

# 保存完整的TE分析结果（带基因名）
save_data_with_gene_names(R_data, "TE_analysis_results_with_symbols.csv", add_gene_symbols = FALSE)

# 保存显著变化的基因列表（带基因名）
te_up_genes <- R_data[R_data$Category == "TE_Up", ]
te_down_genes <- R_data[R_data$Category == "TE_Down", ]

save_data_with_gene_names(te_up_genes, "TE_up_genes_with_symbols.csv", add_gene_symbols = FALSE)
save_data_with_gene_names(te_down_genes, "TE_down_genes_with_symbols.csv", add_gene_symbols = FALSE)

# 同时保存原始版本（兼容性）
write.csv(R_data, "TE_analysis_results.csv", row.names = FALSE)
write.csv(te_up_genes, "TE_up_genes.csv", row.names = FALSE)
write.csv(te_down_genes, "TE_down_genes.csv", row.names = FALSE)

# 保存去批次后的数据（带基因名）
save_data_with_gene_names(ribo, "riboseq_count_nobatch_with_symbols.csv")
save_data_with_gene_names(rna, "rnaseq_count_nobatch_with_symbols.csv")

# 同时保存原始版本
write.csv(ribo, "riboseq_count_nobatch.csv", row.names = TRUE)
write.csv(rna, "rnaseq_count_nobatch.csv", row.names = TRUE)

# 保存统计摘要
stats_summary <- data.frame(
  Category = c("TE_Up", "TE_Down", "Unchanged", "Total"),
  Count = c(sum(R_data$Category == "TE_Up"), 
            sum(R_data$Category == "TE_Down"), 
            sum(R_data$Category == "Unchanged"),
            nrow(R_data))
)
write.csv(stats_summary, "TE_analysis_summary.csv", row.names = FALSE)

# 显示显著基因的基因符号示例
if(nrow(te_up_genes) > 0) {
  cat("\n📈 翻译效率上调基因示例 (前10个):\n")
  print(head(te_up_genes[, c("gene_symbol", "log2FC_TE_final", "pvalue_final")], 10))
}

if(nrow(te_down_genes) > 0) {
  cat("\n📉 翻译效率下调基因示例 (前10个):\n")
  print(head(te_down_genes[, c("gene_symbol", "log2FC_TE_final", "pvalue_final")], 10))
}

cat("\n🎉 分析完成！生成的文件：\n")
cat("- xtail_result.rds: xTail分析结果\n")
cat("- TE_plot.pdf/jpg: 翻译效率变化图\n")
cat("- TE_analysis_results_with_symbols.csv: 完整的TE分析结果（带基因名）\n")
cat("- TE_up_genes_with_symbols.csv: 翻译效率上调基因列表（带基因名）\n")
cat("- TE_down_genes_with_symbols.csv: 翻译效率下调基因列表（带基因名）\n")
cat("- TE_analysis_summary.csv: 统计分析摘要\n")
cat("- riboseq_count_nobatch_with_symbols.csv: 去批次后的Ribo-seq数据（带基因名）\n")
cat("- rnaseq_count_nobatch_with_symbols.csv: 去批次后的RNA-seq数据（带基因名）\n")
cat("- plot/: PCA分析图目录\n")
cat("\n📝 同时生成了不带'_with_symbols'后缀的原始版本文件用于兼容性\n")

cat("\n✅ 所有分析已完成！\n")
