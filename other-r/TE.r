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
  if (add_gene_symbols && "gene_id" %in% colnames(data)) {
    gene_symbols <- convert_ensembl_to_symbol(data$gene_id)
    data_with_names <- data.frame(
      gene_id = data$gene_id,
      gene_symbol = gene_symbols,
      data[, !colnames(data) %in% c("gene_id", "gene_symbol")],
      row.names = NULL,
      check.names = FALSE
    )
  } else if (add_gene_symbols) {
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

# 计算TPM函数
calculate_tpm <- function(counts, gene_lengths = NULL) {
  # 如果没有提供基因长度，假设所有基因长度相同（适用于标准化表达量）
  if (is.null(gene_lengths)) {
    # 简单标准化：TPM = (counts / sum(counts)) * 1e6
    tpm <- apply(counts, 2, function(x) {
      (x / sum(x, na.rm = TRUE)) * 1e6
    })
  } else {
    # 使用基因长度的标准化
    # 确保基因长度与counts的行名匹配
    gene_lengths_kb <- gene_lengths / 1000
    tpm <- apply(counts, 2, function(x) {
      rpk <- x / gene_lengths_kb
      (rpk / sum(rpk, na.rm = TRUE)) * 1e6
    })
  }
  return(tpm)
}

# 修改后的基因表达过滤函数 - 严格版本：所有样本都必须满足条件
filter_genes_by_expression_strict <- function(rna_counts, ribo_counts, 
                                              rna_threshold = 1,    # TPM阈值
                                              ribo_threshold = 10,  # RPF count阈值
                                              gene_lengths = NULL) {
  
  cat("开始基因表达过滤（严格模式）...\n")
  cat("筛选标准：\n")
  cat("1. Ribo-seq：每个基因在所有样本中count ≥", ribo_threshold, "\n")
  cat("2. RNA-seq：每个基因在所有样本中TPM ≥", rna_threshold, "\n")
  
  # 计算RNA-seq的TPM
  cat("计算RNA-seq TPM...\n")
  rna_tpm <- calculate_tpm(rna_counts, gene_lengths)
  
  # 获取所有样本名称
  all_samples <- colnames(rna_counts)
  
  # 检查每个基因在所有样本中是否满足条件
  cat("检查基因表达条件...\n")
  
  # 对于每个基因，检查所有样本
  keep_genes <- sapply(1:nrow(rna_counts), function(i) {
    # 检查RNA-seq：所有样本TPM ≥ rna_threshold
    rna_pass <- all(rna_tpm[i, ] >= rna_threshold, na.rm = TRUE)
    
    # 检查Ribo-seq：所有样本count ≥ ribo_threshold
    ribo_pass <- all(ribo_counts[i, ] >= ribo_threshold, na.rm = TRUE)
    
    # 必须同时满足RNA和Ribo的条件
    return(rna_pass & ribo_pass)
  })
  
  gene_names <- rownames(rna_counts)
  filtered_genes <- gene_names[keep_genes]
  
  # 打印详细的过滤统计
  cat("\n=== 过滤统计 ===\n")
  cat("总样本数:", length(all_samples), "\n")
  cat("原始基因数:", length(gene_names), "\n")
  cat("过滤后基因数:", length(filtered_genes), "\n")
  cat("过滤掉的基因数:", length(gene_names) - length(filtered_genes), "\n")
  cat("过滤比例:", round((1 - length(filtered_genes)/length(gene_names)) * 100, 2), "%\n")
  
  # 显示一些不符合条件的基因示例（如果有的话）
  if (length(filtered_genes) < length(gene_names)) {
    removed_genes <- gene_names[!keep_genes]
    cat("\n=== 不符合条件的基因示例（前10个）===\n")
    for (i in 1:min(10, length(removed_genes))) {
      gene <- removed_genes[i]
      rna_values <- round(rna_tpm[gene, ], 2)
      ribo_values <- ribo_counts[gene, ]
      
      # 找出不满足条件的样本
      rna_failed <- colnames(rna_tpm)[rna_tpm[gene, ] < rna_threshold]
      ribo_failed <- colnames(ribo_counts)[ribo_counts[gene, ] < ribo_threshold]
      
      cat(paste0(gene, ":\n"))
      cat("  RNA TPM: ", paste(names(rna_values), "=", rna_values, collapse=", "), "\n")
      cat("  Ribo count: ", paste(names(ribo_values), "=", ribo_values, collapse=", "), "\n")
      if (length(rna_failed) > 0) cat("  RNA不达标样本: ", paste(rna_failed, collapse=", "), "\n")
      if (length(ribo_failed) > 0) cat("  Ribo不达标样本: ", paste(ribo_failed, collapse=", "), "\n")
      cat("\n")
    }
  }
  
  return(filtered_genes)
}

# 去除批次效应函数
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

# 添加新函数：为基因列表添加样本counts数据
add_sample_counts_to_gene_list <- function(gene_list, rna_counts, ribo_counts, add_gene_symbols = TRUE) {
  # 确保gene_list是基因ID向量
  if (is.data.frame(gene_list)) {
    # 如果是数据框，提取gene_id列
    if ("gene_id" %in% colnames(gene_list)) {
      gene_ids <- gene_list$gene_id
    } else {
      gene_ids <- rownames(gene_list)
    }
  } else {
    gene_ids <- gene_list
  }
  
  # 只保留在两个counts矩阵中都存在的基因
  common_genes <- intersect(gene_ids, intersect(rownames(rna_counts), rownames(ribo_counts)))
  
  # 提取RNA-seq和Ribo-seq的counts数据
  rna_subset <- rna_counts[common_genes, , drop = FALSE]
  ribo_subset <- ribo_counts[common_genes, , drop = FALSE]
  
  # 重命名列以区分RNA和Ribo
  colnames(rna_subset) <- paste0("RNA_", colnames(rna_subset))
  colnames(ribo_subset) <- paste0("Ribo_", colnames(ribo_subset))
  
  # 合并RNA和Ribo数据
  combined_counts <- cbind(rna_subset, ribo_subset)
  
  # 添加基因ID列
  result <- data.frame(
    gene_id = common_genes,
    combined_counts,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  # 如果需要，添加基因符号
  if (add_gene_symbols) {
    result$gene_symbol <- convert_ensembl_to_symbol(common_genes)
    # 调整列顺序，使基因符号在第二列
    result <- result[, c("gene_id", "gene_symbol", setdiff(colnames(result), c("gene_id", "gene_symbol")))]
  }
  
  # 如果gene_list是数据框，合并其他信息
  if (is.data.frame(gene_list) && ncol(gene_list) > 1) {
    # 提取非counts的其他信息
    gene_info <- gene_list
    if ("gene_id" %in% colnames(gene_info)) {
      gene_info <- gene_info[gene_info$gene_id %in% common_genes, ]
      # 合并counts数据
      result <- merge(gene_info, result, by = "gene_id", all = FALSE)
    }
  }
  
  cat("成功为", length(common_genes), "个基因添加了样本counts数据\n")
  
  return(result)
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

# 保存原始共有基因数据
write.csv(rna, "rnaseq_common_genes.csv")
write.csv(ribo, "riboseq_common_genes.csv")

# 可选：如果有基因长度文件，可以读取
# gene_lengths <- read.table("gene_lengths.txt", header = T, row.names = 1)
# gene_lengths <- gene_lengths[common_genes, "length"]

### 2. 去除批次效应 ####
# 创建输出目录
dir.create("./plot", showWarnings = FALSE, recursive = TRUE)

cat("\n=== 去除批次效应 ===\n")
cat("去除Ribo-seq批次效应...\n")
ribo_processed <- sva_remove(ribo, "./plot/ribo")
cat("去除RNA-seq批次效应...\n")
rna_processed <- sva_remove(rna, "./plot/rna")

# 保存去批次后的数据
save_data_with_gene_names(ribo_processed, "./riboseq_nobatch_original.csv")
save_data_with_gene_names(rna_processed, "./rnaseq_nobatch_original.csv")

### 3. 基因表达过滤 ####
cat("\n=== 执行基因表达过滤（在去批次后）===\n")
cat("使用严格过滤标准：\n")
cat("- Ribo-seq: 每个基因在所有6个样本中count ≥ 10\n")
cat("- RNA-seq: 每个基因在所有6个样本中TPM ≥ 1\n")

filtered_genes <- filter_genes_by_expression_strict(
  rna_counts = rna_processed,
  ribo_counts = ribo_processed,
  rna_threshold = 1,     # RNA-seq: 所有样本TPM ≥ 1
  ribo_threshold = 10    # Ribo-seq: 所有样本count ≥ 10
  # gene_lengths = gene_lengths  # 如果有基因长度信息
)

# 应用过滤
rna_filtered <- rna_processed[filtered_genes, ]
ribo_filtered <- ribo_processed[filtered_genes, ]

cat("\n过滤后数据维度:\n")
cat("RNA-seq数据维度:", dim(rna_filtered), "\n")
cat("Ribo-seq数据维度:", dim(ribo_filtered), "\n")

# 验证过滤结果
cat("\n=== 验证过滤结果 ===\n")
# 随机检查一些基因，确保它们满足条件
if (length(filtered_genes) > 0) {
  # 计算过滤后RNA的TPM
  rna_tpm_filtered <- calculate_tpm(rna_filtered)
  
  # 随机选择5个基因进行验证
  set.seed(123)
  test_genes <- sample(filtered_genes, min(5, length(filtered_genes)))
  
  cat("随机检查", length(test_genes), "个过滤后的基因：\n")
  for (gene in test_genes) {
    rna_vals <- round(rna_tpm_filtered[gene, ], 3)
    ribo_vals <- ribo_filtered[gene, ]
    
    cat(paste0("\n基因: ", gene, "\n"))
    cat("RNA TPM值: ", paste(names(rna_vals), "=", rna_vals, collapse=", "), "\n")
    cat("最小值: ", min(rna_vals), " (≥1: ", min(rna_vals) >= 1, ")\n")
    cat("Ribo count值: ", paste(names(ribo_vals), "=", ribo_vals, collapse=", "), "\n")
    cat("最小值: ", min(ribo_vals), " (≥10: ", min(ribo_vals) >= 10, ")\n")
  }
  
  # 统计不符合条件的基因（如果有的话）
  all_passed <- TRUE
  for (gene in filtered_genes) {
    if (any(rna_tpm_filtered[gene, ] < 1) || any(ribo_filtered[gene, ] < 10)) {
      cat(paste0("\n警告：基因 ", gene, " 不满足条件！\n"))
      all_passed <- FALSE
    }
  }
  
  if (all_passed) {
    cat("\n✓ 所有过滤后的基因都满足条件！\n")
  }
}

# 保存过滤后的数据
save_data_with_gene_names(rna_filtered, "rnaseq_count_filtered_nobatch_strict.csv")
save_data_with_gene_names(ribo_filtered, "riboseq_count_filtered_nobatch_strict.csv")

### 4. 运行xTail分析 #####
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
cat("\n=== 运行xTail分析（使用去批次后严格过滤的数据）===\n")
xtail_result <- xtail_pipeline(ribo_filtered, rna_filtered)
saveRDS(xtail_result, "xtail_result_nobatch_filtered_strict.rds")

### 5. TE结果分析和可视化 #########
# 提取TE数据
R_data <- xtail_result$R$data

# 计算FDR (adjusted p-value)
R_data$pvalue.adjust <- p.adjust(R_data$pvalue_final, method = "BH")

# 按照新的筛选逻辑定义差异基因：
# 1. |log2FC_TE_final| ≥ 1 (fold change ≥ 1.5)
# 2. pvalue.adjust < 0.05
log2FC_threshold <- 1
FDR_threshold <- 0.05

R_data$Category <- ifelse(R_data$log2FC_TE_final >= log2FC_threshold & R_data$pvalue.adjust < FDR_threshold, "TE_Up",
                          ifelse(R_data$log2FC_TE_final <= -log2FC_threshold & R_data$pvalue.adjust < FDR_threshold, "TE_Down", 
                                 "Unchanged"))

R_data$Category <- factor(R_data$Category, levels = c("TE_Up", "TE_Down", "Unchanged"))

# 添加基因符号列
R_data$gene_id <- rownames(R_data)
R_data$gene_symbol <- convert_ensembl_to_symbol(rownames(R_data))

# 重新排列列顺序，让基因信息在前
column_order <- c("gene_id", "gene_symbol", "Category", 
                  "EV_log2TE", "orf3a_log2TE", "log2FC_TE_final", "pvalue_final", "pvalue.adjust",
                  setdiff(colnames(R_data), c("gene_id", "gene_symbol", "Category", 
                                              "EV_log2TE", "orf3a_log2TE", "log2FC_TE_final", "pvalue_final", "pvalue.adjust")))
R_data <- R_data[, column_order]

# 统计各分类基因数量
cat("\n=== 翻译效率变化基因统计 ===\n")
cat("筛选阈值: |log2FC| ≥", log2FC_threshold, "(fold change ≥ 1.5), pvalue.adjust <", FDR_threshold, "\n")
cat("TE_Up (翻译效率上调):", sum(R_data$Category == "TE_Up"), "\n")
cat("TE_Down (翻译效率下调):", sum(R_data$Category == "TE_Down"), "\n")
cat("Unchanged (无显著变化):", sum(R_data$Category == "Unchanged"), "\n")

# 设置颜色
te_colors <- c("TE_Up" = "#BF1D2D", "TE_Down" = "#283890", "Unchanged" = "grey")

# 绘制TE图 - 带新筛选标准
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
  geom_hline(yintercept = c(-log2FC_threshold, log2FC_threshold), linetype = 2, alpha = 0.5) +
  geom_vline(xintercept = c(-log2FC_threshold, log2FC_threshold), linetype = 2, alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, color = "black", alpha = 0.3) +
  xlab("log2 EV (RFP/mRNA)") +
  ylab("log2 orf3a (RFP/mRNA)") +
  labs(title = "Translational Efficiency Changes",
       subtitle = paste0("|log2FC| ≥ ", log2FC_threshold, " (fold change ≥ 1.5), pvalue.adjust < ", FDR_threshold, "\n",
                         "TE_Up: ", sum(R_data$Category == "TE_Up"), 
                         " genes | TE_Down: ", sum(R_data$Category == "TE_Down"), " genes")) +
  theme_classic() +
  theme(panel.border = element_rect(size = 1, fill = NA),
        axis.text = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(size = 10, hjust = 0.5))

### 6. 保存结果 ####
# 创建结果目录
dir.create("./results", showWarnings = FALSE)

# 保存图片
ggsave("./results/TE_plot_nobatch_filtered_strict.pdf", R_plot, height = 8, width = 9)
ggsave("./results/TE_plot_nobatch_filtered_strict.jpg", R_plot, height = 8, width = 9, dpi = 300)

# 保存筛选标准摘要
summary_text <- c(
  paste("分析流程: 先去除批次效应 → 基因表达过滤 → xTail分析"),
  paste("严格过滤标准:"),
  paste("  Ribo-seq: 每个基因在所有6个样本中count ≥ 10"),
  paste("  RNA-seq: 每个基因在所有6个样本中TPM ≥ 1"),
  paste("差异基因筛选标准: |log2FC| ≥", log2FC_threshold, "(fold change ≥ 1.5), pvalue.adjust <", FDR_threshold),
  paste("总分析基因数:", nrow(R_data)),
  paste("TE_Up基因数:", sum(R_data$Category == "TE_Up")),
  paste("TE_Down基因数:", sum(R_data$Category == "TE_Down")),
  paste("Unchanged基因数:", sum(R_data$Category == "Unchanged")),
  paste("分析时间:", Sys.time())
)

writeLines(summary_text, "./results/analysis_summary_nobatch_filtered_strict.txt")

# 保存完整的TE分析结果（带基因名）
save_data_with_gene_names(R_data, "./results/TE_analysis_results_nobatch_filtered_strict_with_symbols.csv", add_gene_symbols = FALSE)

# 分离上下调基因
te_up_genes <- R_data[R_data$Category == "TE_Up", ]
te_down_genes <- R_data[R_data$Category == "TE_Down", ]

# 为上下调基因添加样本counts数据
cat("\n=== 为上下调基因添加样本counts数据 ===\n")

# 添加RNA-seq和Ribo-seq的counts数据
te_up_with_counts <- add_sample_counts_to_gene_list(
  gene_list = te_up_genes,
  rna_counts = rna_filtered,
  ribo_counts = ribo_filtered,
  add_gene_symbols = FALSE  # 因为te_up_genes已经有gene_symbol列
)

te_down_with_counts <- add_sample_counts_to_gene_list(
  gene_list = te_down_genes,
  rna_counts = rna_filtered,
  ribo_counts = ribo_filtered,
  add_gene_symbols = FALSE  # 因为te_down_genes已经有gene_symbol列
)

# 保存带样本counts的上下调基因列表
write.csv(te_up_with_counts, "./results/TE_up_genes_with_sample_counts_strict.csv", row.names = FALSE)
write.csv(te_down_with_counts, "./results/TE_down_genes_with_sample_counts_strict.csv", row.names = FALSE)

cat("已保存带样本counts的上下调基因列表:\n")
cat("TE_up_genes_with_sample_counts_strict.csv (包含", nrow(te_up_with_counts), "个基因)\n")
cat("TE_down_genes_with_sample_counts_strict.csv (包含", nrow(te_down_with_counts), "个基因)\n")

# 同时保存原始版本（不带counts数据，用于兼容性）
write.csv(te_up_genes, "./results/TE_up_genes_nobatch_filtered_strict.csv", row.names = FALSE)
write.csv(te_down_genes, "./results/TE_down_genes_nobatch_filtered_strict.csv", row.names = FALSE)

# 保存中间数据
save_data_with_gene_names(ribo_filtered, "./results/riboseq_count_nobatch_filtered_strict_with_symbols.csv")
save_data_with_gene_names(rna_filtered, "./results/rnaseq_count_nobatch_filtered_strict_with_symbols.csv")

# 保存过滤后的TPM数据
rna_tpm_filtered <- calculate_tpm(rna_filtered)
save_data_with_gene_names(rna_tpm_filtered, "./results/rnaseq_tpm_nobatch_filtered_strict_with_symbols.csv")

# 额外：保存所有样本counts的合并数据（方便参考）
all_genes_with_counts <- add_sample_counts_to_gene_list(
  gene_list = rownames(rna_filtered),
  rna_counts = rna_filtered,
  ribo_counts = ribo_filtered,
  add_gene_symbols = TRUE
)

write.csv(all_genes_with_counts, "./results/all_filtered_genes_with_sample_counts_strict.csv", row.names = FALSE)
cat("已保存所有过滤基因的样本counts数据: all_filtered_genes_with_sample_counts_strict.csv\n")

# 验证最终的counts数据中没有小于10的值
cat("\n=== 验证最终数据 ===\n")
cat("检查最终数据中是否有Ribo count < 10的基因:\n")

# 查找是否有小于10的count
low_counts <- which(ribo_filtered < 10, arr.ind = TRUE)
if (length(low_counts) > 0) {
  cat("警告：发现", length(low_counts), "个小于10的Ribo count值！\n")
  
  # 显示具体信息
  for (i in 1:min(10, nrow(low_counts))) {
    row_idx <- low_counts[i, 1]
    col_idx <- low_counts[i, 2]
    gene <- rownames(ribo_filtered)[row_idx]
    sample <- colnames(ribo_filtered)[col_idx]
    count <- ribo_filtered[row_idx, col_idx]
    cat(paste0("基因: ", gene, ", 样本: ", sample, ", count: ", count, "\n"))
  }
} else {
  cat("✓ 所有Ribo count值都 ≥ 10！\n")
}

# 检查RNA TPM
rna_tpm_final <- calculate_tpm(rna_filtered)
low_tpm <- which(rna_tpm_final < 1, arr.ind = TRUE)
if (length(low_tpm) > 0) {
  cat("警告：发现", length(low_tpm), "个小于1的RNA TPM值！\n")
} else {
  cat("✓ 所有RNA TPM值都 ≥ 1！\n")
}

cat("\n=== 分析完成 ===\n")
cat("所有结果已保存到 ./results/ 目录\n")
cat("分析流程: 去批次效应 → 严格基因表达过滤 → xTail分析\n")
cat("\n重要输出文件:\n")
cat("1. TE_up_genes_with_sample_counts_strict.csv - 上调基因（带样本counts）\n")
cat("2. TE_down_genes_with_sample_counts_strict.csv - 下调基因（带样本counts）\n")
cat("3. TE_analysis_results_nobatch_filtered_strict_with_symbols.csv - 完整TE分析结果\n")
cat("4. all_filtered_genes_with_sample_counts_strict.csv - 所有过滤基因的样本counts\n")
cat("\n过滤标准验证:\n")
cat("- Ribo-seq: 所有样本count ≥ 10\n")
cat("- RNA-seq: 所有样本TPM ≥ 1\n")