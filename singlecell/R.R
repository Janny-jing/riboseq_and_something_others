rm(list = ls())
# 加载必要的包
library(igraph)
library(visNetwork)
library(ggplot2)
library(ggrepel)
library(umap)
library(AUCell)
setwd("/storage/liuxiaodongLab/jiangjing/Projects/YutingFu/PD_YutingFu/scenic")
# 读取数据（替换为实际数据路径）
TF_gene <- read.csv("Step2_regulonTargetsInfo.csv", stringsAsFactors = FALSE)
umap_df <- read.table("umap_coordinates.txt")
expr_mat <- read.table("expression_matrix_counts.txt")
expr_mat <- scale(expr_mat)
# 计算AUC评分（保持原有代码不变）
# --------------------------------------------------
regulon_list <- split(TF_gene$gene, TF_gene$TF)
regulon_list <- lapply(regulon_list, function(genes) intersect(genes, rownames(expr_mat)))
regulon_list <- regulon_list[sapply(regulon_list, length) > 0]

gene_rankings <- AUCell_buildRankings(expr_mat)
auc_matrix <- AUCell_calcAUC(regulon_list, gene_rankings)

# 合并数据（确保细胞顺序一致）
# --------------------------------------------------
auc_values <- t(getAUC(auc_matrix))
umap_df <- cbind(umap_df, auc_values)

write.csv(umap_df,file = "umap_AUCell.csv")
# # 批量可视化函数 -------------------------------------------------------------
# plot_regulon_umap <- function(tf_name, umap_data) {
#   ggplot(umap_data, aes(x = UMAP1, y = UMAP2, color = .data[[tf_name]])) +
#     geom_point(size = 0.5, alpha = 0.8) +
#     scale_color_gradient2(low = "grey90", mid = "orange", high = "red", midpoint = median(umap_data[[tf_name]])) +
#     ggtitle(paste(tf_name, "Regulon Activity")) +
#     theme_classic() +
#     theme(plot.title = element_text(hjust = 0.5))
# }
# 
# 
# if (!dir.exists("regulon_plots")) dir.create("regulon_plots")
# 
# # 循环绘制并保存
# for (tf in names(regulon_list)) {
#   p <- plot_regulon_umap(tf, umap_df)
#   ggsave(filename = file.path("regulon_plots", paste0(tf, "_activity.pdf")),
#          plot = p, width = 6, height = 5)
# }

plot_regulon_binary <- function(tf_name, 
                                umap_data,
                                default_threshold = 0.2) {
  
  # 获取阈值对象
  cell_thresholds <- AUCell_exploreThresholds(
    auc_matrix, 
    assignCells = TRUE
  )
  # 在调用 AUCell_exploreThresholds 后添加保存代码 --------------------------------------------------
  saveRDS(cell_thresholds, 
          file = "debug_cell_thresholds.rds",
          compress = FALSE)  # 关闭压缩确保可读性
  
  # 同时保存阈值摘要信息 (可选)
  sink("threshold_summary.txt")
  cat("=== Threshold Object Structure ===\n")
  str(cell_thresholds)
  cat("\n=== TF Coverage Check ===\n")
  cat("Total TFs processed:", length(cell_thresholds), "\n")
  cat("Missing TFs:", setdiff(names(regulon_list), names(cell_thresholds)), "\n")
  sink()
  
  extract_threshold <- function(tf) {
    if (!tf %in% names(cell_thresholds)) {
      warning(paste(tf, "not in cell_thresholds"))
      return(default_threshold)
    }
    
    # 层级访问修正
    if (is.null(cell_thresholds[[tf]]$aucThr)) {
      warning(paste(tf, "missing aucThr structure"))
      return(default_threshold)
    }
    
    th_obj <- cell_thresholds[[tf]]$aucThr$selected
    
    # 处理带有名称的数值
    if (is.numeric(th_obj) && !is.null(names(th_obj))) {
      return(unname(th_obj))  # 移除名称属性
    }
    
    # 其他类型处理保持不变...
    if (is.list(th_obj)) {
      if ("threshold" %in% names(th_obj)) {
        return(th_obj$threshold)
      } else if (!is.null(th_obj$y)) {
        return(th_obj$y[1])
      }
    }
    
    warning(paste("Unhandled threshold type for", tf))
    default_threshold
  }
  # 获取有效阈值
  tf_threshold <- tryCatch(
    {
      th <- extract_threshold(tf_name)
      if (!is.numeric(th)) stop("Non-numeric threshold")
      round(th, 3)  # 提前进行四舍五入
    },
    error = function(e) {
      warning(paste("Using default threshold for", tf_name, "due to:", e$message))
      default_threshold
    }
  )
  
  # 调试输出（可选）
  cat(paste0("[Debug] ", tf_name, 
             " | Threshold type: ", class(tf_threshold),
             " | Value: ", tf_threshold, "\n"))
  
  # 创建二值标签
  umap_data$activity <- ifelse(
    umap_data[[tf_name]] > tf_threshold,
    "Active",
    "Inactive"
  )
  
  # 可视化
  ggplot(umap_data, aes(x = UMAP1, y = UMAP2)) +
    geom_point(aes(color = activity), 
               size = 0.8, 
               alpha = 0.6) +
    scale_color_manual(
      values = c("Active" = "#377EB8", 
                 "Inactive" = "#E0E0E0"),
      drop = FALSE
    ) +
    labs(
      title = paste0("Regulon: ", tf_name),
      subtitle = paste("Activity Threshold:", tf_threshold)
    ) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
      legend.position = "bottom"
    )
}


# 执行可视化 --------------------------------------------------------
dir.create("regulon_plots_secure", showWarnings = FALSE)

for (tf in names(regulon_list)) {
  tryCatch(
    {
      p <- plot_regulon_binary(tf, umap_df)
      ggsave(
        file.path("regulon_plots_secure", paste0(tf, ".png")),
        p, width = 8, height = 6
      )
    },
    error = function(e) {
      message(paste("Failed to process", tf, ":", e$message))
    }
  )
}
