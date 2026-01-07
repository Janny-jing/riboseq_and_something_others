rm(list = ls())
# 加载必要的包 ----------------------------------------------------------------
library(ggplot2)
library(AUCell)
library(doParallel)  # 用于并行计算
library(ggrastr)     # 栅格化点图加速
library(magrittr)
library(dplyr)
# 初始化设置 -----------------------------------------------------------------
setwd("/storage/liuxiaodongLab/jiangjing/Projects/YutingFu/PD_YutingFu/scenic")
output_dir <- "optimized_regulon_plots"
dir.create(output_dir, showWarnings = FALSE)

# 数据预处理 -----------------------------------------------------------------
TF_gene <- read.csv("Step2_regulonTargetsInfo.csv", stringsAsFactors = FALSE)
umap_df <- read.table("umap_coordinates.txt")
expr_mat <- scale(read.table("expression_matrix_counts.txt"))

# AUCell计算核心部分（仅需运行一次） ------------------------------------------
regulon_list <- split(TF_gene$gene, TF_gene$TF)
regulon_list <- lapply(regulon_list, function(x) intersect(x, rownames(expr_mat)))
regulon_list <- regulon_list[lengths(regulon_list) > 0]

gene_rankings <- AUCell_buildRankings(expr_mat)
auc_matrix <- AUCell_calcAUC(regulon_list, gene_rankings)
auc_values <- t(getAUC(auc_matrix))
umap_df <- cbind(umap_df, auc_values)

# 预计算阈值（关键性能优化点） ------------------------------------------------
cell_thresholds <- AUCell_exploreThresholds(
  auc_matrix, 
  assignCells = TRUE
)

# 阈值提取函数（优化版）-----------------------------------------------------
extract_threshold <- function(tf, default = 0.2) {
  if (!tf %in% names(cell_thresholds)) return(default)
  
  th_obj <- tryCatch(
    cell_thresholds[[tf]]$aucThr$selected,
    error = function(e) NULL
  )
  
  if (is.numeric(th_obj)) return(unname(th_obj))
  if (is.list(th_obj) && "threshold" %in% names(th_obj)) return(th_obj$threshold)
  
  warning(paste("Using default threshold for", tf))
  default
}

# 可视化函数集 ---------------------------------------------------------------
#' 快速绘图模板
base_umap_plot <- function() {
  list(
    theme_void(),
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 10),
      legend.position = "bottom"
    ),
    guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)))
  )
}

#' 二值化活动图（优化渲染）
plot_binary_activity <- function(tf, data, threshold) {
    # 创建逻辑筛选条件
    active_mask <- data[[tf]] > threshold
  
    ggplot(data, aes(UMAP1, UMAP2)) +
    # 先绘制Inactive（底层）
      ggrastr::rasterise(
        geom_point(
          data = data[!active_mask, ], 
          color = "#F0F0F0",  # 更柔和的灰色
          size = 0.8,         # 缩小尺寸
          alpha = 0.5         # 更高的透明度
      ), 
      dpi = 300
    ) +
    # 再绘制Active（顶层）
    ggrastr::rasterise(
      geom_point(
        data = data[active_mask, ],
        color = "#4E79A7",  # 更明亮的蓝色
        size = 0.9,         # 增大尺寸
        alpha = 1         # 更不透明
      ),
      dpi = 300
    ) +
    labs(
      title = paste("Regulon:", tf),
      subtitle = sprintf("Threshold: %.2f", threshold)
    ) +
    base_umap_plot()
}

generate_safe_breaks <- function(values, 
                                 probs = seq(0, 1, 0.2),
                                 min_breaks = 3) {
  # 计算初始分位数
  breaks <- quantile(values, probs = probs, na.rm = TRUE, names = FALSE)
  
  # 去除重复值并排序
  unique_breaks <- sort(unique(breaks))
  
  # 处理特殊情况
  if (length(unique_breaks) < min_breaks) {
    if (all(values == 0)) {
      # 全零数据特殊处理
      return(c(-0.1, 0, 0.1))
    } else {
      # 使用等间距分箱作为备用
      return(seq(min(values), max(values), length.out = length(probs)))
    }
  }
  
  return(unique_breaks)
}

plot_gradient_activity <- function(tf, data) {
  if (all(data[[tf]] == 0)) return(NULL)
  
  # 生成安全断点（带扰动防止重复）
  safe_breaks <- generate_safe_breaks(data[[tf]] + 1e-8 * rnorm(nrow(data)))
  
  gradient_data <- data %>%
    arrange(.data[[tf]]) %>%
    mutate(
      alpha_layer = cut(
        .data[[tf]],
        breaks = safe_breaks,
        include.lowest = TRUE,
        labels = format(seq(0.3, 1, length.out = length(safe_breaks)-1), digits = 2)
      ) %>% 
        as.character() %>% 
        as.numeric()  # 双重转换确保数值类型
    )
  
  # 安全处理 NA (基础 R 方案)
  if (anyNA(gradient_data$alpha_layer)) {
    median_alpha <- median(gradient_data$alpha_layer, na.rm = TRUE)
    # 处理全 NA 极端情况
    if (is.na(median_alpha)) median_alpha <- 0.5
    gradient_data$alpha_layer <- ifelse(
      is.na(gradient_data$alpha_layer),
      median_alpha,
      gradient_data$alpha_layer
    )
  }
  
  ggplot(gradient_data, aes(UMAP1, UMAP2)) +
    rasterise(
      geom_point(
        aes(color = .data[[tf]], alpha = alpha_layer), 
        size = 0.8
      ), 
      dpi = 300
    ) +
    scale_color_gradientn(
      colours = c("#F0F0F0", "#90CAF9", "#42A5F5", "#0066CC"),
      values = scales::rescale(c(0, 0.3, 0.6, 1)),
      name = "AUC Value"
    ) +
    scale_alpha_continuous(
      range = c(0.3, 1),  # 强制alpha范围
      guide = "none"       # 隐藏透明度图例
    ) +
    labs(title = paste(tf, "Activity Gradient")) +
    base_umap_plot()
}


# 并行化批量绘图 -------------------------------------------------------------
registerDoParallel(cores = 6)  # 根据服务器核心数调整

foreach(tf = names(regulon_list), .packages = c("ggplot2", "ggrastr")) %dopar% {
  try({
    threshold <- extract_threshold(tf)
    
    # 二值图
    p1 <- plot_binary_activity(tf, umap_df, threshold)
    ggsave(
      file.path(output_dir, paste0(tf, "_binary.png")),
      p1, width = 6, height = 5, dpi = 300, bg = "white"
    )
    
    # 梯度图
    p2 <- plot_gradient_activity(tf, umap_df)
    ggsave(
      file.path(output_dir, paste0(tf, "_gradient.png")),
      p2, width = 6, height = 5, dpi = 300, bg = "white"
    )
  })
  
  NULL  # 防止不必要的返回
}
