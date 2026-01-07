#!/usr/bin/env Rscript
# ==============================================================================
# 精简版Ribo-seq分析脚本 - riboWaltz版本
# ==============================================================================

# 清除环境
rm(list = ls())
gc()

# 设置工作目录（请修改为您的实际路径）
# setwd("/path/to/your/project")

# ==============================================================================
# 第一部分：安装和加载包
# ==============================================================================
cat("========================================\n")
cat("Ribo-seq分析流程 - riboWaltz精简版\n")
cat("========================================\n\n")

cat("1. 检查并安装必要的R包...\n")

# 检查并安装devtools
if (!require("devtools", quietly = TRUE)) {
  install.packages("devtools", repos = "https://cloud.r-project.org")
}

# 检查并安装riboWaltz
if (!require("riboWaltz", quietly = TRUE)) {
  cat("正在安装riboWaltz...\n")
  devtools::install_github("LabTranslationalArchitectomics/riboWaltz", dependencies = TRUE)
}

# 加载必要包
library(riboWaltz)
library(ggplot2)
library(data.table)

cat("✓ R包加载完成\n\n")

# ==============================================================================
# 第二部分：参数设置
# ==============================================================================
cat("2. 设置分析参数...\n")

# 文件路径
paths <- list(
  gtf_file = "gencode.v48.annotation.gtf",
  bam_folder = ".",
  output_dir = "riboseq_results"
)

# 创建输出目录
dir.create(paths$output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(paths$output_dir, "plots"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(paths$output_dir, "tables"), showWarnings = FALSE, recursive = TRUE)

cat("✓ 输出目录:", paths$output_dir, "\n\n")

# ==============================================================================
# 第三部分：创建注释文件
# ==============================================================================
cat("3. 创建注释文件...\n")

if (file.exists(paths$gtf_file)) {
  tryCatch({
    annotation <- create_annotation(gtfpath = paths$gtf_file)
    cat("✓ 注释文件创建成功\n")
  }, error = function(e) {
    cat("✗ 创建注释文件失败，使用示例数据\n")
    data(mm81cdna)
    annotation <- mm81cdna
  })
} else {
  cat("警告: GTF文件不存在，使用示例数据\n")
  data(mm81cdna)
  annotation <- mm81cdna
}

# ==============================================================================
# 第四部分：读取BAM文件
# ==============================================================================
cat("4. 读取BAM文件...\n")

bam_files <- list.files(paths$bam_folder, pattern = "\\.bam$", full.names = TRUE)
bam_files <- bam_files[!grepl("\\.bai$", bam_files)]

cat("   找到", length(bam_files), "个BAM文件\n")

if (length(bam_files) > 0) {
  # 标准化样本名
  sample_names <- gsub("\\.bam$", "", basename(bam_files))
  sample_names <- gsub("_filter$", "", sample_names)
  sample_names <- gsub("_sorted$", "", sample_names)
  sample_names <- gsub("_aligned$", "", sample_names)
  
  # 标准化命名：统一ORF3a为大写
  sample_names <- gsub("orf3a", "ORF3a", sample_names, ignore.case = TRUE)
  
  cat("   样本名:", paste(sample_names, collapse = ", "), "\n")
  
  tryCatch({
    # 读取BAM文件
    reads_list <- bamtolist(
      bamfolder = paths$bam_folder,
      annotation = annotation
    )
    
    # 重命名为标准化名称
    names(reads_list) <- sample_names
    
    cat("✓ BAM文件读取成功\n")
    cat("   样本数量:", length(reads_list), "\n")
    
    saveRDS(reads_list, file.path(paths$output_dir, "tables", "reads_list.rds"))
    
  }, error = function(e) {
    cat("✗ 读取BAM文件失败:", e$message, "\n")
    cat("使用示例数据继续...\n")
    data(reads_list)
  })
} else {
  cat("警告: 未找到BAM文件，使用示例数据\n")
  data(reads_list)
}

# ==============================================================================
# 第五部分：P-site定位
# ==============================================================================
cat("\n5. 计算P-site偏移...\n")

tryCatch({
  psite_offset <- psite(
    reads_list,
    flanking = 6,
    extremity = "auto"
  )
  
  cat("✓ P-site偏移计算完成\n")
  
  fwrite(psite_offset,
         file.path(paths$output_dir, "tables", "psite_offsets.csv"))
  
}, error = function(e) {
  cat("✗ P-site计算失败:", e$message, "\n")
})

# ==============================================================================
# 第六部分：添加P-site信息
# ==============================================================================
cat("\n6. 添加P-site信息...\n")

tryCatch({
  reads_psite_list <- psite_info(reads_list, psite_offset)
  cat("✓ P-site信息添加完成\n")
  
  saveRDS(reads_psite_list, 
          file.path(paths$output_dir, "tables", "reads_psite_list.rds"))
  
}, error = function(e) {
  cat("✗ 添加P-site信息失败\n")
  reads_psite_list <- reads_list
})

# ==============================================================================
# 第七部分：生成核心图表
# ==============================================================================
cat("\n7. 生成核心分析图表...\n")
cat("----------------------------------------\n")

# 7.1 读长分布图
cat("7.1 生成读长分布图...\n")
tryCatch({
  all_samples <- names(reads_list)
  
  rlen_plot <- rlength_distr(
    reads_list,
    sample = all_samples,
    multisamples = "independent",
    plot_style = "dodge",
    cl = 95,
    colour = scales::brewer_pal(palette = "Set3")(length(all_samples))
  )
  
  if(!is.null(rlen_plot$plot)) {
    ggsave(file.path(paths$output_dir, "plots", "read_length_distribution.pdf"),
           rlen_plot$plot,
           width = 12, height = 8, dpi = 300)
    cat("✓ 读长分布图已保存\n")
  }
}, error = function(e) {
  cat("✗ 读长分布图生成失败\n")
})

# 7.2 区域分布柱状图（重点优化）
cat("7.2 生成区域分布柱状图...\n")
tryCatch({
  # 创建样本列表
  sample_list <- as.list(names(reads_psite_list))
  names(sample_list) <- names(reads_psite_list)
  
  # 生成区域分布图
  region_plot <- region_psite(
    reads_psite_list,
    annotation,
    sample = sample_list,
    multisamples = "average",
    plot_style = "stack",
    cl = 95
  )
  
  if(!is.null(region_plot$plot)) {
    # 优化图表外观
    optimized_plot <- region_plot$plot +
      # 使用明亮的颜色
      scale_fill_manual(
        values = c(
          "5utr" = "#FF6B6B",    # 亮红色
          "cds" = "#4ECDC4",     # 亮青色
          "3utr" = "#FFD166"     # 亮黄色
        ),
        name = "Region",
        labels = c("5' UTR", "CDS", "3' UTR")
      ) +
      # 优化标题和标签
      labs(
        title = "P-site Distribution Across Transcript Regions",
        subtitle = "Ribo-seq analysis of translation activity",
        x = "Sample",
        y = "Percentage (%)"
      ) +
      # 优化主题
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40"),
        axis.text.x = element_text(
          angle = 45, 
          hjust = 1,
          size = 12,
          face = "bold",
          color = "black"
        ),
        axis.text.y = element_text(size = 11),
        axis.title = element_text(size = 13, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 11),
        panel.grid.major = element_line(color = "gray90"),
        panel.grid.minor = element_blank(),
        plot.margin = margin(20, 20, 20, 20)
      )
    
    # 保存优化后的图
    ggsave(file.path(paths$output_dir, "plots", "region_distribution.pdf"),
           optimized_plot,
           width = 14, height = 8, dpi = 300)
    
    # 也保存PNG格式
    ggsave(file.path(paths$output_dir, "plots", "region_distribution.png"),
           optimized_plot,
           width = 14, height = 8, dpi = 300)
    
    # 保存数据
    fwrite(region_plot$plot_dt,
           file.path(paths$output_dir, "tables", "region_distribution.csv"))
    
    cat("✓ 区域分布图已保存（使用明亮颜色）\n")
  }
}, error = function(e) {
  cat("✗ 区域分布图生成失败:", e$message, "\n")
})

# 7.3 Frame分布图
cat("7.3 生成Frame分布图...\n")
tryCatch({
  # 自动分组：EV为Control，ORF3a为实验组
  control_samples <- names(reads_psite_list)[grepl("^EV", names(reads_psite_list))]
  orf3a_samples <- names(reads_psite_list)[grepl("ORF3a", names(reads_psite_list))]
  
  sample_info <- list()
  if(length(control_samples) > 0) sample_info[["Control"]] <- control_samples
  if(length(orf3a_samples) > 0) sample_info[["ORF3a"]] <- orf3a_samples
  
  if(length(sample_info) > 0) {
    frame_plot <- frame_psite(
      reads_psite_list,
      annotation,
      sample = sample_info,
      multisamples = "average",
      plot_style = "facet",
      region = "cds",
      colour = c("0" = "#E41A1C", "1" = "#377EB8", "2" = "#4DAF4A")
    )
    
    if(!is.null(frame_plot$plot)) {
      # 优化Frame分布图
      optimized_frame_plot <- frame_plot$plot +
        scale_fill_manual(
          values = c("0" = "#FF6B6B", "1" = "#4ECDC4", "2" = "#FFD166"),
          labels = c("Frame 0", "Frame 1", "Frame 2"),
          name = "Reading Frame"
        ) +
        labs(
          title = "Frame Distribution in CDS Regions",
          x = "Sample Group",
          y = "Percentage (%)"
        ) +
        theme_minimal(base_size = 14) +
        theme(
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
          strip.text = element_text(size = 12, face = "bold")
        )
      
      ggsave(file.path(paths$output_dir, "plots", "frame_distribution.pdf"),
             optimized_frame_plot,
             width = 12, height = 8, dpi = 300)
      
      fwrite(frame_plot$plot_dt,
             file.path(paths$output_dir, "tables", "frame_distribution.csv"))
      
      cat("✓ Frame分布图已保存\n")
    }
  }
}, error = function(e) {
  cat("✗ Frame分布图生成失败\n")
})

# 7.4 Metagene图
cat("7.4 生成Metagene图...\n")
tryCatch({
  # 如果分组成功，生成分组metagene图
  if(exists("sample_info") && length(sample_info) > 0) {
    metagene_plot <- metaprofile_psite(
      reads_psite_list,
      annotation,
      sample = sample_info,
      multisamples = "average",
      plot_style = "split",
      utr5l = 30,
      cdsl = 60,
      utr3l = 30
    )
    
    # 保存各组的图
    for(group_name in names(sample_info)) {
      plot_name <- paste0("plot_", group_name)
      if(!is.null(metagene_plot[[plot_name]])) {
        ggsave(file.path(paths$output_dir, "plots", 
                        paste0("metagene_", tolower(group_name), ".pdf")),
               metagene_plot[[plot_name]],
               width = 10, height = 6, dpi = 300)
      }
    }
    
    # 也生成组合图
    metagene_combined <- metaprofile_psite(
      reads_psite_list,
      annotation,
      sample = sample_info,
      multisamples = "average",
      plot_style = "overlap",
      utr5l = 30,
      cdsl = 60,
      utr3l = 30,
      colour = c("Control" = "#FF6B6B", "ORF3a" = "#4ECDC4")
    )
    
    if(!is.null(metagene_combined$plot)) {
      ggsave(file.path(paths$output_dir, "plots", "metagene_combined.pdf"),
             metagene_combined$plot,
             width = 10, height = 6, dpi = 300)
    }
    
    cat("✓ Metagene图已保存\n")
  }
}, error = function(e) {
  cat("✗ Metagene图生成失败\n")
})

# ==============================================================================
# 第八部分：生成汇总统计
# ==============================================================================
cat("\n8. 生成汇总统计...\n")

tryCatch({
  summary_data <- data.frame()
  
  for(sample_name in names(reads_psite_list)) {
    sample_data <- reads_psite_list[[sample_name]]
    
    if(!is.null(sample_data) && "psite_region" %in% colnames(sample_data)) {
      total_reads <- nrow(sample_data)
      
      # 区域分布
      region_counts <- table(sample_data$psite_region, useNA = "ifany")
      
      # Frame分布（CDS区域）
      cds_data <- sample_data[sample_data$psite_region == "cds", ]
      frame0_pct <- if(nrow(cds_data) > 0) {
        sum(cds_data$psite_from_start %% 3 == 0, na.rm = TRUE) / nrow(cds_data) * 100
      } else 0
      
      summary_data <- rbind(summary_data, data.frame(
        Sample = sample_name,
        Total_Reads = total_reads,
        Reads_in_5UTR = ifelse("5utr" %in% names(region_counts), 
                               region_counts["5utr"], 0),
        Reads_in_CDS = ifelse("cds" %in% names(region_counts), 
                              region_counts["cds"], 0),
        Reads_in_3UTR = ifelse("3utr" %in% names(region_counts), 
                               region_counts["3utr"], 0),
        Frame0_Percentage = frame0_pct,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  if(nrow(summary_data) > 0) {
    # 计算百分比
    summary_data$Pct_5UTR <- round(summary_data$Reads_in_5UTR / summary_data$Total_Reads * 100, 1)
    summary_data$Pct_CDS <- round(summary_data$Reads_in_CDS / summary_data$Total_Reads * 100, 1)
    summary_data$Pct_3UTR <- round(summary_data$Reads_in_3UTR / summary_data$Total_Reads * 100, 1)
    
    # 保存
    fwrite(summary_data,
           file.path(paths$output_dir, "tables", "summary_statistics.csv"))
    
    cat("✓ 汇总统计已保存\n")
    cat("\n样本统计摘要:\n")
    print(summary_data[, c("Sample", "Total_Reads", "Pct_CDS", "Frame0_Percentage")])
  }
}, error = function(e) {
  cat("✗ 汇总统计生成失败\n")
})




##################################################distribution plot###############################################################
library(ggplot2)
library(dplyr)

# 从 region_plot$plot_dt 提取数据
plot_data <- region_plot$plot_dt

# 过滤掉"RNAs"样本
plot_data <- plot_data %>%
  filter(x != "RNAs")

# 简化样本名称
plot_data <- plot_data %>%
  mutate(
    sample_short = gsub("_ribo_Aligned.toTranscriptome.out", "", x),
    # 进一步简化：只保留样本标识符
    sample_simple = case_when(
      grepl("EV", sample_short) ~ gsub("EV-", "EV", sample_short),
      grepl("ORF3a", sample_short) ~ gsub("ORF3a-", "ORF3a", sample_short),
      TRUE ~ sample_short
    )
  )

# 设置新的颜色映射（按你的要求）
region_colors <- c(
  "5' UTR" = "#6A0DAD",     # 深紫色
  "CDS" = "#8B0000",        # 深红色
  "3' UTR" = "#006400"      # 深绿色
)

# 创建绘图
p <- ggplot(plot_data, aes(x = sample_simple, y = y, fill = z)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  scale_fill_manual(
    values = region_colors,
    name = "Transcript Region",
    labels = c("5' UTR", "CDS", "3' UTR")
  ) +
  labs(
    title = "P-site Distribution Across Transcript Regions",
    subtitle = "Ribo-seq analysis of translation activity",
    x = "Sample",
    y = "Percentage of P-sites (%)",
    fill = "Region"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(face = "bold", size = 12),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10),
    legend.position = "right",
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    plot.margin = margin(20, 20, 20, 20)  # 增加边距
  ) +
  # 添加百分比标签（可选）
  geom_text(
    aes(label = sprintf("%.1f%%", y)),
    position = position_stack(vjust = 0.5),
    size = 3,
    color = "white",
    fontface = "bold"
  )

# 保存图形
ggsave(
  filename = "P_site_distribution_colors_modified.png",
  plot = p,
  width = 10,
  height = 6,
  dpi = 300
)

# 同时保存PDF版本
ggsave(
  filename = "P_site_distribution_colors_modified.pdf",
  plot = p,
  width = 10,
  height = 6
)


################################################################################################################
# ==============================================================================
# 第九部分：aTIS周围读码框分布图
# ==============================================================================
cat("\n9. 生成aTIS周围读码框分布图...\n")
cat("----------------------------------------\n")

tryCatch({
  # 确保有reads_psite_list数据
  if (!exists("reads_psite_list") || length(reads_psite_list) == 0) {
    cat("✗ 缺少reads_psite_list数据，跳过此步骤\n")
  } else {
    # 分离控制组和实验组
    ev_samples <- names(reads_psite_list)[grepl("^EV", names(reads_psite_list))]
    orf3a_samples <- names(reads_psite_list)[grepl("ORF3a", names(reads_psite_list))]
    
    if (length(ev_samples) == 0 && length(orf3a_samples) == 0) {
      cat("✗ 未找到EV或ORF3a样本，跳过此步骤\n")
    } else {
      # 为每个样本添加处理条件标签
      reads_with_condition <- rbindlist(lapply(names(reads_psite_list), function(sample_name) {
        dt <- copy(reads_psite_list[[sample_name]])
        dt[, sample := sample_name]
        # 根据样本名称分配组别
        if (grepl("^EV", sample_name)) {
          dt[, condition := "EV"]
        } else if (grepl("ORF3a", sample_name)) {
          dt[, condition := "ORF3a"]
        } else {
          dt[, condition := "Other"]
        }
        return(dt)
      }), fill = TRUE)
      
      # 只保留EV和ORF3a样本
      reads_with_condition <- reads_with_condition[condition %in% c("EV", "ORF3a")]
      
      # 9.1 计算读码框分布（关注起始位点）
      cat("9.1 计算aTIS周围的读码框分布...\n")
      
      # 首先计算P-site的读码框
      reads_with_condition[, frame := psite_from_start %% 3]
      
      # 筛选起始位点周围的reads（距离起始位点-50到+150nt）
      tss_reads <- reads_with_condition[
        psite_from_start >= -50 & psite_from_start <= 150,
        .(sample, condition, position = psite_from_start, frame)
      ]
      
      # 计算每个位置每个frame的read密度
      frame_density <- tss_reads[
        , .(count = .N), by = .(condition, position, frame)
      ][
        , frequency := count / sum(count), by = .(condition, frame)
      ]
      
      # 标准化：使每个frame的总和为1
      frame_normalized <- frame_density[
        , .(position, frame, 
            normalized_density = frequency / max(frequency, na.rm = TRUE)), 
        by = .(condition)
      ]
      
      # 9.2 创建类似示例图的图形
      cat("9.2 创建aTIS周围读码框分布图...\n")
      
      library(ggplot2)
      
      # 创建淡色调颜色方案
      frame_colors <- c(
        "0" = "#D8BFD8",  # Frame 0 - 淡紫色 (Thistle)
        "1" = "#98FB98",  # Frame 1 - 淡绿色 (PaleGreen)
        "2" = "#ADD8E6"   # Frame 2 - 淡蓝色 (LightBlue)
      )
      
      # 创建分面图
      atis_frame_plot <- ggplot(frame_normalized, 
                                aes(x = position, y = normalized_density, 
                                    color = as.factor(frame))) +
        geom_line(linewidth = 1.2, alpha = 0.9) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
        
        # 分面显示不同条件
        facet_grid(condition ~ ., scales = "free_y") +
        
        # 颜色映射
        scale_color_manual(
          name = "Frame",
          values = frame_colors,
          labels = c("Frame 0", "Frame 1", "Frame 2")
        ) +
        
        # 坐标轴和标签
        scale_x_continuous(
          breaks = seq(-50, 150, 25),
          limits = c(-50, 150),
          expand = expansion(mult = 0.02)
        ) +
        
        labs(
          title = "Frame Distribution around Translation Start Sites",
          x = "Distance from start codon (nt)",
          y = "Normalized read density",
          color = "Reading Frame"
        ) +
        
        # 主题设置 - 简洁风格，无网格线
        theme_bw(base_size = 12) +
        theme(
          plot.title = element_text(
            size = 16, 
            face = "bold", 
            hjust = 0.5,
            margin = margin(b = 15)
          ),
          axis.title = element_text(size = 13, face = "bold"),
          axis.text = element_text(size = 11, color = "black"),
          axis.text.x = element_text(margin = margin(t = 5)),
          axis.text.y = element_text(margin = margin(r = 5)),
          
          # 坐标轴实线
          axis.line = element_line(color = "black", linewidth = 0.5),
          axis.ticks = element_line(color = "black", linewidth = 0.5),
          
          # 移除网格线
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          
          # 分面标签
          strip.text = element_text(
            size = 13, 
            face = "bold",
            margin = margin(t = 10, b = 10)
          ),
          strip.background = element_rect(
            fill = "gray95",
            color = "gray70",
            linewidth = 0.8
          ),
          
          # 图例设置
          legend.title = element_text(size = 12, face = "bold"),
          legend.text = element_text(size = 11),
          legend.position = "top",
          legend.box = "horizontal",
          legend.margin = margin(b = 10),
          legend.key = element_rect(fill = "white", color = NA),
          
          # 面板设置
          panel.border = element_rect(color = "gray70", linewidth = 0.5, fill = NA),
          panel.spacing = unit(15, "points"),
          
          # 边距
          plot.margin = margin(20, 25, 20, 20),
          plot.background = element_rect(fill = "white", color = NA)
        )
      
      # 9.3 保存图形
      cat("9.3 保存图形...\n")
      
      # 保存高质量图形
      ggsave(
        filename = file.path(paths$output_dir, "plots", "frame_distribution_start_site.pdf"),
        plot = atis_frame_plot,
        width = 10,
        height = 8,
        dpi = 300
      )
      
      ggsave(
        filename = file.path(paths$output_dir, "plots", "frame_distribution_start_site.png"),
        plot = atis_frame_plot,
        width = 10,
        height = 8,
        dpi = 300,
        bg = "white"
      )
      
      # 9.4 保存数据
      fwrite(
        frame_normalized,
        file.path(paths$output_dir, "tables", "frame_distribution_start_site.csv")
      )
      
      # 9.5 显示数据摘要
      cat("\n   起始位点周围读码框分布摘要:\n")
      for (cond in unique(frame_normalized$condition)) {
        cat(paste0("\n   ", cond, "组:\n"))
        cond_data <- frame_normalized[condition == cond]
        for (fr in 0:2) {
          fr_data <- cond_data[frame == fr]
          if (nrow(fr_data) > 0) {
            peak_pos <- fr_data[which.max(normalized_density), position]
            cat(paste0("     Frame ", fr, ": 峰值位置 ", peak_pos, "nt\n"))
          }
        }
      }
      
      cat("✓ 起始位点周围读码框分布图已保存\n")
      
      # 9.6 可选：生成简化的统计摘要图（无网格线版本）
      cat("9.6 生成读码框偏好性摘要图...\n")
      
      # 计算每个frame的总密度
      frame_summary <- frame_density[
        , .(total_density = sum(count)), by = .(condition, frame)
      ][
        , percentage := total_density / sum(total_density) * 100, by = condition
      ]
      
      # 创建摘要柱状图 - 简洁无网格版本
      summary_plot <- ggplot(frame_summary, 
                            aes(x = as.factor(frame), y = percentage, 
                                fill = as.factor(frame))) +
        geom_bar(stat = "identity", width = 0.7, alpha = 0.9, 
                color = "gray40", linewidth = 0.3) +
        facet_grid(. ~ condition) +
        
        # 使用相同的淡色调
        scale_fill_manual(values = frame_colors, guide = "none") +
        
        # 添加数值标签
        geom_text(aes(label = sprintf("%.1f%%", percentage)),
                  vjust = -0.5, size = 3.5, fontface = "bold") +
        
        labs(
          title = "Frame Preference at Start Site Region",
          x = "Reading Frame",
          y = "Percentage (%)"
        ) +
        
        # 简洁主题，无网格线
        theme_bw(base_size = 12) +
        theme(
          plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          axis.title = element_text(face = "bold", size = 12),
          axis.text = element_text(size = 11, color = "black"),
          
          # 坐标轴实线
          axis.line = element_line(color = "black", linewidth = 0.5),
          axis.ticks = element_line(color = "black", linewidth = 0.5),
          
          # 移除网格线
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          
          # 分面设置
          strip.text = element_text(face = "bold", size = 12),
          strip.background = element_rect(fill = "gray95", color = "gray70"),
          
          # 面板边框
          panel.border = element_rect(color = "gray70", linewidth = 0.8),
          
          plot.margin = margin(15, 15, 15, 15)
        ) +
        
        # 设置y轴范围，为标签留空间
        ylim(0, max(frame_summary$percentage) * 1.1)
      
      # 保存摘要图
      ggsave(
        filename = file.path(paths$output_dir, "plots", "frame_preference_summary.pdf"),
        plot = summary_plot,
        width = 8,
        height = 6,
        dpi = 300
      )
      
      ggsave(
        filename = file.path(paths$output_dir, "plots", "frame_preference_summary.png"),
        plot = summary_plot,
        width = 8,
        height = 6,
        dpi = 300,
        bg = "white"
      )
      
      # 保存摘要数据
      fwrite(
        frame_summary,
        file.path(paths$output_dir, "tables", "frame_preference_summary.csv")
      )
      
      cat("✓ 读码框偏好性摘要图已保存\n")
    }
  }
}, error = function(e) {
  cat("✗ 起始位点读码框分布图生成失败: ", e$message, "\n")
})

# ==============================================================================
# 第十部分：完成
# ==============================================================================
cat("\n========================================\n")
cat("分析完成！\n")
cat("结果保存在:", paths$output_dir, "\n")
cat("========================================\n")

##############################################################CDS全长############################################################################################

# ==============================================================================
# 第九部分-D：CDS全长读码框分布折线图
# ==============================================================================
cat("\n9-D. 生成CDS全长读码框分布折线图...\n")
cat("----------------------------------------\n")

tryCatch({
  if (exists("reads_psite_list") && length(reads_psite_list) > 0) {
    # 创建样本分组
    sample_groups <- list()
    ev_samples <- names(reads_psite_list)[grepl("^EV", names(reads_psite_list))]
    orf3a_samples <- names(reads_psite_list)[grepl("ORF3a", names(reads_psite_list))]
    
    if (length(ev_samples) > 0) sample_groups[["EV"]] <- ev_samples
    if (length(orf3a_samples) > 0) sample_groups[["ORF3a"]] <- orf3a_samples
    
    if (length(sample_groups) == 0) {
      cat("✗ 未找到EV或ORF3a样本\n")
    } else {
      cat("9-D.1 计算CDS全长的frame特异性metaprofile...\n")
      
      # 我们需要分别计算每个frame的分布
      # 首先提取所有CDS区域的reads
      cds_reads_list <- lapply(reads_psite_list, function(sample_data) {
        # 只保留CDS区域
        cds_data <- sample_data[psite_region == "cds", ]
        
        # 计算frame
        cds_data[, frame := psite_from_start %% 3]
        
        # 计算相对位置（标准化到0-100）
        cds_data[, relative_pos := (psite_from_start / 
                                     (psite_from_start + abs(psite_from_stop))) * 100]
        
        return(cds_data)
      })
      
      # 合并所有样本并按组处理
      all_frame_data <- rbindlist(lapply(names(cds_reads_list), function(sample_name) {
        sample_data <- cds_reads_list[[sample_name]]
        if (nrow(sample_data) > 0) {
          # 确定分组
          group_name <- ifelse(grepl("^EV", sample_name), "EV", "ORF3a")
          
          # 创建bins（将CDS分成50个等分）
          sample_data[, bin := cut(relative_pos, breaks = 50, labels = FALSE)]
          
          # 计算每个bin中每个frame的计数
          frame_counts <- sample_data[
            , .(count = .N), by = .(bin, frame)
          ][
            , proportion := count / sum(count), by = .(bin)
          ]
          
          # 添加样本和组信息
          frame_counts[, sample := sample_name]
          frame_counts[, group := group_name]
          
          return(frame_counts)
        }
      }), fill = TRUE)
      
      # 按组平均
      frame_summary <- all_frame_data[
        , .(mean_proportion = mean(proportion, na.rm = TRUE),
            se = sd(proportion, na.rm = TRUE) / sqrt(.N)),
        by = .(group, bin, frame)
      ]
      
      # 转换为实际位置（bin 1-50 对应 2%-98%）
      frame_summary[, position := (bin - 0.5) * 2]  # 每个bin代表2%的CDS长度
      
      cat("9-D.2 创建CDS全长读码框分布折线图...\n")
      
      # 创建淡色调颜色方案
      frame_colors <- c(
        "0" = "#D8BFD8",  # Frame 0 - 淡紫色
        "1" = "#98FB98",  # Frame 1 - 淡绿色
        "2" = "#ADD8E6"   # Frame 2 - 淡蓝色
      )
      
      # 创建折线图（分面显示不同组）
      cds_frame_lineplot <- ggplot(frame_summary, 
                                   aes(x = position, y = mean_proportion * 100, 
                                       color = as.factor(frame),
                                       group = as.factor(frame))) +
        geom_line(linewidth = 1.2, alpha = 0.9) +
        
        # 添加误差带（标准误）
        geom_ribbon(aes(ymin = (mean_proportion - se) * 100, 
                       ymax = (mean_proportion + se) * 100, 
                       fill = as.factor(frame)),
                   alpha = 0.15, color = NA) +
        
        # 分面显示
        facet_grid(group ~ ., scales = "free_y") +
        
        # 颜色和填充映射
        scale_color_manual(
          name = "Reading Frame",
          values = frame_colors,
          labels = c("Frame 0", "Frame 1", "Frame 2")
        ) +
        
        scale_fill_manual(
          name = "Reading Frame",
          values = frame_colors,
          labels = c("Frame 0", "Frame 1", "Frame 2"),
          guide = "none"  # 隐藏填充图例
        ) +
        
        # 坐标轴设置
        scale_x_continuous(
          breaks = seq(0, 100, 25),
          labels = c("Start\n(0%)", "25%", "50%", "75%", "Stop\n(100%)"),
          limits = c(0, 100),
          expand = expansion(mult = 0.02)
        ) +
        
        scale_y_continuous(
          expand = expansion(mult = 0.05)
        ) +
        
        labs(
          title = "Frame Distribution Across Entire CDS",
          subtitle = "Position-dependent frame preference in coding sequence",
          x = "Relative Position in CDS",
          y = "Frame Proportion (%)",
          color = "Reading Frame"
        ) +
        
        # 简洁主题
        theme_bw(base_size = 12) +
        theme(
          plot.title = element_text(
            size = 16, 
            face = "bold", 
            hjust = 0.5,
            margin = margin(b = 10)
          ),
          plot.subtitle = element_text(
            size = 12,
            hjust = 0.5,
            color = "gray40",
            margin = margin(b = 15)
          ),
          axis.title = element_text(size = 13, face = "bold"),
          axis.text = element_text(size = 11, color = "black"),
          axis.text.x = element_text(margin = margin(t = 5)),
          
          # 坐标轴实线
          axis.line = element_line(color = "black", linewidth = 0.5),
          axis.ticks = element_line(color = "black", linewidth = 0.5),
          
          # 移除网格线
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          
          # 分面标签
          strip.text = element_text(
            size = 13, 
            face = "bold",
            margin = margin(t = 10, b = 10)
          ),
          strip.background = element_rect(
            fill = "gray95",
            color = "gray70",
            linewidth = 0.8
          ),
          
          # 图例
          legend.title = element_text(size = 12, face = "bold"),
          legend.text = element_text(size = 11),
          legend.position = "top",
          legend.box = "horizontal",
          legend.margin = margin(b = 10),
          
          # 面板边框
          panel.border = element_rect(color = "gray70", linewidth = 0.8, fill = NA),
          
          # 边距
          plot.margin = margin(20, 25, 20, 20)
        )
      
      # 9-D.3 保存图形
      cat("9-D.3 保存CDS全长读码框分布折线图...\n")
      
      ggsave(
        filename = file.path(paths$output_dir, "plots", "frame_distribution_cds_lineplot.pdf"),
        plot = cds_frame_lineplot,
        width = 10,
        height = 8,
        dpi = 300
      )
      
      ggsave(
        filename = file.path(paths$output_dir, "plots", "frame_distribution_cds_lineplot.png"),
        plot = cds_frame_lineplot,
        width = 10,
        height = 8,
        dpi = 300,
        bg = "white"
      )
      
      # 保存数据
      fwrite(
        frame_summary,
        file.path(paths$output_dir, "tables", "frame_distribution_cds_lineplot.csv")
      )
      
      cat("✓ CDS全长读码框分布折线图已保存\n")
      
      # 显示关键统计
      cat("\n   CDS全长frame分布关键统计:\n")
      for (grp in unique(frame_summary$group)) {
        cat(paste0("\n   ", grp, "组:\n"))
        for (fr in 0:2) {
          fr_data <- frame_summary[group == grp & frame == fr]
          avg_prop <- mean(fr_data$mean_proportion) * 100
          cat(paste0("     Frame ", fr, ": 平均占比 ", round(avg_prop, 1), "%\n"))
        }
      }
    }
  }
}, error = function(e) {
  cat("✗ CDS全长读码框分布折线图生成失败: ", e$message, "\n")
})