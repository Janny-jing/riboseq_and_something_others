setwd("Y:\\temp\\20251122-riboseq_new_3sample\\single_bam")
rm(list=ls())
library(riboSeqR)
library(rtracklayer)
library(GenomicAlignments)
library(ggplot2)
library(dplyr)
library(ChIPseeker)
library(GenomicFeatures)
# 自动获取bam文件
bam_files <- list.files(pattern = "\\.bam$")
sample_names <- gsub("\\.bam$", "", bam_files)

# 自动识别分组
groups <- ifelse(grepl("orf3a", sample_names), "orf3a", "EV")

cat("找到", length(bam_files), "个bam文件:\n")
print(data.frame(Sample = sample_names, Group = groups))

# 设置GTF文件路径
gtf_file <- "Y:/temp/20251122-riboseq_new_3sample/gencode.v48.annotation.gtf"

## 第一步：从GTF文件提取CDS区域
cat("\n正在从GTF文件提取CDS区域...\n")
gtf <- rtracklayer::import(gtf_file)
cds_regions <- gtf[gtf$type == "CDS"]
cat("成功提取", length(cds_regions), "个CDS区域\n")

## 第二步：使用riboSeqR进行专业frame分析
cat("\n开始Ribo-seq Frame分析...\n")
frame_results <- list()

for (i in 1:length(bam_files)) {
  cat("正在分析:", bam_files[i], "...\n")
  
  # 确保bam文件有索引
  if (!file.exists(paste0(bam_files[i], ".bai"))) {
    cat("为", bam_files[i], "创建索引...\n")
    indexBam(bam_files[i])
  }
  
  # 使用riboSeqR读取ribo数据
  riboDat <- readRibodata(bam_files[i], replicates = groups[i])
  
  # 进行frame counting
  fCs <- frameCounting(riboDat, cds_regions)
  
  # 获取reading frame结果
  fS <- readingFrame(rC = fCs)
  
  frame_results[[sample_names[i]]] <- list(
    frame_counts = fCs,
    reading_frame = fS
  )
  cat("✓ 完成", bam_files[i], "的frame分析\n")
}

## 第三步：正确提取riboSeqR的frame结果并可视化
cat("\n处理frame分析结果...\n")
frame_plot_data <- data.frame()

for (sample in names(frame_results)) {
  fS <- frame_results[[sample]]$reading_frame
  
  # 直接使用reading_frame矩阵，前3行就是frame 0,1,2的counts
  # fS是一个矩阵，行1-3对应frame 0-2，列对应read length
  frame_matrix <- fS[1:3, ]  # 提取前3行（frame 0,1,2）
  
  # 汇总所有read length的frame counts
  total_frame0 <- sum(frame_matrix[1, ], na.rm = TRUE)  # frame 0
  total_frame1 <- sum(frame_matrix[2, ], na.rm = TRUE)  # frame 1  
  total_frame2 <- sum(frame_matrix[3, ], na.rm = TRUE)  # frame 2
  
  total_reads <- total_frame0 + total_frame1 + total_frame2
  
  # 计算百分比
  frame0_pct <- (total_frame0 / total_reads) * 100
  frame1_pct <- (total_frame1 / total_reads) * 100
  frame2_pct <- (total_frame2 / total_reads) * 100
  
  # 创建数据框
  temp_df <- data.frame(
    Sample = sample,
    Group = groups[which(sample_names == sample)],
    Frame = c("Frame0", "Frame1", "Frame2"),
    Count = c(total_frame0, total_frame1, total_frame2),
    Percentage = c(frame0_pct, frame1_pct, frame2_pct),
    stringsAsFactors = FALSE
  )
  
  frame_plot_data <- rbind(frame_plot_data, temp_df)
}

# 查看前几行数据确认
cat("Frame数据预览:\n")
print(head(frame_plot_data))

# 修改样本显示名称（与之前一致）
frame_plot_data$Sample <- ifelse(grepl("EV-1", frame_plot_data$Sample), "EV-1",
                                 ifelse(grepl("EV-2", frame_plot_data$Sample), "EV-2",
                                        ifelse(grepl("EV-3", frame_plot_data$Sample), "EV-3",
                                               ifelse(grepl("orf3a-1", frame_plot_data$Sample), "ORF3a-1",
                                                      ifelse(grepl("orf3a-2", frame_plot_data$Sample), "ORF3a-2",
                                                             ifelse(grepl("orf3a-3", frame_plot_data$Sample), "ORF3a-3", frame_plot_data$Sample))))))

# 设置正确的顺序
sample_order <- c("EV-1", "EV-2", "EV-3", "ORF3a-1", "ORF3a-2", "ORF3a-3")
frame_plot_data$Sample <- factor(frame_plot_data$Sample, levels = sample_order)

# 绘制frame分布图
p_frame <- ggplot(frame_plot_data, aes(x = Sample, y = Percentage, fill = Frame)) +
  geom_bar(stat = "identity", position = "fill", width = 0.8) +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  scale_fill_manual(values = c("Frame0" = "#E41A1C", "Frame1" = "#377EB8", "Frame2" = "#4DAF4A")) +
  labs(title = "Ribo-seq Reading Frame Distribution",
       subtitle = "Professional analysis using riboSeqR package",
       y = "Percentage", x = "Sample") +
  theme_classic() +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    axis.line.x = element_blank(),
    axis.line.y = element_line(color = "grey50", linewidth = 0.3),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.2),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    legend.position = "right",
    legend.background = element_rect(fill = "white", colour = NA)
  )

if (length(unique(frame_plot_data$Group)) > 1) {
  p_frame <- p_frame + facet_grid(~ Group, scales = "free_x", space = "free")
}


# 保存frame分布结果
write.csv(frame_plot_data, "professional_riboseq_frame_results.csv", row.names = FALSE)
cat("Frame分布结果已保存到: professional_riboseq_frame_results.csv\n")

# 保存图表
ggsave("frame_distribution.pdf", p_frame, width = 10, height = 6, bg = "white")
ggsave("frame_distribution.png", p_frame, width = 10, height = 6, dpi = 300, bg = "white")

## 第四步：生成正确的统计摘要
cat("\nFrame分布统计摘要:\n")
frame_summary <- frame_plot_data %>%
  group_by(Frame, Group) %>%
  summarise(
    Mean_Percentage = mean(Percentage, na.rm = TRUE),
    SD = sd(Percentage, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  arrange(Group, Frame)

print(frame_summary)

# 数据质量评估
frame0_data <- frame_plot_data[frame_plot_data$Frame == "Frame0", ]
cat("\n数据质量评估:\n")
cat("Frame 0 比例范围:", round(min(frame0_data$Percentage), 1), "-", 
    round(max(frame0_data$Percentage), 1), "%\n")
cat("Frame 0 平均值:", round(mean(frame0_data$Percentage), 1), "%\n")

if (mean(frame0_data$Percentage) > 50) {
  cat("✓ 优秀: 强烈的triplet periodicity信号\n")
} else if (mean(frame0_data$Percentage) > 40) {
  cat("○ 良好: 明显的triplet periodicity信号\n")
} else if (mean(frame0_data$Percentage) > 33) {
  cat("△ 一般: 有一定的triplet periodicity信号\n")
} else {
  cat("⚠ 注意: triplet periodicity信号较弱\n")
}

# 查看每个样本的具体frame分布
cat("\n各样本Frame分布详情:\n")
sample_frame_details <- frame_plot_data %>%
  group_by(Sample, Group) %>%
  summarise(
    Frame0_Pct = Percentage[Frame == "Frame0"],
    Frame1_Pct = Percentage[Frame == "Frame1"], 
    Frame2_Pct = Percentage[Frame == "Frame2"],
    .groups = 'drop'
  )
print(sample_frame_details)







## 第四步：专业区域注释分析
cat("\n开始专业区域注释分析...\n")

# 创建TxDb对象用于ChIPseeker
txdb <- makeTxDbFromGRanges(gtf)

annotation_results <- list()

for (i in 1:length(bam_files)) {
  cat("正在注释:", bam_files[i], "...\n")
  
  # 读取比对结果
  reads <- readGAlignments(bam_files[i])
  reads_gr <- granges(reads)
  
  # 使用ChIPseeker进行专业注释
  peak_anno <- annotatePeak(reads_gr, 
                            tssRegion = c(0, 0),
                            TxDb = txdb,
                            level = "transcript",
                            assignGenomicAnnotation = TRUE,
                            genomicAnnotationPriority = c("5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic"))
  
  annotation_results[[sample_names[i]]] <- peak_anno
  cat("✓ 完成", bam_files[i], "的专业注释\n")
}

# 处理注释结果并重新分类
annotation_plot_data <- data.frame()

for (sample in names(annotation_results)) {
  anno_df <- as.data.frame(annotation_results[[sample]])
  anno_summary <- table(anno_df$annotation)
  total_reads <- sum(anno_summary)
  
  # 创建新的分类汇总
  new_anno_summary <- c(
    CDS = sum(anno_summary[grepl("Exon", names(anno_summary))]),  # 合并所有Exon为CDS
    `5UTR` = ifelse("5' UTR" %in% names(anno_summary), anno_summary[["5' UTR"]], 0),
    `3UTR` = ifelse("3' UTR" %in% names(anno_summary), anno_summary[["3' UTR"]], 0),
    Intron = sum(anno_summary[grepl("Intron", names(anno_summary))]),  # 合并所有Intron
    Downstream = ifelse("Downstream" %in% names(anno_summary), anno_summary[["Downstream"]], 0),
    Intergenic = ifelse("Intergenic" %in% names(anno_summary), anno_summary[["Intergenic"]], 0)
  )
  
  # 移除值为0的类别
  new_anno_summary <- new_anno_summary[new_anno_summary > 0]
  
  for (anno_type in names(new_anno_summary)) {
    temp_df <- data.frame(
      Sample = sample,
      Group = groups[which(sample_names == sample)],
      Annotation = anno_type,
      Count = as.numeric(new_anno_summary[anno_type]),
      Percentage = (as.numeric(new_anno_summary[anno_type]) / total_reads) * 100,
      stringsAsFactors = FALSE
    )
    annotation_plot_data <- rbind(annotation_plot_data, temp_df)
  }
}

# 设置Annotation因子的顺序，让CDS在最前面
annotation_order <- c("CDS", "5UTR", "3UTR", "Intron", "Downstream", "Intergenic")
annotation_plot_data$Annotation <- factor(annotation_plot_data$Annotation, levels = annotation_order)

# 使用scale_x_discrete直接修改x轴标签
p_anno <- ggplot(annotation_plot_data, aes(x = Sample, y = Percentage, fill = Annotation)) +
  geom_bar(stat = "identity", position = "fill", width = 0.8) +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  scale_x_discrete(labels = c("EV-1", "EV-2", "EV-3", "ORF3a-1", "ORF3a-2", "ORF3a-3")) +
  scale_fill_brewer(palette = "Set1") +
  labs(title = "Genomic Region Annotation Distribution",
       subtitle = "Exon regions merged as CDS, Intron regions merged",
       y = "Percentage", x = "Sample") +
  theme_classic() +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    axis.line.x = element_blank(),
    axis.line.y = element_line(color = "grey50", linewidth = 0.3),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.2),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    legend.position = "right",
    legend.background = element_rect(fill = "white", colour = NA)
  )

if (length(unique(annotation_plot_data$Group)) > 1) {
  p_anno <- p_anno + facet_grid(~ Group, scales = "free_x", space = "free")
}




# 保存为PDF和PNG文件
ggsave("annotation_distribution.pdf", p_anno, width = 10, height = 6)
ggsave("annotation_distribution.png", p_anno, width = 10, height = 6, dpi = 300)

# 保存注释结果
write.csv(annotation_plot_data, "professional_annotation_results.csv", row.names = FALSE)
cat("区域注释结果已保存到: professional_annotation_results.csv\n")
cat("图表已保存为: annotation_distribution.pdf 和 annotation_distribution.png\n")


