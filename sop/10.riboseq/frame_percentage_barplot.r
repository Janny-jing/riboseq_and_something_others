library(ggplot2)
library(dplyr)

# 读取所有TSV文件
files <- c("EV-1_plot_data.tsv", "EV-2_plot_data.tsv", "EV-3_plot_data.tsv",
           "orf3a-1_plot_data.tsv", "orf3a-2_plot_data.tsv", "orf3a-3_plot_data.tsv")

# 检查文件存在性
existing_files <- files[file.exists(files)]
cat("找到以下文件:", paste(existing_files, collapse = ", "), "\n\n")

# 初始化数据框
all_data <- data.frame()

for (file in existing_files) {
  # 读取数据
  data <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # 从文件名提取样本名
  sample_name <- gsub("_plot_data\\.tsv", "", file)
  
  # 标准化样本名
  if (grepl("^orf3a", sample_name, ignore.case = TRUE)) {
    sample_name <- gsub("^orf3a", "ORF3a", sample_name, ignore.case = TRUE)
  } else if (grepl("^ev", sample_name, ignore.case = TRUE)) {
    sample_name <- toupper(sample_name)
  }
  
  cat("处理文件:", file, "-> 样本:", sample_name, "\n")
  cat("  数据格式: Frame =", paste(data$Frame, collapse = ","), 
      "Percentage =", paste(data$Percentage, collapse = ","), "\n")
  
  # 添加样本名列
  data$Sample <- sample_name
  
  # 合并数据
  all_data <- rbind(all_data, data)
}

# 查看数据结构
cat("\n合并后的数据结构:\n")
print(str(all_data))

# 重命名列（如果需要）
colnames(all_data) <- c("Frame", "Count", "Percentage", "Sample")

# 确保百分比是数值（如果有百分号，去掉）
if (is.character(all_data$Percentage)) {
  all_data$Percentage <- as.numeric(gsub("%", "", all_data$Percentage))
}

# 将Frame转换为因子并添加"Frame"前缀
all_data$Frame <- factor(all_data$Frame, 
                         levels = c(0, 1, 2),
                         labels = c("Frame0", "Frame1", "Frame2"))


# 设置样本顺序
sample_order <- c("EV-1", "EV-2", "EV-3", "ORF3a-1", "ORF3a-2", "ORF3a-3")
all_data$Sample <- factor(all_data$Sample, levels = sample_order)

# 保存合并数据
write.csv(all_data, "combined_frame_data.csv", row.names = FALSE)
cat("\n合并数据已保存到: combined_frame_data.csv\n")

# 创建摘要统计
summary_data <- all_data %>%
  group_by(Sample, Frame) %>%
  summarise(Percentage = mean(Percentage), .groups = "drop")

cat("\n数据摘要:\n")
print(summary_data)

# 方法1：堆叠柱状图（百分比）
p1 <- ggplot(all_data, aes(x = Sample, y = Percentage, fill = Frame)) +
  geom_bar(stat = "identity", position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values = c("Frame0" = "#E41A1C", 
                               "Frame1" = "#377EB8", 
                               "Frame2" = "#4DAF4A")) +
  labs(title = "Ribo-seq Frame Distribution (Counts)",
       y = "Read Count Percentage", 
       x = "Sample") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    axis.text.y = element_text(size = 10),
    legend.position = "right"
  )

# 方法2：填充柱状图（标准化到100%）
p2 <- ggplot(all_data, aes(x = Sample, y = Percentage, fill = Frame)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = c("Frame0" = "#E41A1C", 
                               "Frame1" = "#377EB8", 
                               "Frame2" = "#4DAF4A")) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "Ribo-seq Frame Distribution (Normalized)",
       y = "Proportion", 
       x = "Sample") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    axis.text.y = element_text(size = 10),
    legend.position = "right"
  )

# 显示图形
print(p1)
print(p2)

# 保存图形
ggsave("frame_distribution_counts.pdf", p1, width = 10, height = 6, bg = "white")
ggsave("frame_distribution_counts.png", p1, width = 10, height = 6, dpi = 300, bg = "white")

ggsave("frame_distribution_normalized.pdf", p2, width = 10, height = 6, bg = "white")
ggsave("frame_distribution_normalized.png", p2, width = 10, height = 6, dpi = 300, bg = "white")
