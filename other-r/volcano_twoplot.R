rm(list=ls())
library(ggplot2)
library(ggrepel)
setwd("G:\\temp\\23.volcano_data")
data1 <- read.table("volcano_1.txt")
# 计算-log10(pvalue)
data1$pvalue <- -log10(data1$pvalue)
# 分别选出上调和下调前15个pvalue最小的基因
top_downregulated_by_pvalue1 <- head(data1[order(data1$log2FoldChange), ], 15) # 下调基因按log2FC排序，然后按pvalue排序
top_upregulated_by_pvalue1 <- head(data1[order(-data1$log2FoldChange), ], 15) # 上调基因按log2FC排序，然后按pvalue排序
# 合并两个数据框，用于在图中标注
to_label1 <- rbind(top_upregulated_by_pvalue1, top_downregulated_by_pvalue1)
# 标记需要标注的基因点
data1$highlight <- ifelse(rownames(data1) %in% rownames(to_label1), "Highlighted", "Normal")

data2 <- read.csv("G:\\temp\\valcano\\my_newdata.csv",header = TRUE,sep = ",", quote = "\"",dec = ".",row.names = 1)
data2$pvalue <- -log10(data2$pvalue)
top_upregulated_by_pvalue2 <- head(data2[order(data2$log2FoldChange), ], 15) # 上调基因按log2FC排序，然后按pvalue排序
top_downregulated_by_pvalue2 <- head(data2[order(-data2$log2FoldChange), ], 15) # 下调基因按log2FC排序，然后按pvalue排序
# 合并两个数据框，用于在图中标注
to_label2 <- rbind(top_upregulated_by_pvalue2, top_downregulated_by_pvalue2)
# 标记需要标注的基因点
data2$highlight <- ifelse(rownames(data2) %in% rownames(to_label2), "Highlighted", "Normal")

library(patchwork)
pdf(file="volcano_all.pdf",width=20,height = 11)
# 绘制火山图
p1 <- ggplot(data2, aes(x = log2FoldChange, y = pvalue)) +
  # 上调基因用红色点表示
  geom_point(data = subset(data2, log2FoldChange > 1 & highlight == "Normal"), aes(color = "Upregulated"), size = 2.5) +
  # 下调基因用蓝色点表示
  geom_point(data = subset(data2, log2FoldChange < -1 & highlight == "Normal"), aes(color = "Downregulated"), size = 2.5) +
  # log2FC在-1到1之间的基因用灰色点表示
  geom_point(data = subset(data2, log2FoldChange >= -1 & log2FoldChange <= 1), aes(color = "Non-significant"), size = 2.5) +
  # 高亮前15个基因点为橙色
  geom_point(data = subset(data2, highlight == "Highlighted"), aes(color = "Highlighted"), size = 2.5) +
  # 标注上调和下调的前15个基因，使用橙色
  geom_text_repel(
    data = to_label2,
    aes(label = rownames(to_label2)),
    color = "black",
    size = 10,
    show.legend = FALSE,
    max.overlaps = Inf,
    fontface="bold"
  ) +
  xlim(-11,11)+
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Highlighted" = "orange", "Non-significant" = "grey")) + # 设置颜色
  theme_minimal() +
  labs(title = "Top genes based on log2Fold change",
       x = "log2FC",
       y = "-log10(pvalue)",
       color = "Regulation") +
  theme_classic() +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey", size = 0.5) +
  geom_vline(xintercept = -1, linetype = "dashed", color = "grey", size = 0.5) +
  # 调整主题以使标题居中并增大字体
  theme(
    plot.title = element_text(hjust = 0.5, size = 30,face = "bold"),  # 居中并增大字体
    axis.text = element_text(size = 30,face = "bold"),
    axis.title = element_text(size=30,face = "bold"),
    legend.position = "none",
    axis.line = element_line(color = "black",size = 2.5),
    axis.ticks.length=unit(0.4, "cm"),
    axis.ticks=element_line(size = 2.5)
  ) 
  
p2 <- ggplot(data1, aes(x = log2FoldChange, y = pvalue)) +
  # 上调基因用红色点表示
  geom_point(data = subset(data1, log2FoldChange > 1 & highlight == "Normal"), aes(color = "Upregulated"), size = 2.5) +
  # 下调基因用蓝色点表示
  geom_point(data = subset(data1, log2FoldChange < -1 & highlight == "Normal"), aes(color = "Downregulated"), size = 2.5) +
  # log2FC在-1到1之间的基因用灰色点表示
  geom_point(data = subset(data1, log2FoldChange >= -1 & log2FoldChange <= 1), aes(color = "Non-significant"), size = 2.5) +
  # 高亮前15个基因点为橙色
  geom_point(data = subset(data1, highlight == "Highlighted"), aes(color = "Highlighted"), size = 2.5) +
  # 标注上调和下调的前15个基因，使用橙色
  geom_text_repel(
    data = top_downregulated_by_pvalue1,
    aes(label = rownames(top_downregulated_by_pvalue1)),
    color = "black",
    size = 10,
    show.legend = FALSE,
    max.overlaps = Inf,
    fontface="bold"
  ) +
  xlim(-11,5)+
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Highlighted" = "orange", "Non-significant" = "grey")) + # 设置颜色
  theme_minimal() +
  labs(title = "Top genes based on log2Fold change",
       x = "log2FC",
       y = "-log10(pvalue)",
       color = "Regulation") +
  theme_classic() +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey", size = 0.5) +
  geom_vline(xintercept = -1, linetype = "dashed", color = "grey", size = 0.5) +
  # 调整主题以使标题居中并增大字体
  theme(
    plot.title = element_text(hjust = 0.5, size = 30,face = "bold"),  # 居中并增大字体
    axis.text = element_text(size = 30,face = "bold"),
    axis.title = element_text(size=30,face = "bold"),
    legend.position = "none",
    axis.line = element_line(color = "black",size = 2.5),
    axis.ticks.length=unit(0.4, "cm"),
    axis.ticks=element_line(size = 2.5)
  ) 

p1+p2

dev.off()






