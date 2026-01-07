rm(list=ls())
library(ggplot2)
library(ggrepel)
setwd("Z:\\16.up_down_go_reactome\\01.WT_IFNR\\heatmap_volcano")
data <- read.csv("DEanalysis.csv",header = TRUE,row.names = 1)
# 计算-log10(pvalue)
data$pvalue <- -log10(data$pvalue)
# 分别选出上调和下调前15个pvalue最小的基因
up.lst <- data[!is.na(data$log2FoldChange) & !is.na(data$padj) & data$log2FoldChange > 1 & data$padj < 0.05 & data$baseMean >100,]
dn.lst <- data[!is.na(data$log2FoldChange) & !is.na(data$padj) & data$log2FoldChange < -1 & data$padj < 0.05 & data$baseMean >100,]
top_upregulated_by_pvalue <- head(up.lst[order(-up.lst$log2FoldChange), ], 15) # 上调基因按log2FC排序，然后按pvalue排序
top_downregulated_by_pvalue <- head(dn.lst[order(dn.lst$log2FoldChange), ], 15) # 下调基因按log2FC排序，然后按pvalue排序
# 合并两个数据框，用于在图中标注
to_label <- rbind(top_upregulated_by_pvalue, top_downregulated_by_pvalue)
# 标记需要标注的基因点
data$highlight <- ifelse(rownames(data) %in% rownames(to_label), "Highlighted", "Normal")

pdf(file="volcano_1.pdf",width=10,height = 12)
# 绘制火山图
ggplot(data, aes(x = log2FoldChange, y = pvalue)) +
  # 上调基因用红色点表示
  geom_point(data = subset(data, log2FoldChange > 1 & highlight == "Normal"), aes(color = "Upregulated"), size = 2.5) +
  # 下调基因用蓝色点表示
  geom_point(data = subset(data, log2FoldChange < -1 & highlight == "Normal"), aes(color = "Downregulated"), size = 2.5) +
  # log2FC在-1到1之间的基因用灰色点表示
  geom_point(data = subset(data, log2FoldChange >= -1 & log2FoldChange <= 1), aes(color = "Non-significant"), size = 2.5) +
  # 高亮前15个基因点为橙色
  geom_point(data = subset(data, highlight == "Highlighted"), aes(color = "Highlighted"), size = 2.5) +
  # 标注上调和下调的前15个基因，使用橙色
  geom_text_repel(
    data = top_downregulated_by_pvalue,
    aes(label = rownames(top_downregulated_by_pvalue)),
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
dev.off()


