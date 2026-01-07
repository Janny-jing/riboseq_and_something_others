rm(list = ls())
# 加载必要的库
library(pheatmap)
library(ggplot2)
library(openxlsx)
setwd("G:\\temp\\16.riboseq\\genes\\04.riboseq_head_tail_10_codons\\05\\02.isg_gene")
# 假设 heatmap_data 是一个数据框，color_map 是颜色映射名称（如 "RdYlBu"），font_size 是字体大小
# 并且 heatmap_file 是保存热图的文件路径

# 示例数据（请替换为您的实际数据）
heatmap_data1 <- read.table("ctl_isg_merge.txt",header = TRUE)
heatmap_data2 <- read.table("pic_isg_merge.txt",header = TRUE)
heatmap_data3 <- heatmap_data2-heatmap_data1
heatmap_file <- "G:\\temp\\16.riboseq\\genes\\04.riboseq_head_tail_10_codons\\05\\02.isg_gene\\pic-ctl_heatmap.png"
write.xlsx(heatmap_data3,file = "pic-ctl_combined_data.xlsx")
custom_color <- colorRampPalette(c("blue", "white", "red"))(100)
font_size <- 20
paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
myBreaks <- c(seq(min(heatmap_data3), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(heatmap_data3)/paletteLength, max(heatmap_data3), length.out=floor(paletteLength/2)))
pheatmap(
  as.matrix(heatmap_data3),  # 将数据框转换为矩阵
  # color =c(colorRampPalette(c("#0000FF", "white"))(length(bk)/2),
  #          colorRampPalette(colors = c("white","#CC0000"))(length(bk)/2)),
  color =myColor,
  breaks = myBreaks,
  # legend_breaks=seq(-10000,10000,5000),
  annotation_col = NULL,  # 如果有列注释，请提供相应的数据框
  annotation_row = NULL,  # 如果有行注释，请提供相应的数据框
  show_rownames = TRUE,   # 显示行名
  show_colnames = TRUE,   # 显示列名
  cluster_rows =FALSE,
  cluster_cols =FALSE,
  fontsize = font_size,   # 设置标签字体大小
  main = "Heatmap of pic-ctl Data",  # 设置标题
  filename = heatmap_file,  # 直接保存到文件
  cellwidth = 20,          # 设置单元格宽度
  cellheight = 15,         # 设置单元格高度
  border_color = "grey60",  # 单元格边框颜色
  angle_col = "90",
)

# 打印保存信息
cat(paste("Heatmap saved to", heatmap_file, "\n"))
