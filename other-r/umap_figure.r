# 获取所有细胞类型的列表
sobj_dataset1$celltype_20251016 <- as.factor(sobj_dataset1$celltype_20251016)
sobj_dataset2$celltype_20251016 <- as.factor(sobj_dataset2$celltype_20251016)

all_celltypes1 <- levels(sobj_dataset1$celltype_20251016)
all_celltypes2 <- levels(sobj_dataset2$celltype_20251016)

# 获取两个数据集中所有唯一的细胞类型
all_unique_celltypes <- unique(c(all_celltypes1, all_celltypes2))

# 使用高对比度的颜色方案
get_high_contrast_colors <- function(celltypes) {
  n <- length(celltypes)
  if(n == 0) return(character(0))
  
  # 预定义一组高对比度的鲜艳颜色
  high_contrast_palette <- c(
    "#FF0000", "#00FF00", "#0000FF", "#FFFF00", "#FF00FF", "#00FFFF",
    "#FF4500", "#32CD32", "#1E90FF", "#FFD700", "#DA70D6", "#20B2AA",
    "#FF6347", "#7CFC00", "#4169E1", "#FFA500", "#BA55D3", "#00CED1",
    "#DC143C", "#9ACD32", "#0000CD", "#FF8C00", "#9932CC", "#00BFFF",
    "#B22222", "#ADFF2F", "#000080", "#FF4500", "#8A2BE2", "#1E90FF"
  )
  
  # 如果需要的颜色多于预定义的颜色，用hue_pal补充
  if(n <= length(high_contrast_palette)) {
    return(high_contrast_palette[1:n])
  } else {
    additional_colors <- scales::hue_pal()(n - length(high_contrast_palette))
    return(c(high_contrast_palette, additional_colors))
  }
}

# 创建统一的颜色映射
unified_colors <- get_high_contrast_colors(all_unique_celltypes)
names(unified_colors) <- all_unique_celltypes

# 降低其他cluster饱和度的函数 - 现在降低更多以增强对比
reduce_saturation <- function(color, factor = 0.2) {  # 降低到20%饱和度
  rgb_val <- col2rgb(color)
  hsv_val <- rgb2hsv(rgb_val)
  hsv_val[2] <- hsv_val[2] * factor
  # 同时稍微提高亮度以避免太暗
  hsv_val[3] <- min(1.0, hsv_val[3] * 1.1)
  return(hsv(hsv_val[1], hsv_val[2], hsv_val[3]))
}

# 为每个数据集创建颜色向量
colors1 <- unified_colors[all_celltypes1]
colors2 <- unified_colors[all_celltypes2]

# 突出显示指定的cluster - 使用更鲜艳的红色
highlight_color1 <- "#FF0000"  # 纯红色
highlight_color2 <- "#FF0000"  # 纯红色

target_cluster1 <- "Cholinergic, monoaminergic, peptidergic neurons"
target_cluster2 <- "Cholinergic"

# 对sobj_dataset1处理 - 大幅降低其他cluster饱和度
for(i in 1:length(colors1)) {
  cluster <- names(colors1)[i]
  if(cluster == target_cluster1) {
    colors1[i] <- highlight_color1
  } else {
    colors1[i] <- reduce_saturation(colors1[i], factor = 0.2)  # 降低到20%饱和度
  }
}

# 对sobj_dataset2处理 - 大幅降低其他cluster饱和度
for(i in 1:length(colors2)) {
  cluster <- names(colors2)[i]
  if(cluster == target_cluster2) {
    colors2[i] <- highlight_color2
  } else {
    colors2[i] <- reduce_saturation(colors2[i], factor = 0.2)  # 降低到20%饱和度
  }
}

# 确保颜色名称正确
names(colors1) <- all_celltypes1
names(colors2) <- all_celltypes2

plot1 <- DimPlot(sobj_dataset1,
                group.by = "celltype_20251016",
                pt.size = 0.3,
                label = TRUE,
                repel = TRUE,
                cols = colors1,
                raster=FALSE) +
  ggtitle("Lei Han et al. 2025") +
  theme(plot.title = element_text(size = 14, face = "bold"),
        legend.position = "right",
        legend.text = element_text(size = 8))

plot2 <- DimPlot(sobj_dataset2,
                group.by = "celltype_20251016",
                pt.size = 0.3,
                label = TRUE,
                repel = TRUE,
                cols = colors2) +
  ggtitle("Our data") +
  theme(plot.title = element_text(size = 14, face = "bold"),
        legend.position = "right",
        legend.text = element_text(size = 8))

# 组合图形并保存
combined_plot <- plot1 + plot2 +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.box = "vertical",
        legend.key.size = unit(0.4, 'cm'))

# 保存高分辨率图片
pdf("celltype_20251017_high_contrast.pdf", width=20, height=8)
print(combined_plot)
dev.off()

png("celltype_20251017_high_contrast.png", width = 5000, height = 2000, res = 300)
print(combined_plot)
dev.off()

# 极端对比版本 - 其他细胞类型变成灰色
colors1_gray <- ifelse(all_celltypes1 == target_cluster1, "#FF0000", "#D3D3D3")
colors2_gray <- ifelse(all_celltypes2 == target_cluster2, "#FF0000", "#D3D3D3")

names(colors1_gray) <- all_celltypes1
names(colors2_gray) <- all_celltypes2

plot1_gray <- DimPlot(sobj_dataset1,
                     group.by = "celltype_20251016",
                     pt.size = 0.3,
                     label = TRUE,
                     repel = TRUE,
                     cols = colors1_gray) +
  ggtitle("Lei Han et al. 2025") +
  theme(plot.title = element_text(size = 14, face = "bold"),
        legend.position = "right",
        legend.text = element_text(size = 8))

plot2_gray <- DimPlot(sobj_dataset2,
                     group.by = "celltype_20251016",
                     pt.size = 0.3,
                     label = TRUE,
                     repel = TRUE,
                     cols = colors2_gray) +
  ggtitle("Our data") +
  theme(plot.title = element_text(size = 14, face = "bold"),
        legend.position = "right",
        legend.text = element_text(size = 8))

combined_plot_gray <- plot1_gray + plot2_gray +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.box = "vertical",
        legend.key.size = unit(0.4, 'cm'))

# 保存这个版本
pdf("celltype_20251017_gray_contrast.pdf", width=20, height=8)
print(combined_plot_gray)
dev.off()