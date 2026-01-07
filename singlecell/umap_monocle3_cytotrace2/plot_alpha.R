pdf("train_monocle2_celltype2_1.pdf", width=10, height=7)

# 获取降维坐标并正确格式化
reduced_coords <- t(reducedDimS(cds))  # 转置矩阵
colnames(reduced_coords) <- c("Component_1", "Component_2")  # 添加列名

# 获取细胞类型信息
cell_types <- pData(cds)[["human_celltypes_250418"]]

# 创建数据框
plot_data <- data.frame(
  Component_1 = reduced_coords[, "Component_1"],
  Component_2 = reduced_coords[, "Component_2"],
  CellType = cell_types
)

# 手动创建散点图
ggplot(plot_data, aes(x = Component_1, y = Component_2, color = CellType)) +
  geom_point(size = 1, alpha = 0.1) +
  labs(color = "Cell Type") +
  theme_minimal()+ facet_wrap("~CellType", ncol = 3)

dev.off()
