# 获取降维后的坐标（正确方式）
reduced_coords <- t(reducedDimS(cds))
# 确保 cell names 在 pData(cds) 和 reduced_coords 中一致
cell_names_reduced <- rownames(reduced_coords)
cell_names_pdata <- rownames(pData(cds))
if (!all(cell_names_reduced == cell_names_pdata)) {
  stop("细胞顺序不一致，请先排序或子集筛选")
}

# 筛选 DA progenitor-CHRNB3 low 的细胞
chrnb3_low_cells <- pData(cds)$human_celltypes_250418 == "DA progenitor-CHRNB3 low"
# 提取这些细胞的坐标
chrnb3_low_coords <- reduced_coords[chrnb3_low_cells, , drop = FALSE]
# 筛选左上角的细胞
left_upper_cells <- chrnb3_low_coords[, 1] < -15

# 提取左上角细胞的坐标和名字
left_upper_coords <- chrnb3_low_coords[left_upper_cells, , drop = FALSE]
left_upper_cell_names <- rownames(left_upper_coords)

# 导出到 CSV 文件，方便导入 Loupe 或其他工具
write.csv(left_upper_coords, file = "left_upper_chrnb3_low_cells.csv", row.names = TRUE)

print(paste("共找到", sum(left_upper_cells), "个位于左上角的细胞"))


