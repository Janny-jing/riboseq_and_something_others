rm(list=ls())
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(scales)
library(slingshot)
library(SingleCellExperiment)
library(zellkonverter)
setwd("/storage/liuxiaodongLab/jiangjing/Projects/YutingFu/PD_YutingFu/slingshot")

sce <- readH5AD("scanpy_UMAP_n_neighbors10_npcs23_mindist0.5_spread0.5_copy.h5ad")
saveRDS(sce, "scanpy_UMAP_n_neighbors10_npcs23_mindist0.5_spread0.5_copy.rds")

data <- readRDS("scanpy_UMAP_n_neighbors10_npcs23_mindist0.5_spread0.5_copy.rds")
# # 提取 X_diffmap 的第2和第3列（注意 R 的列索引从1开始）
# new_rd <- reducedDim(data, "X_diffmap")[, 2:3]
# 
# # 将新降维结果添加到 SingleCellExperiment 对象
# reducedDim(data, "X_diffmap_2_3") <- new_rd

# #data.sce <- as.SingleCellExperiment(data)
sce_slingshot1 <- slingshot(data,      #输入单细胞对象
                            reducedDim = 'X_umap',  #降维方式
                            clusterLabels = data$human_celltypes_250418,  #cell类型
                            start.clus = 'Rgl1-like',
                            end.clus = 'DA progenitor-terminal',
                            approx_points = 150,
                            reweight = FALSE)
SlingshotDataSet(sce_slingshot1)
saveRDS(sce_slingshot1, file = "slingshot_results.rds")

sce_slingshot1 <- readRDS("slingshot_results.rds")
cell_pal <- function(cell_vars, pal_fun,...) {
  if (is.numeric(cell_vars)) {
    pal <- pal_fun(100, ...)
    return(pal[cut(cell_vars, breaks = 100)])
  } else {
    categories <- sort(unique(cell_vars))
    pal <- setNames(pal_fun(length(categories), ...), categories)
    return(pal[cell_vars])
  }
}
cell_colors <- cell_pal(sce_slingshot1$human_celltypes_250418, brewer_pal("qual", "Set2"))
plot(reducedDims(sce_slingshot1)$X_umap, col = cell_colors, pch=16, asp = 1, cex = 0.8)
lines(SlingshotDataSet(sce_slingshot1), lwd=2, col='black')

library(tradeSeq)
# fit negative binomial GAM
sce <- fitGAM(sce_slingshot1)
# test for dynamic expression
ATres <- associationTest(sce_slingshot1)
topgenes <- rownames(ATres[order(ATres$pvalue), ])[1:250]
pst.ord <- order(sce$slingPseudotime_1, na.last = NA)
heatdata <- assays(sce)$counts[topgenes, pst.ord]
heatclus <- sce$GMM[pst.ord]

heatmap(log1p(heatdata), Colv = NA,
        ColSideColors = brewer.pal(9,"Set1")[heatclus])
















# 步骤1：生成初始 Lineages
lineages <- getLineages(
  reducedDim(data, "X_diffmap_2_3"),
  clusterLabels = data$human_celltypes_250418,
  start.clus = "Rgl1-like"
)

# 提取初始邻接矩阵
adj_matrix <- as.matrix(lineages@metadata$mst)
# 清除所有自动生成的连接
adj_matrix[] <- 0  # 重置为全0

# 定义 Lineage1 的精确连接路径（无分支）
clusters1 <- c("Rgl1-like", 
               "DA progenitor-FABP7",
               "DA progenitor-CHRNB3 low", 
               "DA progenitor-terminal")

# 定义 Lineage2 的精确连接路径（无交叉）
clusters2 <- c("Rgl1-like", 
               "DA progenitor-FABP7",
               "DA progenitor-CHRNB3 high", 
               "Mature DA",
               "Immature DA")

# 填充邻接矩阵 (仅允许单链连接)
set_connections <- function(clusters, adj) {
  for(i in 1:(length(clusters)-1)){
    from <- clusters[i]
    to <- clusters[i+1]
    adj[from, to] <- 1
    adj[to, from] <- 1  # 确保双向连接
  }
  return(adj)
}

adj_matrix <- set_connections(clusters1, adj_matrix)
adj_matrix <- set_connections(clusters2, adj_matrix)


# 更新 MST 结构
lineages@metadata$mst <- adj_matrix

curves <- getCurves(
  lineages,
  approx_points = 150,
  
  # 关键参数：关闭所有自动调整
  shrink = FALSE,      # 禁止分叉收缩
  extend = "n",        # 禁止延伸端点
  reweight = FALSE,    # 禁止重新分配权重
  reassign = FALSE     # 禁止重新分配细胞到轨迹
)

sce_slingshot1 <- curves

# 查看最终 Lineages
slingLineages(sce_slingshot1)

# 如果仍然错误，直接覆盖路径（仅限 SlingshotDataSet 对象）
if ("SlingshotDataSet" %in% class(sce_slingshot1)) {
  sce_slingshot1@metadata[["lineages"]] <- list(
    clusters1,
    clusters2
  )
}


# 提取降维坐标和聚类标签
plot_data <- data.frame(
  `Dim-1` = reducedDim(data, "X_diffmap_2_3")[,1],
  `Dim-2` = reducedDim(data, "X_diffmap_2_3")[,2],
  Cluster = data$human_celltypes_250418
)
# 重命名列（避免特殊符号问题）
colnames(plot_data)[1:2] <- c("Dim_1", "Dim_2")  # 转换为合法变量名

# 提取轨迹曲线并统一列名
curves <- slingCurves(sce_slingshot1)
for (i in seq_along(curves)) {
  colnames(curves[[i]]$s) <- c("Dim_1", "Dim_2")  # 统一命名为 Dim_1/Dim_2
}
# 创建颜色映射
lineage_colors <- c("black", "black")  # 红/蓝对应 Lineage1/2
names(lineage_colors) <- names(curves)

# 绘制基础散点图
p <- ggplot(plot_data, aes(x = Dim_1, y = Dim_2)) +
  geom_point(aes(color = Cluster), size = 1.5, alpha = 0.6) +
  scale_color_manual(values = brewer.pal(8, "Set2")) +  # 自定义聚类颜色
  theme_classic(base_size = 14)

# 叠加轨迹曲线
for (i in seq_along(curves)) {
  curve_df <- as.data.frame(curves[[i]]$s)
  p <- p + 
    geom_path(data = curve_df, 
              aes(x = Dim_1, y = Dim_2), 
              color = lineage_colors[i], 
              linewidth = 1.5, 
              linetype = "solid") +
    geom_point(data = curve_df[c(1, nrow(curve_df)), ],  # 标记起点和终点
               aes(x = Dim_1, y = Dim_2), 
               color = lineage_colors[i],
               size = 3, shape = 18)
}

# 显示图形
print(p)

library(gam)
library(tradeSeq)
pseudotime <- sce_slingshot1@assays@data@listData[["pseudotime"]]

pt <- pseudotime[,1] # Using first lineage
pt <- pt[!is.na(pt)] # Remove NA values


if (!"logcounts" %in% assayNames(data)) {
  assay(data, "logcounts") <- log1p(assay(data, "X"))
}

# Verify gene names are properly set
rownames(data)[1:5]  # Check first few genes

# 2.3 Prepare expression matrix
expr <- assay(data, "logcounts") 
expr <- expr[, names(pt)]

# 4.2 Fit GAM models for all genes
gam_results <- apply(expr, 1, function(y) {
  d <- data.frame(y = y, t = pt)
  tryCatch({
    gam.fit <- gam(y ~ lo(t), data = d)
    pval <- summary(gam.fit)$anova$`Pr(F)`[2]
    return(c(pval = pval, dev.expl = summary(gam.fit)$dev.expl))
  }, error = function(e) return(c(pval = NA, dev.expl = NA)))
})

# Convert to data frame
gam_df <- as.data.frame(t(gam_results))
gam_df$gene <- rownames(gam_df)
gam_df <- gam_df[complete.cases(gam_df), ]

# 4.3 Select top 200 significant genes
top_genes <- gam_df[order(gam_df$pval), "gene"][1:200]

# 4.4 Smooth expression (20-cell window as per paper)
smooth_expr <- t(apply(expr[top_genes, order(pt)], 1, function(x) {
  predict(loess(x ~ seq_along(x), span = 20/length(x)))
}))

# 4.5 Plot heatmap
library(pheatmap)
pheatmap(
  smooth_expr,
  cluster_cols = FALSE,
  show_colnames = FALSE,
  color = viridis(100),
  main = "Temporal Gene Expression (Top 200 Genes)"
)
