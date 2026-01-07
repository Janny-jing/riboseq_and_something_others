rm(list=ls())
setwd("G:\\temp\\18.pretRNA-maturetRNA-correlation")
library(GGally)
library(corrplot)
data <- read.table("maturetRNA_isaoaccepter.txt")
ggpairs(data, upper = list(continuous = "cor"),
        lower = list(continuous = "points")) +
  theme_bw() # 使用黑色背景主题


# 计算相关矩阵
cor_matrix <- cor(data)

# 只保留上三角部分的相关系数
upper_tri <- cor_matrix
upper_tri[lower.tri(upper_tri)] <- NA

# 绘制相关系数的热力图，仅显示上三角部分
corrplot(upper_tri, method = "number", type = "upper",
         tl.col = "black", number.cex = 0.8,
         diag = FALSE, # 不显示对角线上的元素
         addCoefasPercent = TRUE, # 显示百分比形式的相关系数
         cl.pos = "n") # 隐藏颜色条
