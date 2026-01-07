library(Seurat)
library(ggplot2)
library(patchwork)

data <- readRDS("preprocess_merge_data_last_harmony.rds")
rename_map <- list(
  "Astrocytes" = "Astrocyte",
  "AC" = "Astrocyte",
  "EPEN" = "Ependymal",
  "OLG" = "Oligodendrocyte",
  "Vascular cells" = "Endothelial",
  "TEINI-ChaT" = "Cholinergic",
  "TEINI-Gpc3" = "GABA-1",
  "TEINI-Hpse" = "GABA-2",
  "TEINI-Il1rapl2" = "GABA-3",
  "TEINI-Meis2" = "GABA-4",
  "TEINI-Sst" = "GABA-5",
  "TEINI-Vip" = "GABA-6",
  "TEGLU-Abi3bp" = "Glut-1",
  "TEGLU-ETV1" = "Glut-2",
  "TEGLU-Grp" = "Glut-3",
  "TEGLU-Nr4a2" = "Glut-4",
  "TEGLU-Nxph3" = "Glut-5",
  "TEGLU-Rxfp1" = "Glut-6",
  "TEGLU-Slc30a3" = "Glut-7",
  "MSN-Drd1" = "D1-SPN",
  "MSN-Drd2" = "D2-SPN",
  "Endothelial cells" = "Endothelial",
  "Ependymal cells" = "Ependymal",
  "Hypendymal cells" = "Ependymal",
  "Oligodendrocytes" = "Oligodendrocyte",
  "Oligodendrocyte precursor cells" = "OPC",
  "Telencephalon excitatory neurons" = "Glut",
  "Telencephalon inhibitory neurons" = "GABA",
  "Vascular and leptomeningeal cells" = "VLMC"
)

data@meta.data$celltype_20250930 <- sapply(as.character(data@meta.data$celltype), function(x) {
  if (x %in% names(rename_map)) {
    return(rename_map[[x]])
  } else {
    return(x)  # 保留未在映射表中的原始名称
  }
})

neuron_clusters <- c(
  "Astrocyte", "Ependymal", "Oligodendrocyte", "Endothelial", "Cholinergic","GABA-1","GABA-2","GABA-3","GABA-4","GABA-5",
  "GABA-6", "Glut-1", "Glut-2", "Glut-3","Glut-4","Glut-5","Glut-6","Glut-7","D1-SPN","D2-SPN","Glut",
  "GABA", "OPC", "VLMC", "Microglia"
)

selected_cells <- data@meta.data[data@meta.data$celltype_20250930 %in% neuron_clusters, ]

# 创建子集
sobj_sub <- subset(data, cells = rownames(selected_cells))

# 验证结果
table(sobj_sub@meta.data$celltype_20250930)

# 安装必要包
if (!require("patchwork")) install.packages("patchwork")
library(Seurat)
library(ggplot2)
library(patchwork)

# 步骤1：提取两个数据集
sobj_dataset1 <- subset(sobj_sub, subset = dataset == "mouseBrain_308Clusters")
sobj_dataset2 <- subset(sobj_sub, subset = dataset == "Our_data")
shared_colors <- scales::hue_pal()(length(unique(sobj_sub$celltype_20250930)))

# 步骤2：为每个数据集创建UMAP图
plot1 <- DimPlot(sobj_dataset1,
                group.by = "celltype_20250930",
                pt.size = 0.6,
                label = TRUE,
                repel = TRUE) +
  ggtitle("Lei Han et al. 2025") +
  theme(plot.title = element_text(size = 14, face = "bold"),
        legend.position = "right",
        legend.text = element_text(size = 8))+
        scale_color_manual(values = shared_colors)

plot2 <- DimPlot(sobj_dataset2,
                group.by = "celltype_20250930",
                pt.size = 0.6,
                label = TRUE,
                repel = TRUE) +
  ggtitle("Our data") +
  theme(plot.title = element_text(size = 14, face = "bold"),
        legend.position = "right",
        legend.text = element_text(size = 8))+
        scale_color_manual(values = shared_colors)

# 步骤3：组合图形并保存
combined_plot <- plot1 + plot2 +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.box = "vertical",
        legend.key.size = unit(0.4, 'cm'))

# 保存高分辨率图片
pdf("celltype_20250930.pdf",width=20,height=8)
png("celltype_20250930.png", width = 5000, height = 1500, res = 300)
print(combined_plot)
dev.off()

png("umap_plot_celltype_20250930_1.png", width = 3000, height = 1500, res = 300)
DimPlot(sobj_sub,
        group.by = "celltype_20250930",
        label = FALSE,
        pt.size = 0.6)
dev.off()
