rm(list = ls())
# 加载所需的包
library(ggplot2)
library(dplyr)
library(pacman)
library(devtools)
library(DBI)
library(ggpubr)
library(readr)
library(writexl)
setwd("E:\\temp\\12.RNAseq\\04.BMDM_LPS\\up_isg\\collectivity")
# 定义输入文件夹和输出文件夹
input_folder <- "E:\\temp\\12.RNAseq\\04.BMDM_LPS\\up_isg\\collectivity"
output_folder <- "E:\\temp\\12.RNAseq\\04.BMDM_LPS\\up_isg\\collectivity"
# 列出文件夹中的所有文件
files <- list.files(path = input_folder, pattern = "*.txt", full.names = TRUE)

# 循环处理每个文件
for (file in files) {
  # 读取文件
  data <- read.table(file,header = T)
  dat <- cor.test(data$trna_codons_proportion,data$mrna_codons_proportion)
  # 创建散点图
  p <- ggplot(data,aes(x=trna_codons_proportion,y=mrna_codons_proportion))+ 
    geom_point(size=1,shape=15,colour = "steelblue")+
    labs(x = "trna_codons_proportion", y = "mrna_codons_proportion")+
    geom_smooth(method = lm, colour ="black")+
    theme_bw()+
    theme(axis.title = element_text(size=18),axis.text=element_text(size=13))+
    stat_cor(method = "pearson") 
  
  # 保存图形
  ggsave(
    filename = file.path(output_folder, paste0(basename(file), ".png")),
    plot = p,
    device = "png",
    width = 8,
    height = 6,
    units = "in",
    dpi = 300
  )
  ggsave(
    filename = file.path(output_folder, paste0(basename(file), ".pdf")),
    plot = p,
    device = "pdf",
    width = 6,
    height = 5,
    units = "in",
    dpi = 300
  )
}

