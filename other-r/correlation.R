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
setwd("E:\\temp\\12.RNAseq\\03.BMDM_polyic\\de_cor\\dn_de\\eachgene")
# 定义输入文件夹和输出文件夹
dir.create("E:\\temp\\12.RNAseq\\03.BMDM_polyic\\de_cor\\dn_de\\eachgene\\eachgene_ctl_result")
dir.create("E:\\temp\\12.RNAseq\\03.BMDM_polyic\\de_cor\\dn_de\\eachgene\\eachgene_ctl_result_P")
dir.create("E:\\temp\\12.RNAseq\\03.BMDM_polyic\\de_cor\\dn_de\\eachgene\\eachgene_polyic_result")
dir.create("E:\\temp\\12.RNAseq\\03.BMDM_polyic\\de_cor\\dn_de\\eachgene\\eachgene_polyic_result_P")
input_folder <- "E:\\temp\\12.RNAseq\\03.BMDM_polyic\\de_cor\\dn_de\\eachgene\\eachgene_ctl"
output_folder <- "E:\\temp\\12.RNAseq\\03.BMDM_polyic\\de_cor\\dn_de\\eachgene\\eachgene_ctl_result"
output_folder_0.05 <- "E:\\temp\\12.RNAseq\\03.BMDM_polyic\\de_cor\\dn_de\\eachgene\\eachgene_ctl_result_P"

# 列出文件夹中的所有文件
files <- list.files(path = input_folder, pattern = "*.txt", full.names = TRUE)

# 创建一个空的数据框来存储结果
results_df <- data.frame(
  gene = character(),
  Correlation = numeric(),
  P_Value = numeric(),
  stringsAsFactors = FALSE
)

# 循环遍历文件列表中的每对文件组合
for (i in 1:length(files) ) {
  # 读取两个文件的内容
  file_data <- read.table(files[i],header = T)
  # 计算相关性及p值
  dat <- cor.test(file_data$trna_codons_proportion,file_data$mrna_codons_proportion)
  # 存储结果
  new_row <- data.frame(
    gene = gsub("\\.txt$", "", basename(files[i])),
    Correlation = dat$estimate,
    P_Value = dat$p.value
  )
  
  results_df <- rbind(results_df, new_row)
}

# 输出结果到Excel文件
write_xlsx(results_df, "Correlation_Results_polyic_ctl.xlsx")

# 查看结果
print(results_df)

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

# 循环处理每个文件
for (file in files) {
  # 读取文件
  data <- read.table(file,header = T)
  dat <- cor.test(data$trna_codons_proportion,data$mrna_codons_proportion)
  if (dat$p.value < 0.05) {
    output_file <- file.path(output_folder_0.05,basename(file))
    write.table(data,output_file)
    p <- ggscatter(data, x = "trna_codons_proportion", y = "mrna_codons_proportion",
                   add = "reg.line", cor.coef = TRUE, cor.method = "pearson",
                   xlab = "trna_codons_proportion", ylab = "mrna_codons_proportion",
                   title = paste0("Scatter Plot for ", basename(file)),
                   conf.int = TRUE, cor.coef.label = TRUE, cor.coef.args = list(size = 3),
                   test.name = "pearson", test.args = list(alternative = "two.sided"),
                   test.label = TRUE, test.label.x = 0.4, test.label.y = 0.9)
    
    # 保存图形
    ggsave(
      filename = file.path(output_folder_0.05, paste0(basename(file), ".png")),
      plot = p,
      device = "png",
      width = 8,
      height = 6,
      units = "in"
    )
  }
}


