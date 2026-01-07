library(reactome.db)
library(ReactomePA)
library(openxlsx)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)
library(readxl)

setwd("G:\\temp\\23.volcano_data\\GOreactome")
path="G:\\temp\\23.volcano_data\\GOreactome"
filename=dir()
keytypes(org.Mm.eg.db)

data_1 <- read.table(filename[2],header = FALSE)
data_2 <- bitr(data_1$V1, fromType = "ENSEMBL",toType = "ENTREZID",OrgDb = org.Mm.eg.db)
data_ALL <- enrichGO(gene  = data_2$ENTREZID,
                       OrgDb         = org.Mm.eg.db,
                       ont           = "ALL" , 
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05,
                       readable      = TRUE)
res <- enrichPathway(data_2$ENTREZID, 
                       organism = "mouse",
                       minGSSize = 10, maxGSSize = 200,
                       pvalueCutoff = 0.05, pAdjustMethod = "BH")
res_genename <- setReadable(res,OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
write.xlsx(data_ALL@result,file =file.path(path, paste0(gsub("\\.xlsx$", "", basename(filename[2])), "_GO.xlsx")))
write.xlsx(res_genename@result,file =file.path(path, paste0(gsub("\\.xlsx$", "", basename(filename[2])), "_reactome.xlsx")))

xlsx_files <- list.files(path, pattern = "*.xlsx", full.names = TRUE)

# 创建一个空列表来存储每个文件的数据框
data_list <- list()

for (file in xlsx_files) {
  # 尝试读取每个Excel文件的数据
  tryCatch({
    data_list[[file]] <- read_excel(file)
    first_file_data <- data_list[[file]]
    
    # 取前15行数据并调整Description列的顺序
    subset_data <- first_file_data[1:15, ]
    subset_data$Description <- factor(subset_data$Description, levels = rev(subset_data$Description))
    
    # 设置PDF输出路径和文件名
    pdf_file <- file.path(path, paste0(tools::file_path_sans_ext(basename(file)), ".pdf"))
    
    # 开始PDF设备
    pdf(file=pdf_file, width=10, height=7)
    
    # 创建并打印气泡图
    p <- ggplot(subset_data, aes(x = FoldEnrichment, y = Description, size = Count, color = p.adjust)) +
      geom_point() +
      scale_color_gradient(low = "red", high = "green") +
      labs(title = "", x = "Fold enrichment", y = "") +
      theme_minimal()
    print(p)
    # 关闭PDF设备
    dev.off()
    cat("成功生成:", pdf_file, "\n")
  }, error = function(e) {
    cat("处理文件时出错:", file, "错误信息:", e$message, "\n")
  })
}

