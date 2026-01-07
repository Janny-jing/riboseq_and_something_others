rm(list = ls())
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
require(plyr)
library(edgeR)
library(clusterProfiler)
library(org.Hs.eg.db)  # 根据您的物种修改，如 org.Mm.eg.db 用于小鼠
library(enrichplot)
library(DOSE)
library(ggrepel)  # 用于火山图标签

setwd("C:\\Users\\jing.jiang\\Desktop\\20251017_ct_rnaseq_dif")
filename <- dir()
data<-read.table(filename[1],header = T,row.names = 1)
data <- data[,6:27]
expr <- na.omit(data)
expr <- expr[apply(expr,1,max)>=10, ]
dup_gene<-data.frame(gene=expr$gene_name,mean=apply(expr[,-grep("gene_name",colnames(expr))],1,mean))
expr <-expr[order(dup_gene$gene,dup_gene$mean,decreasing = T),]
expr<- expr[!duplicated(expr$gene_name),]
rownames(expr) <- expr$gene_name
expr <- expr[,-grep("gene_name",colnames(expr))]
colnames(expr) <- c("0_R1","0_R2","0_R3","10_R1","10_R2","10_R3","2_R1","2_R2","2_R3","4_R1","4_R2","4_R3",
                    "6_R1","6_R2","6_R3","8_R1","8_R2","8_R3","mock_R1","mock_R2","mock_R3")
col.data = data.frame(Sample = rep(c("0","10","2","4","6","8","mock"), c(3,3,3,3,3,3,3)))
dds <- DESeqDataSetFromMatrix(countData = expr, colData = col.data, design = ~ Sample)

####### DE analysis - 循环处理不同时间点 ########
design(dds) = ~Sample
dds = DESeq(dds)

# 定义要比较的时间点
time_points <- c("0", "2", "4", "6", "8", "10")
control_group <- "mock"

# 创建输出目录
dir.create("../deseq2_timecourse", showWarnings = FALSE)
setwd("../deseq2_timecourse")

# 循环分析每个时间点
for(time_point in time_points) {
  cat("正在分析时间点:", time_point, "vs", control_group, "\n")
  
  # 创建时间点特定的输出目录
  time_dir <- paste0("time_", time_point, "_vs_", control_group)
  dir.create(time_dir, showWarnings = FALSE)
  setwd(time_dir)
  
  ####### 差异表达分析 #######
  res <- results(dds, contrast = c("Sample", time_point, control_group))
  
  # 保存原始结果
  write.csv(res, file = paste0("DEanalysis_", time_point, "_vs_", control_group, ".csv"))
  
  # 获取DEGs
  up.lst <- row.names(res[!is.na(res$log2FoldChange) & !is.na(res$padj) & 
                           res$log2FoldChange > 1 & res$padj < 0.05,])
  dn.lst <- row.names(res[!is.na(res$log2FoldChange) & !is.na(res$padj) & 
                           res$log2FoldChange < -1 & res$padj < 0.05,])
  
  cat("上调基因数量:", length(up.lst), "\n")
  cat("下调基因数量:", length(dn.lst), "\n")
  
  # 保存DEGs列表
  write.table(sub(":.*$","",up.lst), 
              file=paste0("DEGs_Up_", time_point, "_vs_", control_group, ".txt"), 
              row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
  write.table(sub(":.*$","",dn.lst), 
              file=paste0("DEGs_Down_", time_point, "_vs_", control_group, ".txt"), 
              row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
  
  # 按p值排序并保存
  res.srt = res[order(res$pvalue),]
  write.csv(res.srt, file = paste0("DE_sorted_", time_point, "_vs_", control_group, ".csv"))
  write.csv(res[!is.na(res$log2FoldChange) & !is.na(res$padj) & 
                 res$log2FoldChange > 1 & res$padj < 0.05,], 
            file = paste0("DEGs_Up_sorted_", time_point, "_vs_", control_group, ".csv"))
  write.csv(res[!is.na(res$log2FoldChange) & !is.na(res$padj) & 
                 res$log2FoldChange < -1 & res$padj < 0.05,], 
            file = paste0("DEGs_down_sorted_", time_point, "_vs_", control_group, ".csv"))
  
  de.lst = unique(c(up.lst,dn.lst))
  
  ####### 火山图 #######
  # 准备火山图数据
  volcano_data <- as.data.frame(res)
  volcano_data$gene <- rownames(volcano_data)
  volcano_data$gene_symbol <- sub(":.*$", "", volcano_data$gene)  # 提取基因符号
  
  # 添加显著性分类
  volcano_data$diffexpressed <- "NO"
  volcano_data$diffexpressed[volcano_data$log2FoldChange > 1 & volcano_data$padj < 0.05] <- "UP"
  volcano_data$diffexpressed[volcano_data$log2FoldChange < -1 & volcano_data$padj < 0.05] <- "DOWN"
  
  # 添加标签（选择最显著的基因进行标注）
  volcano_data$delabel <- NA
  # 选择padj最小的前10个上调基因和前10个下调基因进行标注
  top_up <- volcano_data[volcano_data$diffexpressed == "UP", ]
  top_down <- volcano_data[volcano_data$diffexpressed == "DOWN", ]
  
  if(nrow(top_up) > 0) {
    top_up <- top_up[order(top_up$padj), ]
    top_up_genes <- head(top_up$gene_symbol, 10)
    volcano_data$delabel[volcano_data$gene_symbol %in% top_up_genes] <- top_up_genes
  }
  
  if(nrow(top_down) > 0) {
    top_down <- top_down[order(top_down$padj), ]
    top_down_genes <- head(top_down$gene_symbol, 10)
    volcano_data$delabel[volcano_data$gene_symbol %in% top_down_genes] <- top_down_genes
  }
  
  # 设置颜色
  my_colors <- c("DOWN" = "blue", "NO" = "grey", "UP" = "red")
  
  # PDF格式火山图
  pdf(file = paste0("Volcano_", time_point, "_vs_", control_group, ".pdf"), width = 10, height = 8)
  p <- ggplot(data = volcano_data, 
              aes(x = log2FoldChange, 
                  y = -log10(padj), 
                  col = diffexpressed, 
                  label = delabel)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = my_colors) +
    geom_vline(xintercept = c(-1, 1), col = "black", linetype = "dashed", alpha = 0.5) +
    geom_hline(yintercept = -log10(0.05), col = "black", linetype = "dashed", alpha = 0.5) +
    labs(
      title = paste("火山图 -", time_point, "vs", control_group),
      x = expression("Log"[2]*" Fold Change"),
      y = expression("-Log"[10]*" Adjusted P-value"),
      color = "Expression"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      legend.position = "right"
    )
  
  # 添加基因标签（如果有显著基因）
  if(sum(!is.na(volcano_data$delabel)) > 0) {
    p <- p + geom_text_repel(max.overlaps = 20, size = 3)
  }
  
  print(p)
  dev.off()
  
  # TIFF格式火山图
  tiff(file = paste0("Volcano_", time_point, "_vs_", control_group, ".tiff"), 
       width = 20, height = 16, units = "cm", res = 600, compression = "lzw")
  print(p)
  dev.off()
  
  # 简化版火山图（不带标签，用于快速查看）
  pdf(file = paste0("Volcano_simple_", time_point, "_vs_", control_group, ".pdf"), width = 8, height = 6)
  p_simple <- ggplot(data = volcano_data, 
                     aes(x = log2FoldChange, 
                         y = -log10(padj), 
                         col = diffexpressed)) +
    geom_point(alpha = 0.6, size = 1) +
    scale_color_manual(values = my_colors) +
    geom_vline(xintercept = c(-1, 1), col = "black", linetype = "dashed", alpha = 0.5) +
    geom_hline(yintercept = -log10(0.05), col = "black", linetype = "dashed", alpha = 0.5) +
    labs(
      title = paste("火山图 -", time_point, "vs", control_group),
      x = expression("Log"[2]*" Fold Change"),
      y = expression("-Log"[10]*" Adjusted P-value"),
      color = "Expression"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      legend.position = "right"
    )
  print(p_simple)
  dev.off()
  
  ####### 富集分析 - 只有当有足够DEGs时才进行 #######
  if(length(de.lst) >= 5) {
    # 转换基因ID (假设是ENSEMBL ID)
    gene_ids <- sub(":.*$", "", de.lst)
    
    # GO富集分析
    tryCatch({
      ego <- enrichGO(gene          = gene_ids,
                      OrgDb         = org.Hs.eg.db,
                      keyType       = "ENSEMBL",
                      ont           = "ALL",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.2,
                      readable      = TRUE)
      
      if(nrow(ego) > 0) {
        # 保存GO结果
        write.csv(ego, file = paste0("GO_enrichment_", time_point, "_vs_", control_group, ".csv"))
        
        # GO点图
        pdf(file = paste0("GO_dotplot_", time_point, "_vs_", control_group, ".pdf"), width=10, height=8)
        print(dotplot(ego, showCategory=15))
        dev.off()
        
        # GO条形图
        pdf(file = paste0("GO_barplot_", time_point, "_vs_", control_group, ".pdf"), width=10, height=8)
        print(barplot(ego, showCategory=15))
        dev.off()
      }
    }, error = function(e) {
      cat("GO富集分析出错:", e$message, "\n")
    })
    
    # KEGG富集分析
    tryCatch({
      # 需要将ENSEMBL ID转换为ENTREZID
      gene_df <- bitr(gene_ids, fromType = "ENSEMBL", 
                      toType = "ENTREZID", 
                      OrgDb = org.Hs.eg.db)
      
      if(nrow(gene_df) > 0) {
        kk <- enrichKEGG(gene         = gene_df$ENTREZID,
                         organism     = 'hsa',  # 根据物种修改，小鼠用'mmu'
                         pvalueCutoff = 0.05)
        
        if(nrow(kk) > 0) {
          # 保存KEGG结果
          write.csv(kk, file = paste0("KEGG_enrichment_", time_point, "_vs_", control_group, ".csv"))
          
          # KEGG点图
          pdf(file = paste0("KEGG_dotplot_", time_point, "_vs_", control_group, ".pdf"), width=10, height=8)
          print(dotplot(kk, showCategory=15))
          dev.off()
        }
      }
    }, error = function(e) {
      cat("KEGG富集分析出错:", e$message, "\n")
    })
    
    # GSEA分析
    tryCatch({
      # 准备GSEA输入数据
      res_df <- as.data.frame(res)
      res_df <- res_df[!is.na(res_df$padj), ]
      res_df <- res_df[order(res_df$log2FoldChange, decreasing = TRUE), ]
      
      # 提取基因列表和log2FC
      geneList <- res_df$log2FoldChange
      names(geneList) <- sub(":.*$", "", rownames(res_df))
      
      # 去除重复基因名
      geneList <- geneList[!duplicated(names(geneList))]
      
      # GO GSEA
      gsea_go <- gseGO(geneList     = geneList,
                       OrgDb        = org.Hs.eg.db,
                       ont          = "BP",
                       keyType      = "ENSEMBL",
                       minGSSize    = 10,
                       maxGSSize    = 500,
                       pvalueCutoff = 0.05,
                       verbose      = FALSE)
      
      if(nrow(gsea_go) > 0) {
        write.csv(gsea_go, file = paste0("GSEA_GO_", time_point, "_vs_", control_group, ".csv"))
        
        # GSEA图
        pdf(file = paste0("GSEA_plot_", time_point, "_vs_", control_group, ".pdf"), width=10, height=8)
        if(nrow(gsea_go) >= 3) {
          print(dotplot(gsea_go, showCategory=15))
        } else {
          print(dotplot(gsea_go, showCategory=nrow(gsea_go)))
        }
        dev.off()
        
        # 绘制特定通路的GSEA图
        if(nrow(gsea_go) >= 1) {
          pdf(file = paste0("GSEA_enrichment_plot_", time_point, "_vs_", control_group, ".pdf"), 
              width=8, height=6)
          print(gseaplot2(gsea_go, geneSetID = 1, title = gsea_go$Description[1]))
          dev.off()
        }
      }
    }, error = function(e) {
      cat("GSEA分析出错:", e$message, "\n")
    })
  } else {
    cat("DEGs数量不足，跳过富集分析\n")
  }
  
  cat("完成时间点:", time_point, "的分析\n\n")
  setwd("..")  # 返回上级目录
}

cat("所有时间点分析完成！\n")