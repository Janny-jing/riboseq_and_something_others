
library(GenomicFeatures)
library(GenomicAlignments)
library(rtracklayer)
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)

# 设置文件路径
gtf_file <- "Y:/temp/20251122-riboseq_new_3sample/gencode.v48.annotation.gtf"
bam_files <- list.files(pattern = "\\.bam$")
sample_names <- gsub("\\.bam$", "", bam_files)

# 自动识别分组
groups <- ifelse(grepl("orf3a", sample_names), "ORF3a", "EV")

## 第一步：修复的起始密码子提取函数
cat("提取起始密码子位置...\n")

extract_start_codons <- function(gtf_file) {
  # 使用txdbmaker包来避免警告
  if (!require("txdbmaker", quietly = TRUE)) {
    BiocManager::install("txdbmaker")
  }
  library(txdbmaker)
  
  # 创建TxDb对象
  txdb <- makeTxDbFromGFF(gtf_file, format = "gtf")
  
  # 获取所有编码转录本
  cds <- cdsBy(txdb, by = "tx", use.names = TRUE)
  
  # 过滤掉没有CDS的转录本
  cds <- cds[lengths(cds) > 0]
  
  cat("找到", length(cds), "个有CDS的转录本\n")
  
  # 提取起始密码子
  start_codons <- GRangesList()
  gene_info <- data.frame()
  
  for (tx_id in names(cds)) {
    tx_cds <- cds[[tx_id]]
    
    # 根据链的方向确定起始密码子
    if (as.character(strand(tx_cds[1])) == "+") {
      start_codon <- resize(tx_cds[1], width = 3, fix = "start")
    } else {
      start_codon <- resize(tx_cds[length(tx_cds)], width = 3, fix = "end")
    }
    
    # 获取基因信息
    tx_info <- transcripts(txdb, filter = list(tx_id = tx_id))
    
    start_codons[[tx_id]] <- start_codon
    gene_info <- rbind(gene_info, data.frame(
      transcript_id = tx_id,
      gene_id = ifelse(length(tx_info$gene_id) > 0, tx_info$gene_id, tx_id),
      chr = as.character(seqnames(start_codon)),
      start = start(start_codon),
      end = end(start_codon),
      strand = as.character(strand(start_codon)),
      stringsAsFactors = FALSE
    ))
  }
  
  return(list(start_codons = start_codons, gene_info = gene_info))
}

start_data <- extract_start_codons(gtf_file)
start_codons <- start_data$start_codons
gene_info <- start_data$gene_info

## 第二步：简化的meta-gene分析方法
cat("使用简化但更稳定的方法...\n")

simple_metagene_analysis <- function(bam_files, sample_names, groups, start_codons) {
  # 为所有起始密码子创建统一的窗口
  upstream <- 50
  downstream <- 100
  window_size <- upstream + downstream + 1
  positions <- -upstream:downstream
  
  # 初始化结果矩阵
  all_results <- data.frame()
  
  for (i in 1:length(bam_files)) {
    cat("处理样本:", sample_names[i], "...\n")
    
    # 读取BAM文件
    reads <- readGAlignments(bam_files[i])
    reads_gr <- granges(reads)
    
    # 初始化该样本的结果
    sample_results <- data.frame()
    
    # 对每个起始密码子计算覆盖度
    for (tx_id in names(start_codons)) {
      start_codon <- start_codons[[tx_id]]
      
      # 创建分析窗口
      if (as.character(strand(start_codon)) == "+") {
        window_start <- start(start_codon) - upstream
        window_end <- start(start_codon) + downstream
      } else {
        window_start <- start(start_codon) - downstream
        window_end <- start(start_codon) + upstream
      }
      
      window_gr <- GRanges(
        seqnames = seqnames(start_codon),
        ranges = IRanges(start = window_start, end = window_end),
        strand = strand(start_codon)
      )
      
      # 找到窗口内的reads
      overlaps <- findOverlaps(reads_gr, window_gr)
      if (length(overlaps) > 0) {
        window_reads <- reads_gr[queryHits(overlaps)]
        
        # 计算相对位置
        if (as.character(strand(start_codon)) == "+") {
          rel_pos <- start(window_reads) - start(start_codon)
        } else {
          rel_pos <- start(start_codon) - start(window_reads)
        }
        
        # 过滤有效位置
        valid_idx <- rel_pos >= -upstream & rel_pos <= downstream
        rel_pos <- rel_pos[valid_idx]
        
        if (length(rel_pos) > 0) {
          # 计算frame
          frames <- (rel_pos - 1) %% 3
          
          # 汇总计数
          for (pos in unique(rel_pos)) {
            pos_frames <- frames[rel_pos == pos]
            frame_counts <- table(pos_frames)
            
            for (frame in names(frame_counts)) {
              temp_df <- data.frame(
                sample = sample_names[i],
                group = groups[i],
                transcript_id = tx_id,
                position = pos,
                frame = paste0("Frame", frame),
                count = as.numeric(frame_counts[frame]),
                stringsAsFactors = FALSE
              )
              sample_results <- rbind(sample_results, temp_df)
            }
          }
        }
      }
    }
    
    all_results <- rbind(all_results, sample_results)
  }
  
  return(all_results)
}

# 执行简化的meta-gene分析
metagene_data <- simple_metagene_analysis(bam_files, sample_names, groups, start_codons)

## 第三步：数据标准化和汇总
cat("标准化和汇总数据...\n")

# 计算每个样本的总reads用于标准化
sample_totals <- metagene_data %>%
  group_by(sample) %>%
  summarise(total_reads = sum(count), .groups = 'drop')

metagene_normalized <- metagene_data %>%
  left_join(sample_totals, by = "sample") %>%
  mutate(normalized_count = (count / total_reads) * 1e6)  # RPM标准化

# 按组别、位置和frame汇总
metagene_summary <- metagene_normalized %>%
  group_by(group, position, frame) %>%
  summarise(
    mean_reads = mean(normalized_count, na.rm = TRUE),
    se_reads = sd(normalized_count, na.rm = TRUE) / sqrt(n()),
    n_genes = n_distinct(transcript_id),
    .groups = 'drop'
  )

## 第四步：绘制meta-gene图
cat("绘制meta-gene图...\n")

# 修改分组名称用于显示
metagene_summary$group_display <- ifelse(metagene_summary$group == "EV", "EV", "ORF3a")

p_metagene <- ggplot(metagene_summary, aes(x = position, y = mean_reads, color = frame)) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  facet_wrap(~ group_display, ncol = 1) +
  scale_color_manual(values = c("Frame0" = "#E41A1C", "Frame1" = "#377EB8", "Frame2" = "#4DAF4A")) +
  labs(title = "Meta-gene Analysis of Ribo-seq Density around Start Codon",
       subtitle = paste("Analysis of", length(start_codons), "coding transcripts"),
       x = "Distance from aTIS (nt)",
       y = "Normalized mean reads (RPM)",
       color = "Reading Frame") +
  theme_classic() +
  theme(
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold", size = 12),
    axis.text = element_text(color = "black"),
    legend.position = "top"
  ) +
  scale_x_continuous(breaks = seq(-50, 100, by = 25))

print(p_metagene)

# 保存图表
ggsave("metagene_analysis.pdf", p_metagene, width = 10, height = 8, bg = "white")
ggsave("metagene_analysis.png", p_metagene, width = 10, height = 8, dpi = 300, bg = "white")

## 第五步：保存结果和统计摘要
write.csv(metagene_summary, "metagene_analysis_results.csv", row.names = FALSE)

# 生成统计摘要
cat("\n=== Meta-gene分析统计摘要 ===\n")
cat("分析样本:", paste(sample_names, collapse = ", "), "\n")
cat("分组:", paste(unique(groups), collapse = " vs "), "\n")
cat("分析的转录本数量:", length(start_codons), "\n")
cat("分析窗口: -50 to +100 nt around start codon\n")

# 检查起始密码子周围的frame分布
start_region <- metagene_summary %>%
  filter(position >= -5 & position <= 5) %>%
  group_by(group, frame) %>%
  summarise(
    avg_density = mean(mean_reads, na.rm = TRUE),
    .groups = 'drop'
  )

cat("\n起始密码子周围frame分布 (-5 to +5 nt):\n")
print(start_region)

cat("\n分析完成！\n")
cat("生成的文件:\n")
cat("- metagene_analysis.pdf\n")
cat("- metagene_analysis.png\n")
cat("- metagene_analysis_results.csv\n")