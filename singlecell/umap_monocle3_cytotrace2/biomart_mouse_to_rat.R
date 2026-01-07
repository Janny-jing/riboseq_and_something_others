library(Seurat)
library(homologene)
public_data <- readRDS("Ximerakis.rds")
mouse_genes <- rownames(public_data)
homolog_table <- homologene(mouse_genes, inTax = 10090, outTax = 10116)
homologs <- data.frame(
  mouse_gene = homolog_table$`10090`,
  rat_gene = homolog_table$`10116`
)

# 去除空值和重复
homologs <- homologs[homologs$rat_gene != "" & !is.na(homologs$rat_gene), ]
homologs_dedup <- homologs[!duplicated(homologs$mouse_gene), ]


gene_conversion <- setNames(
  homologs_dedup$rat_gene,
  homologs_dedup$mouse_gene
)
new_genes <- ifelse(
  mouse_genes %in% names(gene_conversion),
  gene_conversion[mouse_genes],  # 转换后的基因
  mouse_genes                   # 未转换的基因（保留原名称）
)

sc_data_converted <- public_data

# 2. 更新基因名
rownames(sc_data_converted@assays$RNA@counts) <- new_genes
rownames(sc_data_converted@assays$RNA@data) <- new_genes  # 确保所有层面都更新

# 3. 安全移除无效/重复基因（逐步处理）
# 步骤1：移除空基因名
empty_mask <- (new_genes == "") | is.na(new_genes)
if(any(empty_mask)) {
  sc_data_converted <- sc_data_converted[!empty_mask, ]
  new_genes <- new_genes[!empty_mask]  # 同步更新基因名向量
}

# 步骤2：处理重复基因（保留表达量最高的）
duplicate_genes <- unique(new_genes[duplicated(new_genes)])

if(length(duplicate_genes) > 0) {
  # 创建保留索引
  keep_indices <- rep(TRUE, length(new_genes))

  for(gene in duplicate_genes) {
    dup_idx <- which(new_genes == gene)

    # 计算每个副本的总表达量
    expr_sums <- rowSums(GetAssayData(sc_data_converted, slot = "counts")[dup_idx, ])

    # 保留最高表达量的副本
    keep_this <- dup_idx[which.max(expr_sums)]
    discard <- setdiff(dup_idx, keep_this)

    keep_indices[discard] <- FALSE
  }

  sc_data_converted <- sc_data_converted[keep_indices, ]
  new_genes <- new_genes[keep_indices]  # 同步更新
}


saveRDS(sc_data_converted, "rat_converted_sc_data.rds")


