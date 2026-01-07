if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("gridExtra"))

library(msigdbr)
library(msigdbdf)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(gridExtra)
library(openxlsx)
library(purrr)
library(readr)
library(stringr)
library(tools)

# ---- Load Hallmark Database ----
#retrieve Hallmark Gene Sets for Mouse - this shouldn't need editing, just run as is
hallmark_mouse <- msigdbr(species = "Mus musculus", category = "H") %>%
  dplyr::select(gs_name, gene_symbol)

#convert into a list format required for enrichment analysis
hallmark_list <- split(hallmark_mouse$gene_symbol, hallmark_mouse$gs_name)

#create R vectors for each gene set
for (term in names(hallmark_list)) {
  assign(term, hallmark_list[[term]])
  cat(paste0(term, " <- c(\n  \"", paste(hallmark_list[[term]], 
                                         collapse = "\", \""), "\"\n)\n\n"))}


#load working directory
setwd("insert pathname")



# ---- Upregulated Pathway Analysis ----
Upregulated_DEGs <- read.csv("insert data file name.csv", header = TRUE)
Upregulated_DEGs <- Upregulated_DEGs %>% filter(log2FoldChange > 0.6, padj < 0.05) #change this for log2FC of 0.6 or 1
up_genes <- Upregulated_DEGs$genes #column name is genes; log2FoldChange and padj are column names too

# Run enrichment analysis (Hallmark)
up_enriched <- enricher(
  gene = up_genes,
  TERM2GENE = hallmark_mouse,
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.25)
up_enriched@result$Description <- gsub("^HALLMARK_", "", up_enriched@result$Description)

# Export results to a CSV file
write.csv(up_enriched@result, file = "insert output file name.csv", row.names = FALSE)

hallmark_summary_up <- as.data.frame(up_enriched)
up_enriched_filtered <- up_enriched@result %>%
  filter(p.adjust < 0.05)

# Format the Description: remove underscores and capitalize each word
up_enriched_filtered$Formatted_Description <- gsub("_", " ", up_enriched_filtered$Description)
up_enriched_filtered$Formatted_Description <- toTitleCase(tolower(up_enriched_filtered$Formatted_Description))

# Get top 20 pathways by adjusted p-value
top_up_enriched <- up_enriched_filtered %>%
  arrange(p.adjust) %>%
  head(20)

# Calculate upper limit with a small buffer (e.g., 10% of the max)
x_limit <- max(top_up_enriched$FoldEnrichment) * 1.1

# Custom ggplot barplot (currently for top 20 pathways)
ggplot(top_up_enriched) +  
  geom_bar(aes(x = FoldEnrichment, y = reorder(Formatted_Description, FoldEnrichment)), 
           stat = "identity", fill = "firebrick", width = 0.8) +  
  geom_text(aes(x = FoldEnrichment, 
                y = reorder(Formatted_Description, FoldEnrichment), 
                label = Count), 
            hjust = -0.25, size = 5, color = "black") +  
  labs(
    title = "Hallmark Enrichment Analysis:\nGenes Upregulated in ???",
    x = "Fold Enrichment",
    y = "") +  
  theme_minimal() +  
  theme(
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.text.x = element_text(size = 15),
    plot.title = element_text(size = 18, face = "bold"),
    legend.position = "none",
    axis.title = element_text(size = 15)) +  
  scale_x_continuous(name = "Fold Enrichment", limits = c(0, x_limit))




# ---- Downregulated Pathway Analysis ----
Downregulated_DEGs <- read.csv("insert data file name.csv", header = TRUE)
Downregulated_DEGs <- Downregulated_DEGs %>% filter(log2FoldChange < -0, padj < 0.05) #change this to be log2fc of -0.6 or -1
down_genes <- Downregulated_DEGs$genes

# Run enrichment analysis (Hallmark)
down_enriched <- enricher(
  gene = down_genes, 
  TERM2GENE = hallmark_mouse, 
  pvalueCutoff = 0.05, 
  qvalueCutoff = 0.25)
down_enriched@result$Description <- gsub("^HALLMARK_", "", down_enriched@result$Description)

write.csv(down_enriched@result, file = "insert output file name.csv", row.names = FALSE)

hallmark_summary_down <- as.data.frame(down_enriched)
down_enriched_filtered <- down_enriched@result %>%
  filter(p.adjust < 0.05)

# Format the Description: remove underscores and capitalize each word
down_enriched_filtered$Formatted_Description <- gsub("_", " ", down_enriched_filtered$Description)
down_enriched_filtered$Formatted_Description <- toTitleCase(tolower(down_enriched_filtered$Formatted_Description))

# Get top 20 pathways by adjusted p-value
top_down_enriched <- down_enriched_filtered %>%
  arrange(p.adjust) %>%
  head(20)

# Plotting for downregulated genes
ggplot(top_down_enriched) + 
  # Bar graph for fold enrichment with solid color and horizontal bars
  geom_bar(aes(x = FoldEnrichment, y = reorder(Formatted_Description, FoldEnrichment)), 
           stat = "identity", fill = "mediumblue", show.legend = FALSE, width = 0.8) + 
  # Gene counts text beside the bars
  geom_text(aes(x = FoldEnrichment, 
                y = reorder(Formatted_Description, FoldEnrichment), 
                label = Count), 
            hjust = -0.25, size = 5, color = "black", show.legend = TRUE) + 
  labs(
    title = "Hallmark Enrichment Analysis:\nGenes Downregulated in ???",
    x = "Fold Enrichment",
    y = "") +  
  theme_minimal() +  
  theme(
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.text.x = element_text(size = 15),  # Change size of fold enrichment values on x-axis
    plot.title = element_text(size = 18, face = "bold"),  # Larger title
    legend.position = "none",  # No legend
    axis.title = element_text(size = 15)) +  
  scale_x_continuous(
    name = "Fold Enrichment",   
    limits = c(0, max(top_down_enriched$FoldEnrichment) + 0.5))


