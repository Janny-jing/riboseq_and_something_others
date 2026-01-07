# devtools::install_github("mariodosreis/tai")
# require("tAI")
rm(list = ls())
library(tAI)
library(ggplot2)
library(dplyr)
library(coRdon)
setwd("G:\\temp\\06.trna\\07.zhang\\tai_2")

#code{sking} indicates the superkingdom, with 0 indicating Eukaryota, and 1 Prokaryota.

calculate_tAI <- function(trna_file_path, m_file_path, ffn_file_path) {
  # 从文件读取 tRNA 数据
  eco.trna <- scan(trna_file_path)
  eco.ws <-  get.ws(tRNA=eco.trna, sking=0)
  
  # 从文件读取密码子频率矩阵
  eco.m <- matrix(scan(m_file_path), ncol=61, byrow=TRUE)
  eco.m <- eco.m[,-33]  # 删除甲硫氨酸
  get.tai_processs <- function(x,w) {
    n = apply(x,1,'*',w)    #multiply each row of by the weights
    n = t(n)                #transpose
    return(n)
  }
  # 计算 tAI 值
  eco.tai <- get.tai_processs(eco.m, eco.ws)
  
  # 从 DNA 文件读取数据
  dna_file <- readSet(file =ffn_file_path) %>% 
    codonTable(.)
  
  # 存储结果到数据框
  out_data <- data.frame(
    gene_name = dna_file@ID,
    tAI = eco.tai
  )
  
  return(out_data)
}

# 使用函数并返回结果
result_data <- calculate_tAI(trna_file_path = "G:\\temp\\06.trna\\07.zhang\\tai_2\\ctrinfected_isodecoder.txt",
                             m_file_path = "G:\\temp\\06.trna\\07.zhang\\tai_2\\virus.m",
                             ffn_file_path = "G:\\temp\\06.trna\\07.zhang\\tai_2\\virus.fa")
write.csv(result_data, "./ctrinfected_isodecoder_tAI.csv", row.names = F)


# plot(eco.tai, ecolik12$w$Nc, xlab="tAI", ylab="Nc")
