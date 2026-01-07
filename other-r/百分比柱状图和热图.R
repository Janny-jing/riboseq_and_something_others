rm(list=ls())
setwd("G:\\temp\\20.codons_propotion")
library(plyr)
data <- read.table("ctl_polyic_wsn_sas.txt",header = TRUE)
# 
# # 将宽数据转换为绘制分组柱状图的长数据
# data2 <- data %>% pivot_longer(cols=c(NO:NO2),
#                                names_to = 'species',
#                                values_to = 'con')

data1 <- ddply(data,'sample',transform,percent_con=propotion/sum(propotion)*100)

data1$sample <- factor(data1$sample,levels=c("ctl","pic","wsn","mock","infect"))

ggplot(data1,aes(x=sample,y=percent_con,fill=base))+
  geom_bar(stat="identity",colour="black")+
  guides(fill=guide_legend(reverse = TRUE))+
  scale_fill_brewer(palette="BuPu")

rm(list=ls())
setwd("G:\\temp\\20.codons_propotion")
library(plyr)
library(ggplot2)
library(reshape2)

data_diff <- data.frame(
  Sample = c("pic_ctl", "wsn_ctl", "infected"),
  A = c(0.0202, -0.1228, -0.0174),
  T = c(-0.3254, -0.1813, 0.0068),
  C = c(0.5131, 0.36, 0.0356),
  G = c(-0.2079, -0.0559, -0.025)
)
data_long <- melt(data_diff, id.vars = "Sample")

ggplot(data_long, aes(x = variable, y = Sample, fill = value)) +
  geom_tile() +
  geom_text(aes(label = round(value, 2)), color = "black", size = 3) +
  scale_fill_gradient2(low = "#132B43", mid = "white", high = "#56B1F7", midpoint = 0) +
  labs(title = "Difference in Proportion of ATCG from Controls", x = "Base", y = "Sample") +
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


