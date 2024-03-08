
#Set working directory
setwd("/data/projects/mqtl/")

#Load packages
library(ggplot2)
library(data.table)

methyl <- fread("methylome.bed")

#Inverse normal transformation
data_trans <- qnorm(t(apply(methyl[,-c(1:4)], 1, rank, ties.method = "average"))/ (ncol(methyl[,-c(1:4)])+1))
data_trans <- cbind(methyl[,c(1,2,3,4)], data_trans)
write.table(x = data_trans, file = "rnt_betas.txt", 
            col.names = T, row.names = F, sep = "\t", quote = F)
