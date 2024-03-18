#Set working directory
setwd("/data/projects/eqtl")

#Load packages
library(data.table)

expr <- fread("/genes_filt_tmm.txt") #samples in columns, genes in rows

#Inverse normal transformation
data_trans <- qnorm(t(apply(expr[,-1], 1, rank, ties.method = "average"))/ (ncol(expr[,-1])+1))
data_trans <- as.data.frame(data_trans)
rownames(data_trans) <- expr$genes

write.table(data_trans, file = "expr_matrix_filtgenes_rnt_tmm.txt", 
            col.names = T, row.names = T, quote = F, sep = "\t")
