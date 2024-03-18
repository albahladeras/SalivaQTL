#Set working directory
setwd("/data/projects/eqtl/")

##Load packages
library(data.table)
library(dplyr)
library(edgeR)
library(stringr)

gene_counts <- fread("/data/preprocess_exo/06_count_matrix/gene_count_matrix.csv")
sample_more5M<- fread("sample_more5M.txt")

#Subset good quality samples (more than 5M reads)

gene_id <- gene_counts_filt[,1]
gene_counts_filt <- as.data.frame(gene_counts_filt)
gene_counts_filt_samples <- gene_counts_filt[,-1][,colnames(gene_counts_filt)[-1] %in% sample_more5M$Sample]
gene_counts_filt_samples <- cbind(gene_id, gene_counts_filt_samples)
rm(gene_counts_filt, gene_id)

#How many genes have more than 6 reads in more than 20% of te samples?
keep_genes <- rowSums(gene_counts_filt_samples[,-1] >= 6) > (0.2 * ncol(gene_counts_filt_samples[,-1]))
number_genes_kept <- sum(keep_genes) 
gene_counts_filt_samples2 <- gene_counts_filt_samples[keep_genes,]

#TMM normalize
gene_counts_tmp <- gene_counts_filt_samples2[,-1]
gene_counts_tmm <- DGEList(counts = gene_counts_tmp, genes = gene_counts_filt_samples2[, "gene_id"])
gene_counts_tmm <- calcNormFactors(gene_counts_tmm, method = "TMM")
gene_counts_tmm <- cpm(gene_counts_tmm)
dim(gene_counts_tmm)
gene_counts_tmm <- as.data.frame(gene_counts_tmm)

#Add gene ids

gene_counts_tmm <- cbind(gene_counts_filt_samples2$gene_id, gene_counts_tmm)

#Format colnames
colnames(gene_counts_tmm)[1] <- "gene_id"

write.table(gene_counts_tmm, file = "genes_filt_tmm.txt", col.names = T,
            row.names = F, sep = "\t", quote = F)

