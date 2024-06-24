#Write files in FastQTL format to perform SMR

#Load packages
setwd("/data/projects/mqtl/smr")
library(data.table)
library(stringr)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

#Load files
res1e03 <- fread("MQTL_LMM_1e03.txt")
cpg_info <- fread("MQTL_LMM.cis_gene_table.txt.gz")

colnames(cpg_info) <- paste(colnames(cpg_info), "cpg",sep = "_")

df_1e03 <- merge(res1e03, cpg_info, by.x = "gene", by.y = "gene_cpg")
df_1e03 <- df_1e03[,-c(9, 11:16)]
colnames(df_1e03) <- c("cpg", "chr", "snp_pos", "ref", "alt", "b", "se", "p_nominal", 
                       "cpg_pos")
df_1e03$distance <- df_1e03$snp_pos - df_1e03$cpg_pos
df_1e03$snp <- paste0(df_1e03$chr, ":",df_1e03$snp_pos)

fast_qtl <- df_1e03[, c(1, 11, 10,8, 6)]

write.table(fast_qtl, file = "/data/projects/mqtl/smr/fastqtlnomi.txt", 
            col.names = F, row.names = F, quote = F, sep = "\t")

#Create the epi file for SMR
# Get genes and strand
annot <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
annot <- cbind(annot$strand, annot$Name, annot$UCSC_RefGene_Name, annot$pos, annot$chr) 

df <- merge(fast_qtl[,1:2], annot, by.x = "cpg", by.y = "V2", sort = FALSE)
df$genetic_distance <- 0

df <- df[,c(6,1,7,5,4,3)]
df$V5 <- gsub(pattern = "chr", replacement = "", x = df$V5)
df <- df[!duplicated(df$cpg),]
df$V3[df$V3 == ""] <- "."

write.table(df, "/data/projects/mqtl/smr/mymqtl.epi", col.names = F, row.names = F, 
            quote = F, sep = "\t")


