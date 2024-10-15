## Create the BED file for mQTL analysis 

#Set working directory
setwd("/data/projects/mqtl")

#Load packages
library(stringr)
#Load QC-ed methylation data
outlier <- readRDS("../EPIC/betafinal_outlierremoved.RDS")

#Annotate CpGs

#1. Get annotation
annot <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

#2. Merge annotation and beta data.frames by CpG IDs
annot_beta <- merge(x=annot[,c(1,2,2,4)], y=outlier, by="row.names")
rm(annot)
rownames(annot_beta) <- annot_beta$Row.names #Make CpG names' as row.names of the data.frame 
annot_beta <- annot_beta[,-1]  #Remove column with CpG names'

#3. Change "chr1" by "1" 
annot_beta$chr <- str_remove(annot_beta$chr, "chr")

#4. Change column names
colnames(annot_beta)[1:4] <- c("#Chr", "start", "end", "ID")

#5. Define the cis-window (end = start + 1)
annot_beta$end <- annot_beta$start + 1

#6: start and end columns must be integer
class(annot_beta$start)
# [1] "integer"
class(annot_beta$end)
# [1] "numeric"
annot_beta$end <- as.integer(annot_beta$end)
class(annot_beta$end)

#7: Obtain the final dataframe inside a text file following a BED structure format
write.table(x = annot_beta, file = "methylome_BED_with_all_samples.txt", quote = FALSE, 
            sep = "\t", row.names = FALSE, col.names = TRUE)
