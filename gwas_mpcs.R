#Set working directory
setwd("/data/projects/mqtl/")

#Load packages
library(data.table)
library(dplyr)

#Step 1: Load fam file
fam <- fread("/data/projects/geno/all_chr_imputed_rsq09_maf001_hwe005_chr1_22.fam", col.names = c("FID","IID","fID","mID","sex","pheno"))

#Step 2: Load mPCs
mpcs <- fread("residualized_mPCs.txt")
colnames(mpcs)[1] <- "IID"

#Step 3: Merge fam file with mPCs
fam_modified <- merge(fam, mpcs, by= "IID", sort = F)
fam_modified <- fam_modified[,c(2,1,3:26)]

#Step 4: Creating one fam file per mPC:
for (i in colnames(mpcs)[-1]) {
  print(paste0("working on ", i))
  data <- fam_modified
  data$pheno <- subset(data,select = i) # assign as phenotype column an mPC
  data <- data[, c(1:6)] # select only PLINK fam file columns
  fwrite(x = data, file = paste0("gwas_pcs/all_chr_imputed_rsq09_maf001_hwe005_chr1_22_", i, ".fam"), # write fam file with phenotype being an mPC
         col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
}

#Step 5: Perform GWAS
for (i in colnames(mpcs)[-1]) {
  file <- paste("gwas_pcs/all_chr_imputed_rsq09_maf001_hwe005_chr1_22_", i, ".fam", sep="") # create string with name of the modified fam file
  print(file)
  output_file <- paste("gwas_pcs/gwas_", i, sep="") # create string with name of the output file
  print(output_file)
  gwas_command <- paste(paste(paste("plink --bfile gwas_pcs/all_chr_imputed_rsq09_maf001_hwe005_chr1_22 --fam ", file, sep="")," --assoc --out ",sep=""), output_file, sep="") # write string with GWAS command
  print(gwas_command)
  system(gwas_command) # perform GWAS using the string with the command
}
