### MAtrixEQTL

#Provided by C. Lesseur

#Modified by Alba Hernangomez-Laderas

## Set working directory
setwd('/data/projects/eqtm/')

## Load packages
library(MatrixEQTL)

#data & output files
met_file_name="methylation.txt"
exp_file_name="expression.txt"
covariates_file_name ="covariates.txt"
output_file_name_cis = "eqtms"


#genomic position files
met_location_file_name<-read.table("met_location.txt",header = TRUE, stringsAsFactors = FALSE)
exp_location_file_name<-read.table("exp_location.txt",header = TRUE, stringsAsFactors = FALSE)


#general parameter MatrixQTL
#model
useModel = modelLINEAR;
# Only associations significant at this level will be saved
pvOutputThreshold_cis = 0.99;
pvOutputThreshold_trans = 0;
errorCovariance = as.matrix(read.table("genetic_relatedness_matrix.txt"));
# Distance for local CpG-Gene pairs
cisDist=1e6

##Load methylation data
snps= SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 10,000 rows
snps$LoadFile(met_file_name);

## Load gene expression data
gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 10,000 rows
gene$LoadFile(exp_file_name);

#Load_covars
cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name);
}


## Run the analysis
eqtm<-Matrix_eQTL_main(
  snps = snps, 
  gene = gene, 
  cvrt = cvrt,
  output_file_name= NULL,
  pvOutputThreshold   = pvOutputThreshold_trans,
  useModel = useModel, 
  errorCovariance = errorCovariance, 
  verbose = TRUE, 
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = met_location_file_name, 
  genepos = exp_location_file_name,
  cisDist = cisDist,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);

unlink(output_file_name_cis);

##results#####
eqtmcis=(eqtm$cis$eqtls)

#Add standard error information
eqtmcis$beta_se = eqtmcis$beta / eqtmcis$statistic

write.table(eqtmcis, file = "eQTMs_gene_level.txt", col.names = T, row.names = F,
            quote = F, sep = "\t")
