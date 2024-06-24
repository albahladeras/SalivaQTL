#Set working directory
setwd("/data/projects/mqtl/smr")

#Load packages
library(stringr)
library(data.table)

#Load files
myeqtl.esi <- fread("../imputed_rsq09_maf001_hwe005_chr1_22_common_methyl.bim")
freq <- fread("imputed_rsq09_maf001_hwe005_chr1_22_common_methyl.frq", header = T)
eqtl <- fread("mqtls_nominal.txt")

myeqtl.esi <- cbind(myeqtl.esi, freq[, 5])
myeqtl.esi$V2 <- paste(str_split_i(myeqtl.esi$V2, pattern = ":", i = 1), str_split_i(myeqtl.esi$V2, pattern = ":", i = 2), sep = ":")

myeqtl.esi2 <- merge(eqtl,myeqtl.esi, by = "V2", sort = FALSE)
myeqtl.esi2 <- myeqtl.esi2[,c(6, 1, 7:11)]

myeqtl.esi3<- myeqtl.esi2[!duplicated(myeqtl.esi2$V2),]

write.table(x = myeqtl.esi3, file = "mymqtl.esi", col.names = F, row.names = F, quote = F, sep = "\t")
