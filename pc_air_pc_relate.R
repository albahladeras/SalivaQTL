
#Set working directory
setwd("/data/projects/mqtl/pc_air_pc_relate")

#Loadp ackages
library(SNPRelate)
library(GWASTools)
library(GENESIS)
library(ggplot2)

# The SNPRelate package provides the snpgdsBED2GDS function to convert binary 
# PLINK files into a GDS file.

gdsfile <- snpgdsBED2GDS(bed.fn = "../imputed_rsq09_maf001_hwe005_chr1_22_common_methyl.bed", 
                         bim.fn = "../imputed_rsq09_maf001_hwe005_chr1_22_common_methyl.bim",
                         fam.fn = "../imputed_rsq09_maf001_hwe005_chr1_22_common_methyl.fam",
                         out.gdsfn = "genotype.gds")

# read in GDS data
gds <- snpgdsOpen(gdsfile)
snpset <- snpgdsLDpruning(gds, method="corr", slide.max.bp=10e6, 
                          ld.threshold=sqrt(0.1), verbose=FALSE)
pruned <- unlist(snpset, use.names=FALSE)

# Pairwise Measures of Ancestry Divergence
# create matrix of KING estimates

king <- snpgdsIBDKING(gds)

kingMat <- king$kinship
dimnames(kingMat) <- list(king$sample.id, king$sample.id)
snpgdsClose(gds)

kinship <- snpgdsIBDSelection(king)

# Running PC-AiR
saliva_geno <- GdsGenotypeReader(filename = "genotype.gds")
# create a GenotypeData class object
saliva_genoData <- GenotypeData(saliva_geno)
saliva_genoData

# run PC-AiR on pruned SNPs
mypcair <- pcair(saliva_genoData, kinobj = kingMat, divobj = kingMat,
                 snp.include = pruned)

colnames(mypcair$vectors) <- paste0("PC", seq(1:32))
write.table(x = mypcair$vectors, file = "genetic_pcs.txt", col.names = T, row.names = T, 
            quote = F, sep = "\t")

# Running PC-Relate 
saliva_genoData <- GenotypeBlockIterator(saliva_genoData, snpInclude=pruned)
mypcrelate <- pcrelate(saliva_genoData, pcs = mypcair$vectors[,1:5], 
                       training.set = mypcair$unrels,
                       BPPARAM = BiocParallel::SerialParam())

iids <- as.character(getScanID(saliva_genoData))
GRM <- pcrelateToMatrix(mypcrelate)
GRM_df <- as.matrix(GRM)

write.table(GRM_df, file =  "genetic_relatedness_matrix.txt", col.names = T, row.names = T, 
            quote = F, sep = "\t")


