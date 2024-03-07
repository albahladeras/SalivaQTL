## Quality control of methylation data

#Set working directory
setwd("/data/projects/EPIC")

#Load packages
library(vctrs)
library(minfi)
library(ChAMP)
library(ChAMPdata)
library(impute)
library(sva)
library(data.table)
library(PACEanalysis)

#Filter samples
myLoad <- champ.load(directory = "/data/saliva/meth", method="minfi",
                     methValue="B",
                     autoimpute=FALSE,
                     filterDetP=T,
                     ProbeCutoff=0.05,
                     SampleCutoff=0.05,
                     detPcut=0.01,
                     filterBeads=T,
                     beadCutoff=0.05,
                     filterNoCG=TRUE,
                     filterSNPs=TRUE,
                     population="EUR",
                     filterMultiHit=TRUE,
                     filterXY=TRUE,
                     force=FALSE,
                     arraytype="EPIC")

saveRDS(myLoad, "myLoad.RDS")

#Functional and Noob Normalization 
myNorm2.temp <- preprocessFunnorm(myLoad$rgSet)
beta_temp <- getBeta(myNorm2.temp)

#Removing failed probes
myNorm2 <- beta_temp[rownames(beta_temp) %in% rownames(myLoad$beta),]
saveRDS(myNorm2, file = "myNorm2.RDS")

#Removing CpGs from the Mcartney et al Genome Biology
cpgs1 <- fread("mmc1.txt.gz")
cpgs2 <- fread("mmc2.txt.gz", header = F)
cpgs3 <- fread("mmc3.txt.gz", header = F)

#Remove CpGs with SNPs with a MAF in european samples higher than 0.05
cpgs_eur <- cpgs1[cpgs1$EUR_AF > 0.05,]

myNorm2 <- myNorm2[!rownames(myNorm2) %in% cpgs_eur$IlmnID,]
rm(cpgs1, cpgs_eur)
myNorm2 <- myNorm2[!rownames(myNorm2) %in% cpgs2$V1,]
myNorm2 <- myNorm2[!rownames(myNorm2) %in% cpgs3$V1,]

#The myNorm2 beta matrix contains beta values from failed probes with a detPval > 0.01
#but they were kept because the samplecutoff. We need NA values in this positions

#1. Subset 

beta_subset <- myLoad$beta[rownames(myLoad$beta) %in% rownames(myNorm2),]

#2. Order the beta matrix according the cpgs order in myNorm2
beta_subset <- beta_subset[order(match(rownames(beta_subset), rownames(myNorm2))), ]

#3. Insert Na
myNorm2[is.na(beta_subset)] <- NA 
table(is.na(myNorm2)) 

saveRDS(myNorm2, file = "myNorm_withNA.RDS")

#Impute Na values)
tmp <- impute.knn(myNorm2,k=5)$data

#BMIQ Normalization (within array normalization)
myNorm3 <- champ.norm(beta=tmp,
                      rgSet=NULL,
                      mset=NULL,
                      resultsDir="CHAMP_Normalization",
                      method="BMIQ",
                      plotBMIQ=TRUE,
                      arraytype="EPIC",
                      cores=3)

saveRDS(myNorm3, file = "myNorm3.RDS")

#Singular value decomposition method (SVD) (execute in terminal)
champ.SVD(beta=myNorm3,pd=myLoad$pd)

## Remove batch effects with combat and treat outliers. PACE Analysis R package

## Calculate the variance of each probe and remove any with a variance of 0 prior to Combat
mval <- apply(myNorm3, 2, function(x) log2((x)/(1-x)))
vars = as.matrix(rowVars(mval))

## Replace all probes with no variance with NA and remove them from the FunNorm set
vars[vars == 0] = NA
vars = na.omit(vars)
intersectingCpGs = intersect(rownames(vars), rownames(mval))
cat(nrow(mval) - length(intersectingCpGs), "probes had no variance and were removed","\n")
fn.sub = myNorm3[intersectingCpGs,]
mval = mval[intersectingCpGs,]

expit2 = function(x) 2^x/(1+2^x)
modcombat <- model.matrix(~1, data=myLoad$pd)
combat.adj <- ComBat(dat=mval, batch=myLoad$pd$Slide, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
betafinal = expit2(combat.adj)

destf <- "/data/projects/EPIC/"
analysisdate <- "28_02_2023"
outlier <- outlierprocess(processedBetas = betafinal,
                          quantilemethod = "EmpiricalBeta",
                          trimming = FALSE,
                          pct = 0.005,
                          destinationfolder = destf,
                          cohort = NULL,
                          analysisdate = analysisdate)


saveRDS(object = outlier,"betafinal_outlierremoved.RDS")
