#Set working directory
setwd("/data/projects/mqtl/")

#Load packages
library(data.table)

#Load RNT transform values
rnt <- fread("rnt_betas.txt")

#Transform BED file into data.frame with columns as samples, and rows as CpGs
cpgs <- rnt$ID  # get vector with CpG's ID
rnt <- rnt[,-c(1:4)] # remove chr, start, end and ID columns
rnt <- as.data.frame(rnt) # convert to dataframe
rownames(rnt) <- cpgs # set rownames as CpG's ID

#Load covariates files 
cov <- read.table("covariates/covariates_tmp.txt", header = T)

#Order covariate file as rnt file
samples <- colnames(rnt)
index<- match(samples, colnames(cov))
cov <- cov[, index] 
table(colnames(cov) == colnames(rnt))

#Create a matrix with methylation values and known confounders 
colnames_cov <- rownames(cov) #get vector with covariates name
cov_t <- as.data.frame(t(cov)) # transpose datatable
cov_t <- apply(cov_t, 2, as.numeric)
rownames(cov_t) <- colnames(cov)

#RNT file
rnt_t <- t(rnt) #transpose datatable 
dim(rnt_t)
rnt_t <- as.data.frame(rnt_t)
# columns as CpGs, rows as samples

#Merge covariates and RNT data frames  
merged_df <- merge(x = rnt_t, y = cov_t, by ="row.names")
rownames(merged_df) <- merged_df$Row.names # assign as rownames sample FID_IID

#Perform mulitple linear model (methylaion ~ sex + PC1 + PC2 + PC3 + PC4 + PC5 + age)
#Get dataframe with only methylation values
cpgs <- rownames(rnt) # get vector with CpG's ID
subset_cpgs <- merged_df[, cpgs] # subset dataframe by CpGs
rownames(subset_cpgs) <- rownames(merged_df) # set sample's name as rownames

#Apply as.numeric to all covariates
merged_df[colnames_cov] <- sapply(merged_df[colnames_cov],as.numeric)

#Perform the multiple linear model
fit <- lm(as.matrix(subset_cpgs) ~ sex + PC1_geno + PC2_geno + PC3_geno + PC4_geno + PC5_geno + age, data = merged_df)

#Get residuals from multiple linear model
residuals <- fit$residuals

#Perform PCA on the residuals 
pca <- prcomp((na.omit(residuals)))
colnames(pca$x) <- paste("m", colnames(pca$x), sep="")

#Get variance explained per PC 
prop<-pca$sdev^2/sum(pca$sdev^2)*100
cumprop<-cumsum(pca$sdev^2)/sum(pca$sdev^2)*100
pcn<-paste("mPC", 1:length(prop), sep="")
pc<-data.frame(pcn,prop,cumprop)

#Write text file with the mPCs (select a maximum of 20 mPCs)
pca <- pca$x[,c(1:20)]
write.table(x = pca6, file = "/data/projects/mqtl/residualized_mPCs.txt", col.names = T, row.names = T, sep="\t", quote = F)
