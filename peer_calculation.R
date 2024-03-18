#PEER code from GTEx pipeline https://raw.githubusercontent.com/broadinstitute/gtex-pipeline/master/qtl/src/run_PEER.R

# Being run with GTEX image

setwd("mnt/eqtl/")

library(peer)

# Load covariates data
cov <- read.table("covariates_expr_tmp.txt", header = T)
head(cov)
cov_t <- t(as.matrix(cov))
table(sapply(cov_t, class)) 
cov_t <- apply(cov_t,2, as.numeric)
table(sapply(cov_t, class))

# Load expression data
expr <- read.table("expr_matrix_filtgenes_rnt_tmm.txt", header = T)
#order expr as cov
samples <- colnames(cov)
index <- match(samples, colnames(expr))
expr <- expr[,index]
M <- t(as.matrix(expr))
colnames(M) <- rownames(expr)


#Number of hidden confounders to estimate"
#according to gtex pipeline https://github.com/broadinstitute/gtex-pipeline/tree/master/qtl
#set number of PEER factors
#to 15 factors for N < 150
#to 30 factors for 150 ≤ N < 250 (choose this for 200 samples)
#to 45 factors for 250 ≤ N < 350
#to 60 factors for N ≥ 350


getNumPeer <- function(ss) {
  if (ss<150) return (min(15, ceiling(ss / 5)))
  else if (ss >=150 && ss < 250) return(30)
  else return(35)
}
# n = getNumPeer(nrow(M)) #number of peer factors
n = 30
alphaprior_a=0.001
alphaprior_b=0.01
epsprior_a=0.1
epsprior_b=10
max_iter=1000

# Set model object

model <- PEER() # create object model

PEER_setPhenoMean(model, M) # set the observed data
dim(PEER_getPhenoMean(model)) #

PEER_setNk(model, n) # set number of k factors
PEER_getNk(model)  #

# Set default values
PEER_setPriorAlpha(model, alphaprior_a,alphaprior_b)
PEER_setPriorEps(model, epsprior_a, epsprior_b)
PEER_setNmax_iterations(model, max_iter)

PEER_setCovariates(model, as.matrix(cov_t))
dim(PEER_getCovariates(model))

# Perform the inference
PEER_update(model)

# Get PEER factors
factors <- PEER_getX(model) 
dim(factors)
#247 300
saveRDS(factors, "peer/factors.RDS")

# Get weights
weights <- PEER_getW(model)
dim(weights)

saveRDS(object = weights, file = "peer/weights.RDS")

# Get precision (inverse variance) of the weights
precision <- PEER_getAlpha(model)
dim(precision)

saveRDS(object = precision, file = "peer/precision.RDS")

# Get residuals 
residuals <- PEER_getResiduals(model)
residuals_t <- t(PEER_getResiduals(model))  # genes x samples
saveRDS(object = residuals, file = "peer/residuals.RDS")
saveRDS(object = residuals_t, file = "peer/residuals.RDS")

# add relevant row/column names
c <- paste0("InferredCov",1:30)
rownames(factors) <- rownames(M)
colnames(factors) <- c(colnames(cov_t), c)
rownames(precision) <- c(colnames(cov_t), c)
colnames(precision) <- "Alpha"
precision <- as.data.frame(precision)
precision$Relevance <- 1.0 / precision$Alpha
rownames(residuals_t) <- colnames(M)
colnames(residuals_t) <- rownames(M)

#write results
factors_df<-as.data.frame(factors)
factors_df$ID <- rownames(factors_df)
write.table(factors_df, file = "peer/PEERcovars.txt" , sep = "\t", row.names = TRUE, quote=FALSE, 
            col.names = T)
write.table(precision, file = "peer/PEERprecision.txt", sep = "\t", row.names = TRUE, quote=FALSE, 
            col.names = T)
write.table(residuals_t, file = "peer/PEERresiduals_t.txt", sep = "\t", row.names = TRUE, quote=FALSE, 
            col.names = T)

