rm(list=ls())

source("http://bioconductor.org/biocLite.R")
biocLite("BiocUpgrade")

library(devtools)
library('Aclust') 
library("gee")
library(stringr)

##########################################################################

## GEE function

GEE.clusters <- function(betas, clusters.list, exposure, id, covariates = NULL,  working.cor = "ex", minimum.cluster.size = 2, result.file.name = NULL){
  ## the ids are rownames of the betas, and are ordered by the same ordering of the covariates matrix.
  ## function that gets a matrix of beta values, a list of clusters, exposure variable, and covariates
  ## and performs GEE. Returns a matrix with the results of the GEE model (only the exposure effect and p-value)
  ## for each cluster. If the cluster has a single probe?
  ## The results are ordered according to p-value
  ## Note: methylation is treated here as outcome. 
  ## if file is given, print results to file. 
  
  require(geepack)
  
  n.sites <- unlist(lapply(clusters.list, function(x) return(length(x))))
  inds.rm <- which(n.sites < minimum.cluster.size)
  
  if (length(inds.rm) > 0) clusters.list <- clusters.list[-inds.rm]
  
  
  n.mod <- length(clusters.list)
  
  effect <- rep(0, n.mod)
  se <- rep(0, n.mod)
  pvals <- rep(0, n.mod)
  sites <- rep("", n.mod)
  n.sites <- rep(0, n.mod)
  
  if (!is.null(covariates)){
    
    if (is.null(dim(covariates))) covariates <- as.matrix(covariates)
    
    if (is.null(colnames(covariates)))  colnames(covariates) <- rep("", ncol(covariates))
    for (i in 1:ncol(covariates)){
      if (colnames(covariates)[i] == "") colnames(covariates)[i] <- paste("covariate", i, sep = '_')
      
    }
    
    model.expression <- "model <- geeglm(Beta ~ exposure +"
    for (j in 1:ncol(covariates)){
      model.expression <- paste(model.expression, colnames(covariates)[j], "+")  }
    
    model.expression <- paste(model.expression, "as.factor(probeID), id = as.numeric(id), data = temp.long, corstr = working.cor)")
    
    model.expr.1.site <- paste("model <- geeglm(clus.betas[ind.comp] ~ exposure[ind.comp] + ")  
    for (j in 1:ncol(covariates)){
      if (j < ncol(covariates)) model.expr.1.site <- paste(model.expr.1.site, colnames(covariates)[j], "+")
      else    model.expr.1.site <- paste(model.expr.1.site, colnames(covariates)[j])}
    
    model.expr.1.site <- paste(model.expr.1.site, ", data = covariates[ind.comp,], id = as.numeric(id)[ind.comp])")
    covariates <- as.data.frame(covariates)
    
    
    for (i in 1:n.mod){
      clus.probes <- clusters.list[[i]]
      clus.betas <- betas[clus.probes,]
      
      if (length(clus.probes) == 1){
        
        ind.comp <- which(complete.cases(cbind(clus.betas, exposure, covariates, id)))
        eval(parse(text = model.expr.1.site))
        
        effect[i] <- summary(model)[[6]][2,1]
        se[i] <- summary(model)[[6]][2,2]
        pvals[i] <- summary(model)[[6]][2,4]
        n.sites[i] <- 1
        sites[i] <- clus.probes
        
        
      } else{
        
        ind.comp <- which(complete.cases(cbind(t(clus.betas), exposure, covariates, id)))
        
        temp.betas <- clus.betas[,ind.comp]
        
        temp.covars <- as.data.frame(cbind(exposure[ind.comp], covariates[ind.comp,], id[ind.comp]))
        colnames(temp.covars) <- c("exposure", colnames(covariates), "id")
        
        temp.long <- organize.island.repeated(temp.betas, temp.covars)
        temp.long <- temp.long[complete.cases(temp.long),]
        temp.long <- temp.long[order(temp.long$id),]
        temp.long$probeID <- as.factor(temp.long$probeID)
        #temp.long$exposure <- as.numeric(as.character(temp.long$exposure))
        temp.long$exposure <- as.factor(temp.long$exposure)
        temp.long$Beta <- as.numeric(as.character(temp.long$Beta))
        
        eval(parse(text = model.expression))
        
        
        effect[i] <- summary(model)[[6]][2,1]
        se[i] <- summary(model)[[6]][2,2]
        pvals[i] <- summary(model)[[6]][2,4]
        n.sites[i] <- length(clus.probes)
        for (c in 1:length(clus.probes)) sites[i] <- paste(sites[i], clus.probes[c], sep = ";")
        substr(sites[i], 1, 1) <- ""
        
      }
      
    }
    
  } else{
    
    model.expression <- "model <- geeglm(Beta ~ exposure +"
    model.expression <- paste(model.expression, "as.factor(probeID), id = as.numeric(id), data = temp.long, corstr = working.cor)")
    
    model.expr.1.site <- paste("model <- geeglm(clus.betas[ind.comp] ~ exposure[ind.comp] ")      
    model.expr.1.site <- paste(model.expr.1.site, ", data = covariates[ind.comp,], id = as.numeric(id)[ind.comp])")
    
    for (i in 1:n.mod){
      clus.probes <- clusters.list[[i]]
      clus.betas <- betas[clus.probes,]
      
      
      if (length(clus.probes) == 1){
        
        ind.comp <- which(complete.cases(cbind(clus.betas, exposure, id)))
        eval(parse(text = model.expr.1.site))
        
        effect[i] <- summary(model)[[6]][2,1]
        se[i] <- summary(model)[[6]][2,2]
        pvals[i] <- summary(model)[[6]][2,4]
        n.sites[i] <- 1
        sites[i] <- clus.probes
        
      } else{
        
        ind.comp <- which(complete.cases(cbind(t(clus.betas), exposure,  id)))
        
        temp.betas <- clus.betas[,ind.comp]
        
        temp.covars <- as.data.frame(cbind(exposure[ind.comp], id[ind.comp]))
        colnames(temp.covars) <- c("exposure",  "id")
        
        temp.long <- organize.island.repeated(temp.betas, temp.covars)
        temp.long <- temp.long[complete.cases(temp.long),]
        temp.long <- temp.long[order(temp.long$id),]
        temp.long$probeID <- as.factor(temp.long$probeID)
        #temp.long$exposure <- as.numeric(as.character(temp.long$exposure))
        temp.long$exposure <- as.factor(temp.long$exposure)
        temp.long$Beta <- as.numeric(as.character(temp.long$Beta))
        
        eval(parse(text = model.expression))
        
        
        effect[i] <- summary(model)[[6]][2,1]
        se[i] <- summary(model)[[6]][2,2]
        pvals[i] <- summary(model)[[6]][2,4]
        n.sites[i] <- length(clus.probes)
        for (c in 1:length(clus.probes)) sites[i] <- paste(sites[i], clus.probes[c], sep = ";")
        substr(sites[i], 1, 1) <- ""
        
      }
    }
    
    
  }
  
  result <- data.frame(exposure_effect_size = effect, exposure_effect_se = se, exposure_pvalue= pvals, n_sites_in_cluster = n.sites, cluster_sites =sites)  
  #result <- result[order(result$exposure_pvalue),]
  
  
  if(!is.null(result.file.name)) write.table(result, file = result.file.name, col.names = T , row.names = F,  append = T)
  
  
  return(result)
  
}

####################################################################

## Analysis for all CpGs (present everywhere in the 1st batch)
setwd("Files_and_code")
betas <- read.csv("all_cpgs.10x_leukocytes_everywhere_n_=36.txt", sep = "\t", header = T)
test <- read.csv("~/Bioinformatics/Spring_Project/CpGs_n_=10_max_min_=20.tsv", sep='\t')
str(betas)
str(test)
splitted <- as.data.frame(str_split_fixed(betas$X, ":", n = 2))
betas$CHR <- splitted$V1
betas$Coordinate_37 <- splitted$V2
annot <- betas[,c(38,39)]
annot$IlmnID <- betas$X
rownames(betas) <- betas$X
betas$X <- NULL
betas$CHR <-NULL
betas$Coordinate_37 <- NULL
rows_for_filtration <- rownames(rbind(betas[rowSums(betas) <= 0,], betas[rowSums(betas) == ncol(betas)*100,]))
betas <- betas[rowSums(betas) > 0,]
betas <- betas[rowSums(betas) < ncol(betas)*100,]
annot <- annot %>% filter(!IlmnID %in% c(rows_for_filtration))
betas <- betas/100
annot.betas <- annot[annot$IlmnID %in% rownames(betas),]
annot.betas <- annot
annot.betas <- annot.betas
annot$Coordinate_37 <- as.numeric(levels(annot$Coordinate_37))[annot$Coordinate_37]
annot$IlmnID <- as.character(annot$IlmnID)
annot.betas$Coordinate_37 <- as.numeric(levels(annot.betas$Coordinate_37))[annot.betas$Coordinate_37]
annot.betas$IlmnID <- as.character(annot.betas$IlmnID)
colnames(betas) <- gsub( "_", "", as.character(colnames(betas)))
colnames(betas)[26] <- "X463lym"

######## A-clustering ##########

chroms <- unique(annot.betas$CHR)
chroms <- as.character(chroms[order(as.numeric(as.character(chroms)))])
betas.by.chrom <- vector(mode = "list", length = length(chroms))
sites.by.chrom <- vector(mode = "list", length = length(chroms))
names(betas.by.chrom) <- names(sites.by.chrom) <- chroms

for (i in 1:length(chroms)) {
  cpg.chrom <- as.character(annot.betas[annot.betas$CHR == chroms[i],]$IlmnID)
  betas.by.chrom[[i]] <- as.matrix(betas[cpg.chrom, ])
  if (ncol(betas.by.chrom[[i]]) == 1) {
    betas.by.chrom[[i]] <- t(betas.by.chrom[[i]])
    rownames(betas.by.chrom[[i]]) <- cpg.chrom
  }
  sites.by.chrom[[i]] <- annot.betas[annot.betas$CHR == chroms[i], ][,c(1:3)]
}

chrom.list = list(betas.by.chrom = betas.by.chrom, sites.locations.by.chrom = sites.by.chrom)
clusters.by.chrom <- vector(mode = "list", length = length(chrom.list[[1]]))


for (i in 1:length(chrom.list[[1]])) {
  betas.temp <- chrom.list[[1]][[i]]
  locations.temp <- chrom.list[[2]][[i]]
  which.clust <- Dbp.merge(t(betas.temp), thresh.dist = 0.25, bp.thresh.dist = 999, as.numeric(locations.temp$Coordinate_37), dist.type = "spearman")
  clust.vec <- Acluster(t(betas.temp), thresh.dist = 0.25, which.clust = which.clust, location.vec = as.numeric(locations.temp$Coordinate_37), max.dist = 1000, type = "average", dist.type = "spearman")
  clusters.by.chrom[[i]] <- lapply(clust.vec, function(x) return(locations.temp$IlmnID[which(clust.vec == x)]))
  clusters.by.chrom[[i]] <- clusters.by.chrom[[i]][which(!duplicated(clusters.by.chrom[[i]]))]
  if (i == 1) 
    clusters.all <- clusters.by.chrom[[1]]
  else clusters.all <- c(clusters.all, clusters.by.chrom[[i]])
}
clusters.list.1st_batch = clusters.all

max_1st <- 2
clusters.list.1st_batch_ <- list()
for (i in 1:length(clusters.all)) {
  if (length(clusters.list.1st_batch[[i]]) > 1) {
    if (i == 1) 
      clusters.list.1st_batch_ <- clusters.list.1st_batch[1]
    else clusters.list.1st_batch_ <- c(clusters.list.1st_batch_, clusters.list.1st_batch[i])
    if (length(clusters.list.1st_batch[[i]]) > max_1st)
      max_1st = length(clusters.list.1st_batch[[i]])
  }
}


library(stringi)
res <- as.data.frame(t(stri_list2matrix(clusters.list.1st_batch_)))
write.csv(res, "Clusters.csv")

############GEE bivariates (expositions only)#############

my_exp_and_cov <- read.csv("smoking_last6months.csv")
list_36_samples <- read.csv("list_36_samples.txt", stringsAsFactors = F, header = F)
my_exp_and_cov <- my_exp_and_cov[,-c(1,3)]
my_exp_and_cov$ID <- paste0(my_exp_and_cov$ID, "lym")
my_exp_and_cov <- my_exp_and_cov[my_exp_and_cov$ID %in% list_36_samples$V1,]
my_exp_and_cov$exp_id <- sub("^", "X", my_exp_and_cov$ID)
my_exp_and_cov <- my_exp_and_cov[match(colnames(betas), my_exp_and_cov$exp_id),]
identical(colnames(betas), my_exp_and_cov$exp_id) ##must be TRUE before proceeding to next stage

#for clusters
exposure <- my_exp_and_cov$smoking_last6months

GEE.results.clusters <- GEE.clusters(betas, clusters.list.1st_batch, exposure, id = colnames(betas), working.cor = "ex")
GEE_results_smoking_last6months <- GEE.results.clusters
GEE_results_smoking_last6months$exposure_padjusted <- p.adjust(GEE_results_smoking_last6months$exposure_pvalue, "BH")

GEE_results_last6months_sig <- subset(GEE_results_smoking_last6months, exposure_padjusted <= 0.05)
GEE_results_last6months_nan <- GEE_results_smoking_last6months[is.na(GEE_results_smoking_last6months$exposure_pvalue) == T,]
GEE_results_last6months_zero <- GEE_results_smoking_last6months[GEE_results_smoking_last6months$exposure_pvalue == 0,]

fileName <- paste("./", "last6months", "_GEEclusters_Bvalues.tsv", sep="")
write.table(GEE_results_smoking_last6months,file=fileName,col.names=TRUE,row.names=F,sep="\t")
fileName <- paste("./", "last6months", "_GEEclusters_Bvalues_sig.tsv", sep="")
write.table(GEE_results_last6months_sig,file=fileName,col.names=TRUE,row.names=F,sep="\t")