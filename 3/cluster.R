rm(list = ls())

library(dplyr)
# library(BiocManager)
# BiocManager::install("ConsensusClusterPlus")
library(ConsensusClusterPlus)
library(data.table)
library(tidyverse)
dat <- data.table::fread("GSE87466_Matrix_norm.csv") %>% 
  as.data.frame() 
rownames(dat) <- dat$V1
dat <- dat[,-1]

# dat <- data.table::fread("GSE87466_Matrix_norm.csv") %>% 
#   as.data.frame() %>% 
#   remove_rownames() %>% 
#   tibble::column_to_rownames(.,var = "V1")
# 
# 
# range(dat)
# dat <- log2(dat+1)
# range(dat)
# write.csv(x = dat,file = "OV_norm.csv")


gene <- read.table("intersect.txt",header = T)
dat <- dat[gene$id,]





cluster <- ConsensusClusterPlus(d = dat %>% as.matrix(), 
                                maxK = 6, 
                                reps = 100,  
                                pItem = 0.8, 
                                pFeature = 1, 
                                
                                clusterAlg = 'pam', 
                                title = 'consensus_cluster_pdf',
                                innerLinkage = 'average', finalLinkage = 'average',
                                
                                distance = 'pearson', 
                                plot = "pdf", writeTable = T, seed = 500,
                                verbose = T, 
                                
                                corUse = 'everything')


icl <- calcICL(cluster, title = "cluster-consensus_pdf", plot = "pdf", writeTable = TRUE)

