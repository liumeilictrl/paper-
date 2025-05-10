rm(list = ls())
library(data.table)
library(tidyverse)

dat <- fread("input/LUSC_TPM_tumor.csv", data.table = F)
dat <- dat %>% column_to_rownames(var ='V1')
genes <- fread("input/PRDEGs.csv", data.table = F)
genes <- genes$V1
dat <- dat[genes,] %>% as.matrix()
# BiocManager::install('ConsensusClusterPlus')
library(ConsensusClusterPlus)
results = ConsensusClusterPlus(dat,
                               maxK = 9,
                               reps = 500,
                               pItem = 0.8,
                               pFeature = 1,
                               title = "output",
                               clusterAlg = "km",
                               distance = "euclidean",
                               seed = 123,
                               plot = "pdf")

cluster <- results[[2]][["consensusClass"]] %>% as.data.frame()
colnames(cluster) <- "cluster"
cluster$cluster <- ifelse(cluster$cluster == 1, "Cluster1", "Cluster2")
cluster$sample_id <- rownames(cluster)
cluster <- cluster %>% 
  dplyr::select(sample_id, cluster)
write.csv(cluster, "output/cluster.csv", row.names = F)
