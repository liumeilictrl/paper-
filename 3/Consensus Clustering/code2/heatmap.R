rm(list = ls())
library(data.table)
library(pheatmap)
library(magrittr)
library(tibble)
library(ggplot2)
dat <- fread("input/LUSC_TPM_tumor.csv", data.table = F)
dat <- dat %>% column_to_rownames(var ='V1')
genes <- fread("input/PRDEGs.csv", data.table = F)
genes <- genes$V1

clinical <- fread("input/LUSC_clinical_cleaned.csv", data.table = F)
cluster <- fread("output/cluster.csv", data.table = F)
colnames(cluster)[1] <- "sample_id"
clinical <- dplyr::left_join(clinical, cluster)
write.csv(clinical, "output/clinical_cluster.csv")
clinical1 <- clinical
clinical1$Age <- ifelse(clinical1$Age <= 60, "<=60", ">60") %>% 
  factor(levels = c("<=60", ">60"))
clinical1$Status <- ifelse(clinical1$OS == 1, "Alive", "Dead") %>% 
  factor(levels = c("Dead", "Alive"))
clinical1$Stage <- clinical1$Stage %>% 
  factor(levels = c("Stage I", "Stage II", "Stage III", "Stage IV"))
clinical1 <- dplyr::arrange(clinical1, cluster)

dat <- dat[genes,clinical1$sample_id]
col1 <- colorRampPalette(colors = c("blue","white","red"))(50) # 
annotation_col = data.frame(Cluster = clinical1$cluster,
                            Stage = clinical1$Stage,
                            Age = clinical1$Age,
                            Status = clinical1$Status)
rownames(annotation_col) <- clinical1$sample_id
ann_colors = list(Cluster = c(Cluster1 = "#4DBBD5", Cluster2 = "#E64B35"),
                  Stage = c(`Stage I` = "#00A087", `Stage II` = "#3C5488", `Stage III` = "#F39B7F", `Stage IV` = "#8491B4"),
                  Age = c(`<=60` = "#91D1C2", `>60` = "#DC0000"),
                  Status = c(`Alive` = "#7E6148", `Dead` = "#B09C85")) # 

p <- pheatmap(dat, cluster_col = F, cluster_rows = T, # 
              scale = "row",
              color = col1, # 
              border_color = NA, # 
              show_rownames = T, show_colnames = F,#  
              fontsize = 6, # 
              annotation_col = annotation_col,
              annotation_colors = ann_colors) # 

ggsave(plot = p, file = "output/Heatmap.pdf", width = 16, height = 10, units = "cm")

