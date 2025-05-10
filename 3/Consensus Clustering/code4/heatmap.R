rm(list = ls())
library(data.table)
library(pheatmap)
library(magrittr)
library(tibble)
library(ggplot2)
dat <- fread("input/LUSC_TPM_tumor.csv", data.table = F)
dat <- dat %>% column_to_rownames(var ='V1')
genes <- fread("input/geneCoef.csv", data.table = F)
genes <- genes$Gene

clinical <- fread("input/clinical_lasso.csv", data.table = F)
clinical <- clinical[,-1]
clinical1 <- clinical
clinical1$Age <- ifelse(clinical1$Age <= 60, "<=60", ">60") %>% 
  factor(levels = c("<=60", ">60"))
clinical1$Status <- ifelse(clinical1$OS == 1, "Alive", "Dead") %>% 
  factor(levels = c("Dead", "Alive"))
clinical1$Stage <- clinical1$Stage %>% 
  factor(levels = c("Stage I", "Stage II", "Stage III", "Stage IV"))
clinical1 <- dplyr::arrange(clinical1, risk)

dat <- dat[genes,clinical1$sample_id]
col1 <- colorRampPalette(colors = c("blue","white","red"))(50) 
annotation_col = data.frame(risk = clinical1$risk,
                            Stage = clinical1$Stage,
                            Age = clinical1$Age,
                            Status = clinical1$Status)
rownames(annotation_col) <- clinical1$sample_id
ann_colors = list(risk = c(low = "#4DBBD5", high = "#E64B35"),
                  Stage = c(`Stage I` = "#00A087", `Stage II` = "#3C5488", `Stage III` = "#F39B7F", `Stage IV` = "#8491B4"),
                  Age = c(`<=60` = "#91D1C2", `>60` = "#DC0000"),
                  Status = c(`Alive` = "#7E6148", `Dead` = "#B09C85")) 

p <- pheatmap(dat, cluster_col = F, cluster_rows = T, 
              scale = "row",
              color = col1, 
              border_color = NA, 
              show_rownames = T, show_colnames = F,
              fontsize = 6, 
              annotation_col = annotation_col,
              annotation_colors = ann_colors) 

ggsave(plot = p, file = "output/Heatmap.pdf", width = 16, height = 10, units = "cm")

