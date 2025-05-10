rm(list = ls())
library(data.table)
library(magrittr)
clinical <- fread("output/LUSC_clinical.csv", data.table = F)
clinical1 <- clinical %>% 
  dplyr::filter(stringr::str_sub(sample_id,14,14) == 0)
clinical1 <- clinical1[,1:6]
colnames(clinical1)[4:6] <- c("Age","Gender","Stage")
table(clinical1$Stage)
clinical1$Stage[clinical1$Stage == "Stage IA"] <- "Stage I"
clinical1$Stage[clinical1$Stage == "Stage IB"] <- "Stage I"
clinical1$Stage[clinical1$Stage == "Stage IIA"] <- "Stage II"
clinical1$Stage[clinical1$Stage == "Stage IIB"] <- "Stage II"
clinical1$Stage[clinical1$Stage == "Stage IIIA"] <- "Stage III"
clinical1$Stage[clinical1$Stage == "Stage IIIB"] <- "Stage III"
clinical1$Stage[clinical1$Stage == "'--"] <- NA
table(clinical1$Stage)
str(clinical1) # 
clinical1 <- clinical1 %>% 
  dplyr::filter(!is.na(OS)) %>% 
  dplyr::filter(OS.time > 0)
write.csv(clinical1, 'output/LUSC_clinical_cleaned.csv', row.names = F)

dat_counts <- fread("output/LUSC_Counts_matrix.csv", data.table = F)
dat_counts <- dat_counts %>% 
  dplyr::distinct(gene, .keep_all = T)
rownames(dat_counts) <- dat_counts$gene
dat_counts <- dat_counts[,-c(1:2)]
sample_normal <- colnames(dat_counts)[stringr::str_sub(colnames(dat_counts) ,14,14) == 1]
sample <- c(clinical1$sample_id, sample_normal)
dat_counts <- dat_counts[,sample]
write.csv(dat_counts, "output/LUSC_Counts.csv")
dat_counts_tumor <- dat_counts[,clinical1$sample_id]
write.csv(dat_counts_tumor, "output/LUSC_Counts_tumor.csv")

dat_tpm <- fread("output/LUSC_TPM_matrix.csv", data.table = F)
dat_tpm <- dat_tpm %>% 
  dplyr::distinct(gene, .keep_all = T)
rownames(dat_tpm) <- dat_tpm$gene
dat_tpm <- dat_tpm[,-c(1:2)]
sample_normal <- colnames(dat_tpm)[stringr::str_sub(colnames(dat_tpm) ,14,14) == 1]
sample <- c(clinical1$sample_id, sample_normal)
dat_tpm <- dat_tpm[,sample]
dat_tpm <- log2(dat_tpm + 1)
write.csv(dat_tpm, "output/LUSC_TPM.csv")

dat_tpm_tumor <- dat_tpm[,clinical1$sample_id]
write.csv(dat_tpm_tumor, "output/LUSC_TPM_tumor.csv")
