library(data.table)
library(tidyverse)
# install.packages("rjson")
library(rjson)
# 
json <- jsonlite::fromJSON("./input/TCGA_RNAseq/metadata.cart.2024-05-21.json")
View(json)
# 
sample_id <- sapply(json$associated_entities, function(x){x[,1]})
case_id <- sapply(json$associated_entities,function(x){x[,3]})
sample_case <- t(rbind(sample_id,case_id))

clinical <- read.delim('input/TCGA_clinical/clinical.tsv',header = T)
clinical <- clinical[!duplicated(clinical$case_id),]
rownames(clinical) <- NULL
View(clinical)
# 
clinical_matrix <- merge(sample_case, clinical, by="case_id", all.x = T)
clinical_matrix <- clinical_matrix[,-1]
write.csv(clinical_matrix,'./output/clinical_matrix.csv',row.names = FALSE)
View(clinical_matrix)
colnames(clinical_matrix)[grep("age|gender|pathologic|stage", colnames(clinical_matrix))]
clinical_matrix$OS <- ifelse(grepl("Dead", clinical_matrix$vital_status), 1, 0) %>% as.numeric()
clinical_matrix$OS.time <- ifelse(
  clinical_matrix$vital_status == "Dead",
  clinical_matrix$days_to_death,
  clinical_matrix$days_to_last_follow_up
) %>% as.numeric()
clinical_matrix1 <- clinical_matrix[,c("sample_id", "OS", "OS.time", "age_at_index", "gender", "ajcc_pathologic_stage", "ajcc_pathologic_t", "ajcc_pathologic_m", "ajcc_pathologic_n")]
# clinical_matrix1 <- clinical_matrix[,c("sample_id", "OS", "OS.time", "age_at_index", "gender", "ajcc_pathologic_stage", "ajcc_pathologic_t", "ajcc_pathologic_m", "ajcc_pathologic_n", "masaoka_stage")]
str(clinical_matrix1)
clinical_matrix1$age_at_index <- as.numeric(clinical_matrix1$age_at_index)
write.csv(clinical_matrix1, './output/LUSC_clinical.csv', row.names = F)

