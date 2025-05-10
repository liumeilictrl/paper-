library(data.table)
library(tidyverse)
# install.packages("rjson")
library(rjson)
# 
json <- jsonlite::fromJSON("./input/TCGA_RNAseq/metadata.cart.2024-05-21.json")
View(json)
# 
sample_id <- sapply(json$associated_entities, function(x){x[,1]})
file_sample <- data.frame(sample_id, file_name = json$file_name)

# 
count_file <- list.files('./input/TCGA_RNAseq/', pattern = '*.tsv', recursive = TRUE)
# 
count_file_name <- strsplit(count_file, split='/')
count_file_name <- sapply(count_file_name, function(x){x[2]})

#
matrix_counts = data.frame(matrix(nrow = 60660, ncol = 0))
matrix_fpkm = data.frame(matrix(nrow = 60660, ncol = 0))
matrix_tpm = data.frame(matrix(nrow = 60660, ncol = 0))

# 
for (i in 1:length(count_file)){
  # i = 1
  path = paste0('./input/TCGA_RNAseq/', count_file[i])
  data <- read.delim(path, fill = TRUE, header = FALSE, row.names = 1) %>% as.data.frame()
  colnames(data) <- data[2,]
  data <- data[-c(1:6),]
  rownm <- data$gene_name
  data_counts <- data[3] # 
  data_fpkm <- data[7] # 
  data_tpm <- data[6] # 
  colnames(data_counts) <- file_sample$sample_id[which(file_sample$file_name == count_file_name[i])]
  colnames(data_fpkm) <- file_sample$sample_id[which(file_sample$file_name == count_file_name[i])]
  colnames(data_tpm) <- file_sample$sample_id[which(file_sample$file_name == count_file_name[i])]
  matrix_counts <- cbind(matrix_counts, data_counts)
  matrix_fpkm <- cbind(matrix_fpkm, data_fpkm)
  matrix_tpm <- cbind(matrix_tpm, data_tpm)
}

# 
matrix_counts$gene <- rownm
matrix_fpkm$gene <- rownm
matrix_tpm$gene <- rownm

# 
matrix_counts <- matrix_counts %>% select(gene, everything())
matrix_fpkm <- matrix_fpkm %>% select(gene, everything())
matrix_tpm <- matrix_tpm %>% select(gene, everything())

# 
write.csv(matrix_counts, './output/LUSC_Counts_matrix.csv', row.names = TRUE)
write.csv(matrix_fpkm, './output/LUSC_FPKM_matrix.csv', row.names = TRUE)
write.csv(matrix_tpm, './output/LUSC_TPM_matrix.csv', row.names = TRUE)

