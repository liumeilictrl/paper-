# BiocManager:: install("DESeq2")
library(data.table)
library(DESeq2)
library(tidyverse)
exp_nc <- fread('./output/LUSC_Counts.csv')
# exp_nc <- exp_nc[,-1]
# exp_nc <- exp_nc[!duplicated(exp_nc$gene),]
exp_nc <- exp_nc %>% column_to_rownames(var ='V1')
group <- ifelse(substr(colnames(exp_nc), 14, 14) == "0", "tumor", "normal")
colData <- data.frame(id = colnames(exp_nc), group = group)
table(colData$group)
write.csv(colData, "./output/LUSC_group.csv", row.names = F)

exp_nc <- exp_nc[,colnames(exp_nc) %in% colData$id]
# 

# 
exp_nc <- exp_nc[,colData$id]
identical(colnames(exp_nc), colData$id)
colnames(colData)[2] <- 'Type'
# 
colData$Type <- factor(colData$Type, levels = c('tumor','normal'))
# 
exp_nc <- round(exp_nc,0)
dds <- DESeqDataSetFromMatrix(countData = exp_nc, colData = colData, design = ~Type)
# 
dds <- dds[rowSums(counts(dds)) > 1,]

# 
vsd <- vst(dds, blind = FALSE)
# 
plotPCA(vsd, "Type")
# 
dds <- DESeq(dds)
# 
res <- results(dds)
# 
# 
contrast = c("Type", "tumor", "normal")
# 
dd1 <- results(dds, contrast = contrast, alpha = 0.05)
# 
plotMA(dd1, ylim = c(-5,5))
# 1ogFC
# dd2 <- lfcShrink(dds, contrast = contrast, res = dd1, type = "ashr")
# plotMA(dd2, ylim = c(-5,5))

# 
res <- dd1 %>%
as.data.frame()

write.csv(res,'output/diff_all.csv')


