rm(list = ls())
library(data.table)
library(pheatmap)
library(tibble)
dat <- fread("output/LUSC_TPM.csv", data.table = F)
dat <- dat %>% column_to_rownames(var ='V1')
genes <- fread("input/Pyroptosis related genes.csv", data.table = F)
genes <- genes$genes
res <- fread("output/diff_all.csv", data.table = F)
res <- res %>% column_to_rownames(var ='V1')
dat <- dat[genes,]
res <- res[genes,]
res$psig <- ifelse((abs(res$log2FoldChange) > 1) & (res$padj < 0.001), paste0(rownames(res), "***"),
                   ifelse((abs(res$log2FoldChange) > 1) & (res$padj < 0.01), paste0(rownames(res), "**"),
                          ifelse((abs(res$log2FoldChange) > 1) & (res$padj < 0.05), paste0(rownames(res), "*"),
                                 rownames(res))))
res1 <- res %>% 
  dplyr::filter(abs(log2FoldChange) > 1 & padj < 0.05)
write.csv(res1, "output/PRDEGs.csv") 
 
rownames(dat) <- res$psig
group <- data.table::fread("output/LUSC_group.csv", data.table = F)
group <- dplyr::arrange(group, group) 
dat <- dat[,group$id]
col1 <- colorRampPalette(colors = c("blue","white","red"))(50) 
annotation_col = data.frame(Group = group$group)
rownames(annotation_col) <- group$id
ann_colors = list(Group = c(normal = "blue", tumor = "red")) 
p <- pheatmap(dat, cluster_col = F, cluster_rows = T, 
              scale = "row",
              color = col1, # 
              border_color = NA, 
              show_rownames = T, show_colnames = F,
              fontsize = 6, 
              annotation_col = annotation_col, annotation_colors = ann_colors) 

ggsave(plot = p, file = "output/Heatmap.pdf", width = 16, height = 10, units = "cm")
