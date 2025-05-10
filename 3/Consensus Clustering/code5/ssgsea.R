rm(list = ls())
library(data.table)
library(tibble)
library(GSVA)
library(dplyr)
dat <- fread("input/LUSC_TPM_tumor.csv", data.table = F)
dat <- dat %>% column_to_rownames(var ='V1')
clinical <- fread("input/clinical_lasso.csv", data.table = F)
group <- clinical[,c("sample_id", "risk")]
dat <- dat[,group$sample_id]

cell_gene <- openxlsx::read.xlsx("input/16cells.xlsx")
cell_list <- as.list(cell_gene) %>% lapply(function(x){x[!is.na(x)]})
gsva_cell <- gsva(as.matrix(dat),cell_list,method="ssgsea",abs.ranking=F)
write.csv(gsva_cell, "output/gsva_cell.csv")

cb2 <- as.data.frame(t(gsva_cell))
immu <- rep(colnames(cb2),each = nrow(cb2)) 
immu <- factor(immu) 
a <- c(group$risk)
group <- rep(a, ncol(cb2)) 
group <- factor(group,levels = c('low','high')) 
value <- c() 
for (j in 1:ncol(cb2)) { value<-c(value,cb2[,j])}
value<-as.numeric(value)
Data <- data.frame(immu_cell=immu,group=group,value=value) 
library(ggplot2)
library(corrgram)
library(ggthemes)
library(ggpubr)
#
immu_pvalue <- compare_means(value~group,data = Data,group.by = "immu_cell")
immu_pvalue <- immu_pvalue %>% arrange(immu_cell)
write.csv(immu_pvalue, "output/immu_pvalue.csv")

p<-ggplot(Data,aes(x=immu_cell,y=value,fill=group))+
  geom_boxplot(width=0.7,size=0.3,outlier.color = NA,linewidth=0.1,fatten=1,position=position_dodge(0.85))+
  theme_bw()+scale_fill_manual(values = c("#00468B","#ED0000"))+#
  theme(panel.grid = element_blank())+
  stat_compare_means(symnum.args = list(cutpoints = c(0,0.001, 0.01, 0.05, 1),
                                        symbols = c("***", "**", "*", "ns")),label = "p.signif",size=6,hjust = 0.5,vjust = 0)+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))+
  theme(legend.position = 'top')+xlab('')+ylab('Infiltration Abundance')+labs(fill='Group')

p
ggsave("output/barplot.pdf",plot=p,width = 10,height = 8)
dev.off()
