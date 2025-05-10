#install.packages("colorspace")
#install.packages("stringi")
#install.packages("ggplot2")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("DOSE")
#BiocManager::install("clusterProfiler")
#BiocManager::install("enrichplot")

 
library("clusterProfiler")   
library("org.Hs.eg.db")      
library("enrichplot")        
library("ggplot2")           

 
inputFile="intersect.txt"  
setwd("D:\\")  

 
rt=read.table(inputFile, sep="\t", check.names=F, header=F)  
genes=as.vector(rt[,1])  

 
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)  
entrezIDs <- as.character(entrezIDs)  


out=cbind(rt, entrezID=entrezIDs)  
colnames(out)[1]="Gene"  
write.table(out, file="id.txt", sep="\t", quote=F, row.names=F)  

 
pvalueFilter=0.05  
qvalueFilter=1     

 
rt=read.table("id.txt", sep="\t", header=T, check.names=F)  
rt=rt[is.na(rt[,"entrezID"])==F,]  
colnames(rt)[1]="Gene"  
gene=rt$entrezID  

 
colorSel="qvalue"  
if(qvalueFilter>0.05){
  colorSel="pvalue"  
}


kk <- enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff = 1, qvalueCutoff = 1)  
KEGG=as.data.frame(kk)  

 
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt$Gene[match(strsplit(x,"/")[[1]],as.character(rt$entrezID))],collapse="/")))  # 转换基因ID为基因名

 
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$qvalue<qvalueFilter),]  
write.table(KEGG, file="KEGG.txt", sep="\t", quote=F, row.names = F)  

 
showNum=30  
if(nrow(KEGG)<showNum){
  showNum=nrow(KEGG)  
}

pdf(file="barplot.pdf", width = 9, height = 11)  
barplot(kk, drop = TRUE, showCategory = showNum, color = colorSel)  
dev.off()  

 
pdf(file="bubble.pdf", width = 9, height = 11)  
dotplot(kk, showCategory = showNum, orderBy = "GeneRatio", color = colorSel)  
dev.off() 