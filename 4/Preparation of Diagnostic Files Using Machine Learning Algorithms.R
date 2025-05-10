#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("sva")


library(limma)  
library(sva)    

geneFile = "intersect.txt"  
#setwd("D:\\")  

 
files = dir()  
files = grep("normalize.txt$", files, value = TRUE)  

geneList = list()  

 
for (file in files) {
  rt = read.table(file, header = TRUE, sep = "\t", check.names = FALSE)  
  geneNames = as.vector(rt[, 1])  
  uniqGene = unique(geneNames)   
  header = unlist(strsplit(file, "\\.|\\-"))  
  geneList[[header[1]]] = uniqGene  
}

 
interGenes = Reduce(intersect, geneList)  
allTab = data.frame()  
batchType = c()        

 
for (i in 1:length(files)) {
  inputFile = files[i]  # 
  header = unlist(strsplit(inputFile, "\\.|\\-"))   
  rt = read.table(inputFile, header = TRUE, sep = "\t", check.names = FALSE)  
  rt = as.matrix(rt)  
  rownames(rt) = rt[, 1]  
  exp = rt[, 2:ncol(rt)]  
  dimnames = list(rownames(exp), colnames(exp))  
  data = matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)  
  rt = avereps(data)  
  colnames(rt) = paste0(header[1], "_", colnames(rt))  
  
   
  if (i == 1) {
    allTab = rt[interGenes, ]  
  } else {
    allTab = cbind(allTab, rt[interGenes, ])  
  }
  batchType = c(batchType, rep(i, ncol(rt))) 
}

 
svaTab = ComBat(allTab, batchType, par.prior = TRUE)   


geneRT = read.table(geneFile, header = FALSE, sep = "\t", check.names = FALSE)  
 
geneTab = svaTab[intersect(row.names(svaTab), as.vector(geneRT[, 1])), ]
geneTab = t(geneTab)  

 
train = grepl("^merge", rownames(geneTab), ignore.case = TRUE)  
trainExp = geneTab[train, , drop = FALSE]  
testExp = geneTab[!train, , drop = FALSE]  
 
rownames(trainExp) = gsub("merge_", "Train.", rownames(trainExp))
 
trainType = gsub("(.*)\\_(.*)\\_(.*)", "\\3", rownames(trainExp))
testType = gsub("(.*)\\_(.*)\\_(.*)", "\\3", rownames(testExp))
trainType = ifelse(trainType == "Normal", 0, 1)  #######################################################
testType = ifelse(testType == "Normal", 0, 1)  ########################################################
 
trainExp = cbind(trainExp, Type = trainType)
testExp = cbind(testExp, Type = testType)

 
trainOut = rbind(id = colnames(trainExp), trainExp)  
write.table(trainOut, file = "data.train.txt", sep = "\t", quote = FALSE, col.names = FALSE)  
testOut = rbind(id = colnames(testExp), testExp)  
write.table(testOut, file = "data.test.txt", sep = "\t", quote = FALSE, col.names = FALSE)  
