# install.packages("VennDiagram")  

library(VennDiagram)  
setwd("D:\\")   
files = dir()  
files = grep("txt", files, value = TRUE)  
geneList = list()  

 
for(i in 1:length(files)) {
  inputFile = files[i]
  if(inputFile == "intersect.txt") { next }  
  rt = read.table(inputFile, header = FALSE)  
  header = unlist(strsplit(inputFile, "\\.|\\-"))  
  geneList[[header[1]]] = as.vector(rt[,1])  
  uniqLength = length(unique(as.vector(rt[,1])))  
  print(paste(header[1], uniqLength, sep = " "))  
}

 
venn.plot = venn.diagram(
  geneList,
  filename = NULL,
  fill = adjustcolor(rainbow(length(geneList)), alpha.f = 0.6),  
  cex = 1.5,  
  cat.cex = 1.5,  
  scaled = FALSE  
)
pdf(file = "venn.pdf", width = 8, height = 8)  
grid.draw(venn.plot)
dev.off()  


intersectGenes = Reduce(intersect, geneList)  
write.table(file = "intersect.txt", intersectGenes, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)  
