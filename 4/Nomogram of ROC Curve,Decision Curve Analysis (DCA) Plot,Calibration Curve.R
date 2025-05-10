
#install.packages("rmda")

 
library(rms)
library(rmda)

inputFile="data.train.txt"             
geneFile="gene.txt"     


 
data=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
data=t(data)  
row.names(data)=gsub("-", "_", row.names(data))  

 
geneRT=read.table(geneFile, header=T, sep="\t", check.names=F)
data=data[as.vector(geneRT[,1]),]  

 
data=t(data)  
group=gsub("(.*)\\_(.*)\\_(.*)", "\\3", row.names(data))  
rt=cbind(as.data.frame(data), Type=group)  
paste(colnames(data), collapse="+")  

 
ddist=datadist(rt)
options(datadist="ddist")  

 
lrmModel=lrm(Type~ AQP9+	RPS4Y1+	CXCL9+	DUOX2+	VNN2+	PCK1+	MMP1+	SELL, data=rt, x=T, y=T)
nomo=nomogram(lrmModel, fun=plogis,
              fun.at=c(0.0001,0.1,0.3,0.5,0.7,0.9,0.99),
              lp=F, funlabel="Risk of Disease")  
pdf("Nomo.pdf", width=8, height=6)  
plot(nomo)  
dev.off()

 
cali=calibrate(lrmModel, method="boot", B=1000)  
pdf("Calibration.pdf", width=5.5, height=5.5)  
plot(cali,
     xlab="Predicted probability",
     ylab="Actual probability", sub=F)  
dev.off()
 
rt$Type=ifelse(rt$Type=="Normal", 0, 1)  
dc=decision_curve(Type ~  AQP9+	RPS4Y1+	CXCL9+	DUOX2+	VNN2+	PCK1+	MMP1+	SELL, data=rt,
                  family = binomial(link ='logit'),
                  thresholds= seq(0,1,by = 0.01),
                  confidence.intervals = 0.95)  
pdf(file="DCA.pdf", width=5.5, height=5.5)  
plot_decision_curve(dc,
                    curve.names="Model",
                    xlab="Threshold probability",
                    cost.benefit.axis=T,
                    col="red",
                    confidence.intervals=FALSE,
                    standardize=FALSE)  
dev.off()
