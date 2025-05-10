rm(list = ls())
library(survival)
library(survminer)
library(data.table)
library(magrittr)
rt <- fread("output/clinical_cluster.csv")
rt <- rt[,c("OS", "OS.time", "cluster")]
rt$OS.time <- rt$OS.time/365

diff <- survdiff(Surv(OS.time, OS) ~ cluster, data = rt)
pValue <- 1 - pchisq(diff$chisq, df = 1)
if(pValue < 0.001){
  pValue <- "p < 0.001"
}else{
  pValue <- paste0("p = ", sprintf("%.03f", pValue))
}
fit <- survfit(Surv(OS.time, OS) ~ cluster, data = rt)

# 
surPlot <- ggsurvplot(fit, 
                      data = rt,
                      conf.int = TRUE,
                      pval = pValue,
                      pval.size = 5,
                      legend.labs = c("Cluster1", "Cluster2"),
                      legend.title = "Cluster",
                      xlab = "Time(years)",
                      break.time.by = 1,
                      risk.table.title = "",
                      palette = c("red","blue"),
                      risk.table = T,
                      risk.table.height = .25)
pdf(file = "output/km_cluster.pdf", onefile = F, width = 6, height = 5)
print(surPlot)
dev.off()
