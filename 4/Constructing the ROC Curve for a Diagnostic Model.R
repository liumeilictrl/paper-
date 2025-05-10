# install.packages("pROC")
library(pROC)

 
rsFile = "model.riskMatrix.txt"      
method = "Lasso+Stepglm[both]"          
#setwd("")      


riskRT = read.table(rsFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)

 
CohortID = gsub("(.*)\\_(.*)\\_(.*)", "\\1", rownames(riskRT))
CohortID = gsub("(.*)\\.(.*)", "\\1", CohortID)
riskRT$Cohort = CohortID

 
for (Cohort in unique(riskRT$Cohort)) {
   
  rt = riskRT[riskRT$Cohort == Cohort,]
  y = gsub("(.*)\\_(.*)\\_(.*)", "\\3", row.names(rt))
  y = ifelse(y == "Normal", 0, 1)
  
   
  roc1 = roc(y, as.numeric(rt[, method]))  
  ci1 = ci.auc(roc1, method = "bootstrap") 
  ciVec = as.numeric(ci1)
  
   
  pdf(file = paste0("ROC.", Cohort, ".pdf"), width = 5, height = 4.75)
  #plot(roc1, print.auc = TRUE, col = "#1f77b4", legacy.axes = TRUE, main = Cohort) 
  plot(roc1, print.auc = TRUE, col = "#1f77b4", legacy.axes = TRUE, main = Cohort, print.auc.col = "#000000",font.main = 1)
  
   
  polygon(c(1, roc1$specificities, 0), c(0, roc1$sensitivities, 0), col = rgb(31/255, 119/255, 180/255, 0.2), border = NA) 
  
  text(0.39, 0.43, paste0("95% CI: ", sprintf("%.03f", ciVec[1]), "-", sprintf("%.03f", ciVec[3])), col = "#000000") # 
  dev.off()  
}









#library(pROC)

 
nCohorts <- length(unique(riskRT$Cohort))


par(mfrow = c(nCohorts, 1))  

for (Cohort in unique(riskRT$Cohort)) {
  
  rt = riskRT[riskRT$Cohort == Cohort, ]
  y = gsub("(.*)\\_(.*)\\_(.*)", "\\3", row.names(rt))
  y = ifelse(y == "Normal", 0, 1)
  
  
  roc1 = roc(y, as.numeric(rt[, method]))  
  ci1 = ci.auc(roc1, method = "bootstrap") 
  ciVec = as.numeric(ci1)
  
  
  plot(roc1, print.auc = TRUE, col = "#1f77b4", legacy.axes = TRUE, 
       main = Cohort, print.auc.col = "#000000", font.main = 1)
  
  
  polygon(c(1, roc1$specificities, 0), c(0, roc1$sensitivities, 0), 
          col = rgb(31/255, 119/255, 180/255, 0.2), border = NA) 
  
  
  text(0.39, 0.43, paste0("95% CI: ", sprintf("%.03f", ciVec[1]), "-", sprintf("%.03f", ciVec[3])), 
       col = "#000000")
}