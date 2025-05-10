rm(list = ls())
library(data.table)
library(tibble)
library(magrittr)
library(survival)

dat <- fread("input/LUSC_TPM_tumor.csv", data.table = F)
dat <- dat %>% column_to_rownames(var ='V1')
genes <- fread("input/PRDEGs.csv", data.table = F)
genes <- genes$V1
clinical <- fread("input/clinical_cluster.csv", data.table = F)
clinical <- clinical[,-1]
dat_cox <- dat[genes,] %>% t() %>% data.frame()
dat_cox$sample_id <- rownames(dat_cox)
dat_cox <- dplyr::left_join(dat_cox, clinical) %>% 
  dplyr::select(sample_id, OS, OS.time, everything())
dat_cox$OS.time <- dat_cox$OS.time/365

##1.uni_cox ####
cox_extr <- function(fit){ 
  fit_summary <- summary(fit)
  dat_res <- data.frame(
    row.names = rownames(fit_summary$coef),
    p.value = signif(fit_summary$coef[,"Pr(>|z|)"]),
    mean = signif(fit_summary$coef[,"exp(coef)"]),
    lower = signif(fit_summary$conf.int[,"lower .95"]),
    upper = signif(fit_summary$conf.int[,"upper .95"]),
    coef = signif(fit_summary$coef[,"coef"]),
    check.rows = F
  )
  dat_res <- signif(dat_res,digits = 3)
  return(dat_res)
}
###1.1 
ormulas <- sapply(dat_cox[,genes] %>% colnames(), 
                  function(x) as.formula(paste0('Surv(OS.time,OS)~', paste0("`", x, "`", sep = ""))))
res_uni <- lapply(ormulas, function(x){coxph(x, data = dat_cox)}) %>% 
  lapply(cox_extr) %>% 
  c(use.names = F) %>% 
  do.call(what = rbind)
res_uni1 <- res_uni %>% 
  dplyr::filter(p.value < 0.1) %>% 
  dplyr::mutate(Factor = rownames(.), HR = paste(.$mean, " [", .$lower, " - ", .$upper, "]", sep = "")) %>%
  dplyr::select("Factor", "p.value", "HR", "mean", "lower", "upper")
write.csv(res_uni1, "output/res_uni.csv", row.names = F)

###1.2 forest
labeltext <- rbind(colnames(res_uni1)[1:3], res_uni1[, 1:3])
HR <- rbind(rep(NA, 3), res_uni1[,4:6])
library(forestplot)
pdf("output/forest_uni.pdf", height = 6, width = 10, onefile = F)
forestplot(labeltext = labeltext,
           HR,
           zero = 1, 
           lwd.ci = 1, 
           title ="univariable  Cox regression analysis",
           # colgap = unit(5, 'mm'),
           xlab = "HR", 
           txt_gp = fpTxtGp(ticks = gpar(cex = 0.8), xlab = gpar(cex = 0.8),
                            title = gpar(cex = 1.2), cex = 1) ,
           boxsize = 0.2, 
           xlog = T,
           col = fpColors(box = "#3B4992", line = "black", zero = "black"))
dev.off()

##2. LASSO ####
library(glmnet)
set.seed(123) 
rt <- dat_cox[,c("OS", "OS.time", res_uni1$Factor)]
x = as.matrix(rt[,c(3:ncol(rt))])
y = data.matrix(Surv(rt$OS.time, rt$OS))
fit=glmnet(x, y, family = "cox", maxit = 1000)
pdf("output/lasso.pdf", width = 6, height = 5)
plot(fit, xvar = "lambda", label = TRUE)

cvfit = cv.glmnet(x, y, family="cox", maxit = 1000)
plot(cvfit)
#
abline(v = log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")
dev.off()
###4. 
###4.1 
coef=coef(fit, s = cvfit$lambda.min)
index=which(coef != 0)
actCoef=coef[index]
lassoGene=row.names(coef)[index]
geneCoef=cbind(Gene=lassoGene,Coef=actCoef)
geneCoef   
write.csv(geneCoef, "output/geneCoef.csv")

###2.1 
FinalGeneExp = rt[,lassoGene]
myFun = function(x){crossprod(as.numeric(x),actCoef)}
riskScore = apply(FinalGeneExp,1,myFun)
risk = as.vector(ifelse(riskScore > median(riskScore), "high", "low"))
clinical1 <- clinical
clinical1$riskScore <- riskScore
clinical1$risk <- risk
write.csv(clinical1, "output/clinical_lasso.csv")

##3. risk_map ####
library(pheatmap)
library(ggplotify)
library(ggplot2)
dat_risk <- clinical1[,c("OS", "OS.time", "riskScore", "risk")] %>% 
  dplyr::arrange(riskScore) %>%
  dplyr::mutate(x = 1:nrow(.)) %>% 
  dplyr::mutate(Status = dplyr::if_else(OS == 1, "Dead", 'Alive') %>% factor(levels = c("Dead", 'Alive')))

q1 <- ggplot(data = dat_risk, mapping = aes(x = x, y = riskScore, colour = risk)) + 
  #geom_bar(stat = 'identity', position = 'identity')
  geom_point(size = 2)+
  geom_vline(xintercept = nrow(dat_risk)/2, linetype = 2, color = "#566573")+
  geom_hline(yintercept = median(dat_risk$riskscore), linetype = 2, color = "#566573")+
  scale_color_manual(values = c("#EE0000", "#3B4992"))+
  theme_bw(base_size = 17)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(colour = "black"),
        legend.position = c(0.01,0.99),
        legend.justification = c(0, 1))+
  labs(x = "Patients(increasing risk score)")
q1
ggsave(q1, filename = "output/riskmap1.pdf", width = 6,height = 5)

q2 <- ggplot(data = dat_risk, mapping = aes(x = x, y = OS.time, colour = Status)) + 
  geom_point(size = 1.5)+
  scale_color_manual(values = c("#EE0000","#3B4992"))+
  geom_vline(xintercept = nrow(dat_risk)/2, linetype = 2, color = "#566573")+
  theme_bw(base_size = 17)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(colour = "black"),
        legend.position = c(0.01,0.99),
        legend.justification = c(0, 1))+
  labs(x="Patients(increasing risk score)")
q2
ggsave(q2, filename = "output/riskmap2.pdf", width = 6,height = 5)

##4. km ####
library(survival)
library(survminer)
rt <- clinical1
rt$OS.time <- rt$OS.time/365
diff <- survdiff(Surv(OS.time, OS) ~ risk, data = rt)
pValue <- 1 - pchisq(diff$chisq, df = 1)
if(pValue < 0.001){
  pValue <- "p < 0.001"
}else{
  pValue <- paste0("p = ", sprintf("%.03f", pValue))
}
fit <- survfit(Surv(OS.time, OS) ~ risk, data = rt)

# 
surPlot <- ggsurvplot(fit, 
                      data = rt,
                      conf.int = TRUE,
                      pval = pValue,
                      pval.size = 5,
                      legend.labs = c("high", "low"),
                      legend.title = "risk",
                      xlab = "Time(years)",
                      break.time.by = 1,
                      risk.table.title = "",
                      palette = c("red","blue"),
                      risk.table = T,
                      risk.table.height = .25)
pdf(file = "output/km.pdf", onefile = F, width = 6, height = 5)
print(surPlot)
dev.off()

##5. time_roc ####
library(timeROC)
time_roc_res <- timeROC(
  T = rt$OS.time * 365,
  delta = rt$OS,
  marker = rt$riskScore,
  cause = 1,
  times = c(1 * 365, 2 * 365, 3 * 365),
  ROC = T,
  iid = T
)
time_ROC_df <- data.frame(
  TP_1year = time_roc_res$TP[, 1],
  FP_1year = time_roc_res$FP[, 1],
  TP_2year = time_roc_res$TP[, 2],
  FP_2year = time_roc_res$FP[, 2],
  TP_3year = time_roc_res$TP[, 3],
  FP_3year = time_roc_res$FP[, 3]
)
p <- ggplot(data = time_ROC_df) +
  geom_line(aes(x = FP_1year, y = TP_1year), size = 1, color = "#3B4992") +
  geom_line(aes(x = FP_2year, y = TP_2year), size = 1, color = "#EE0000") +
  geom_line(aes(x = FP_3year, y = TP_3year), size = 1, color = "#008B45") +
  geom_abline(slope = 1, intercept = 0, color = "grey", size = 1, linetype = 2) +
  theme_bw(base_size = 17) +
  annotate("text",
           x = 0.75, y = 0.25, size = 4.5,
           label = paste0("AUC at 1 years = ", sprintf("%.3f", time_roc_res$AUC[[1]])), color = "#3B4992"
  ) +
  annotate("text",
           x = 0.75, y = 0.15, size = 4.5,
           label = paste0("AUC at 2 years = ", sprintf("%.3f", time_roc_res$AUC[[2]])), color = "#EE0000"
  ) +
  annotate("text",
           x = 0.75, y = 0.05, size = 4.5,
           label = paste0("AUC at 3 years = ", sprintf("%.3f", time_roc_res$AUC[[3]])), color = "#008B45"
  ) +
  labs(x = "1-specificity", y = "Sensitivity") +
  theme(
    panel.grid.minor = element_blank(), 
    panel.grid.major = element_blank(),
    axis.text = element_text(colour = "black")
    
  )
p
ggsave("output/timeroc.pdf", width = 5, height = 5)
