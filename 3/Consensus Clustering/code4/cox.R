rm(list = ls())
library(data.table)
library(tibble)
library(magrittr)
library(survival)

clinical <- fread("input/clinical_lasso.csv", data.table = F)
clinical <- clinical[,-1]
dat_cox <- clinical
dat_cox$OS.time <- dat_cox$OS.time/365
dat_cox$Stage[dat_cox$Stage == "Stage I"] <- 1
dat_cox$Stage[dat_cox$Stage == "Stage II"] <- 2
dat_cox$Stage[dat_cox$Stage == "Stage III"] <- 3
dat_cox$Stage[dat_cox$Stage == "Stage IV"] <- 4
dat_cox$Stage <- as.numeric(dat_cox$Stage)
str(dat_cox)
## uni_cox----
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
### 
ormulas <- sapply(dat_cox[,c("Age", "Stage", "riskScore")] %>% colnames(), 
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

### forest
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

## mul_cox----
fit <- as.formula(paste('Surv(', "OS.time", ",", "OS", ')~',  paste(res_uni1$Factor, collapse = "+"))) %>% 
  coxph(data = dat_cox)
res_mul <- fit %>%
  cox_extr
res_mul1 <- res_mul %>% 
  dplyr::mutate(Factor = rownames(.), HR = paste(res_mul$mean, " [", res_mul$lower, " - ", res_mul$upper, "]", sep = "")) %>%
  dplyr::select("Factor", "p.value", "HR", "mean", "lower", "upper")
write.csv(res_mul1, "output/res_mul.csv", row.names = F)

### forest
labeltext <- rbind(colnames(res_mul1)[1:3], res_mul1[, 1:3])
HR <- rbind(rep(NA, 3), res_mul1[, 4:6])
pdf("output/forest_mul.pdf", height = 6, width = 10, onefile = F)
forestplot(labeltext = labeltext, 
           HR,
           zero = 1, 
           lwd.ci = 1,
           title ="multifactor Cox regression analysis",
           colgap = unit(2, 'mm'),
           xlab = "HR", 
           txt_gp=fpTxtGp(ticks = gpar(cex = 0.8), xlab = gpar(cex = 0.8),
                          title = gpar(cex = 1.2), cex = 1) ,
           boxsize = 0.2, 
           graph.pos = 4,
           xlog = T,
           col = fpColors(box = "#3B4992", line = "black", zero = "black"))
dev.off()
