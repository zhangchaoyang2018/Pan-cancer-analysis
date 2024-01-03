
############################## 2 KEGG pathways to 130 nucleotide names/KEGG ids and HMDB ids and mass values ######################
setwd("C:/MyRdata8/Summary")
nucleotide <- read.csv("Nucleotide metabolism1.csv")
head(nucleotide)
dim(nucleotide)
colnames(nucleotide)[6] <- "compound_id"
head(nucleotide)
dim(nucleotide)

setwd("C:/Database1 HMDB annotation")
library(dplyr)
annotation <- read.csv("355372 Annotation on all mass values.csv")
dim(annotation)
head(annotation)
annotation <- annotation[,-c(1,4)]
dim(annotation)
head(annotation)

nucleotide.annotation <- merge(nucleotide,annotation)
dim(nucleotide.annotation)
head(nucleotide.annotation)


############################## Lung primary tumor (n=85) #####################################
library(survival)
library(survminer)

setwd("C:/MyRdata8/Lung primary tumor (n=85)")

table <- read.csv("Lung primary tumor (n=85).csv")
dim(table)
table[1:5,1:7]

table1 <- table[,c(2,3,4,7,22,42,44,50:10024)]
dim(table1)
table1[1:5,1:10]
rownames(table1) <- table1[,1]
table2 <- table1[,-1]
dim(table2)
table2[1:5,1:10]
colnames(table2)[2] <- 'Sex'

table3 <- table2[,-c(1:6)]
dim(table3)
table3[1:5,1:7]
table3.n85 <- table3

# nucleotide profile data
table3 <- table3[,match(intersect(round(as.numeric(substr(colnames(table3),2,12)),4),nucleotide.annotation$query_mass),round(as.numeric(substr(colnames(table3),2,12)),4))]
dim(table3)
table3[1:5,1:7]

# median intensity on each mass
median1 <- apply(table3,MARGIN=2,FUN=median)
head(median1)
length(median1)
median1.n85 <- median1


# save nucleotide profile data
nucleotide.profile.lung.primary.tumor.n85 <- table3
dim(nucleotide.profile.lung.primary.tumor.n85)


# calculate zero number for each m/z value and cutoff
len1 <- c()
for (i in 1:ncol(table3)) {
  len <- length(which(table3[,i]==0))
  len1 <- c(len1,len)
}
len1
table4 <- table3[,which(len1 <= 85)]
dim(table4)

table2[1:5,1:10]
table5 <- cbind(table2[,5:6],table4)
dim(table5)
table5[1:5,1:7]

res.cut <- surv_cutpoint(table5, time = "OS", event = "STATOS",
                         minprop = 0.1,
                         variables = colnames(table5)[3:ncol(table5)])
res.cat.1 <- surv_categorize(res.cut,labels = c(0, 1))
dim(res.cat.1)
res.cat.1[1:5,1:7]

table2[1:5,1:10]
res.cat <- cbind(table2[,c(1,2,4)],res.cat.1)
dim(res.cat)
res.cat[1:5,1:10]
res.cat$AGE <- ifelse(res.cat$AGE<60,0,1)
res.cat[1:5,1:10]

# save profile data with cutoff
res.cat.lung.primary.tumor.n85 <- res.cat
dim(res.cat.lung.primary.tumor.n85)
head(res.cat.lung.primary.tumor.n85)


################ univariate cox analysis to calculate p value, HR
result <- data.frame()
for (i in c(1:3,6:ncol(res.cat))) {
  res.cox <- coxph(Surv(OS,STATOS)~res.cat[,i],data =res.cat)
  summ <- summary(res.cox)
  hypo <- cox.zph(res.cox)
  hypo1 <- hypo$table[-nrow(hypo$table),]
  result1 <- data.frame(mz.value=colnames(res.cat)[i],
                        HR = signif(summ$coef[2], digits=3),
                        HR.confint.lower = signif(summ$conf.int[,"lower .95"], 3),
                        HR.confint.upper = signif(summ$conf.int[,"upper .95"],3),
                        logtest.p=summ$logtest[3],
                        waldtest.p=summ$waldtest[3],
                        sctest.p=summ$sctest[3],
                        hypo.p=hypo1[3])
  result <- rbind(result,result1)
}
head(result)
dim(result)
univariate.lung.primary.tumor.n85 <- result


# screen by sctest p value
result2 <- result[result$sctest.p<0.05&result$hypo.p>0.05,]
dim(result2)
head(result2)
result3 <- result2[-c(1:2),]
head(result3)
dim(result3)

# annotation on significant mass
result3$query_mass <- round(as.numeric(substr(result3$mz.value,2,12)),4)
result4 <- merge(result3,nucleotide.annotation)
result4$name.adduct <- paste(result4$adduct,result4$name,sep = '.')
result4$median.intensity <- median1[match(result4$mz.value,colnames(table3))]
dim(result4)
head(result4)
result5 <- result4[order(result4$ppm,decreasing = F),]
result6 <- result5[which(duplicated(result5$name.adduct)==F),]
dim(result6)
result7 <- result6[order(result6$median.intensity,decreasing = T),]
result8 <- result7[which(duplicated(result7$name)==F),]
result8 <- result8[order(result8$query_mass,decreasing = F),]
dim(result8)
head(result8)
univariate4.lung.primary.tumor.n85 <- result8


############## validate by survival curve

# fit survival curve
fit <- survfit(Surv(OS,STATOS)~Sex,data =res.cat)
fit
survdiff(Surv(OS,STATOS)~Sex ,data =res.cat) #采用log-rank 检验分析生存率差异

##绘制生存曲线##
ggsurvplot(fit, data = res.cat, conf.int = F, # 增加置信区间
           risk.table = TRUE, pval=T) #不同时间点风险人数表，增加P值



################# forest plot
library(forestplot)

dim(result8)
head(result8)
length(unique(result8$query_mass))
head(result2)

result10 <- data.frame(Feature=c(result2$mz.value[1:2],result8$query_mass),
                       Annotation=c("Male:2","3/4",result8$name),
                       HR=c(result2$HR[1:2],result8$HR),
                       HR.lower=c(result2$HR.confint.lower[1:2],result8$HR.confint.lower),
                       HR.upper=c(result2$HR.confint.upper[1:2],result8$HR.confint.upper),
                       pvalue=round(c(result2$sctest.p[1:2],result8$sctest.p),3))
dim(result10)
head(result10)

result10$HR1 <- paste(result10$HR,"(",result10$HR.lower,"-",result10$HR.upper,")",sep = '')
dim(result10)
head(result10)

labeltext <- as.matrix(rbind(c("Feature","Annotation","HR","pvalue"),result10[,c(1,2,7,6)]))
labeltext

forestplot(labeltext, 
           mean = c(NA,result10$HR),
           lower = c(NA,result10$HR.lower),
           upper = c(NA,result10$HR.upper),
           zero = 1,lwd.zero = 2,
           clip = c(0,6),
           graph.pos = 3,
           boxsize = 0.5,
           align = "l",
           xlog = FALSE,
           col=fpColors(line = "black", zero = "black",
                        box="black"),
           lty.ci = 7,
           lwd.ci = 1,
           ci.vertices.height = 0.15,
           lwd.xaxis=2,
           xlab="Hazard Ratio",
           txt_gp = fpTxtGp(ticks = gpar(cex = 0.9), xlab = gpar(cex = 0.8), cex = 0.7),
           lineheight = unit(8,"mm")) 




################################################ multivariate cox analysis identify independent prognostic factor
dim(res.cat)
res.cat[1:5,1:10]

head(result2)
head(result8)

res.cat1 <- res.cat[,c("OS","STATOS",result2$mz.value[1:2],unique(result8$mz.value))]
dim(res.cat1)
res.cat1[1:5,1:7]

# save final profiling data after annotation
res.cat1.lung.primary.tumor.n85 <- res.cat1


################## method1 feature selection from the beginning

### Step regression on Cox model
res.cox1 <- coxph(Surv(OS, STATOS) ~ .,data =res.cat1)
res.cox2 <- step(res.cox1)
res.cox3 <- coxph(Surv(OS, STATOS) ~ UICC + X211.0014 + X225.0493 + X228.0301 + 
                    X272.0778 + X278.0573 + X285.0519 + X287.0062 + X289.0333 + 
                    X289.0526 + X299.0647 + X313.0329 + X322.0447 + X362.0506 + 
                    X365.0509 + X402.011 + X424.918 + X424.977 + X462.9862 + 
                    X505.9817,data =res.cat1)
ggforest(model=res.cox3,data=res.cat1)


###等比例风险假定
hypo <- cox.zph(res.cox3) 
#hypo <- cox.zph(res.cox4) 
hypo1 <- hypo$table[-nrow(hypo$table),]
dim(hypo1)
hypo2 <- hypo1[hypo1[,3]>0.05,]
dim(hypo2)
hypo2

res.cat2 <- res.cat1[,c("OS","STATOS",rownames(hypo2))]
dim(res.cat2)
res.cat2[1:5,1:7]
res.cox4 <- coxph(Surv(OS, STATOS) ~ .,data = res.cat2)
ggforest(model=res.cox4,data=res.cat2)


summ <- as.data.frame(cbind(summary(res.cox4)$coefficients[,c(2,5)],summary(res.cox4)$conf.int[,c(3,4)]))
summ$mz.value <- rownames(summ)
head(summ)
dim(summ)
multivariate.lung.primary.tumor.n85.sum <- summ

summ1 <- summ[-c(1),]
dim(summ1)
head(summ1)

summ2 <- merge(summ1,result8)
multivariate.lung.primary.tumor.n85 <- summ2
dim(multivariate.lung.primary.tumor.n85)
head(multivariate.lung.primary.tumor.n85)
multivariate.lung.primary.tumor.n85.sig <- summ2[summ2[,3]<0.05,]
dim(multivariate.lung.primary.tumor.n85.sig)


#### forest plot
library(forestplot)

result12 <- summ2
dim(result12)
head(result12)

head(summ)
result13 <- summ[1,]
head(result13)

result14 <- data.frame(Feature=c(result13$mz.value,result12$query_mass),
                       Annotation=c("3/4",result12$name),
                       HR=c(result13[,1],result12[,2]),
                       HR.lower=c(result13[,3],result12[,4]),
                       HR.upper=c(result13[,4],result12[,5]),
                       pvalue=round(c(result13[,2],result12[,3]),3))
dim(result14)
head(result14)

result14$HR1 <- paste(round(result14$HR,2),"(",round(result14$HR.lower,2),"-",round(result14$HR.upper,2),")",sep = '')
result14$pvalue1 <- ifelse(result14$pvalue<0.001,"***",ifelse(result14$pvalue<0.01,"**",ifelse(result14$pvalue<0.05,"*"," ")))
result14$pvalue2 <- paste(result14$pvalue,result14$pvalue1,sep = '')
dim(result14)
head(result14)

labeltext <- as.matrix(rbind(c("Feature","Annotation","HR","pvalue"),result14[,c(1,2,7,9)]))
labeltext

forestplot(labeltext, 
           mean = c(NA,result14$HR),
           lower = c(NA,result14$HR.lower),
           upper = c(NA,result14$HR.upper),
           zero = 1,lwd.zero = 2,
           clip = c(0.01,40),
           graph.pos = 3,
           boxsize = 0.5,
           align = "l",
           xlog = T,
           col=fpColors(line = "black", zero = "black",
                        box="black"),
           lty.ci = 7,
           lwd.ci = 1,
           ci.vertices.height = 0.15,
           lwd.xaxis=2,
           xlab="Hazard Ratio",
           txt_gp = fpTxtGp(ticks = gpar(cex = 0.9), xlab = gpar(cex = 0.8), cex = 0.8),
           lineheight = unit(8,"mm")) 



############### use Cox model to predict survival outcome by ROC curve
res.cat3 <- res.cat1[,c("OS","STATOS",rownames(summ))]
dim(res.cat3)
res.cat3[1:5,1:5]

res.cox3 <- coxph(Surv(OS,STATOS)~.,data =res.cat3)
sum2 <- summary(res.cox3)
# C index export
sum2$concordance

# forest plot
ggforest(model = res.cox3,data=res.cat3)

library(pROC)
gfit <- roc(STATOS ~predict(res.cox3), data = res.cat3)

plot(gfit,
     print.auc=TRUE, #输出AUC值
     print.thres=TRUE, #输出cut-off值
     #main = "ROC Curve in Lung NAC tumor (n=77)", #设置图形的标题
     col= "red", #曲线颜色
     print.thres.col="black", #cut-off值字体的颜色
     identity.col="blue", #对角线颜色
     identity.lty=1,identity.lwd=1)


################# use Cox model to predict good/poor prognosis or low/high risk groups by ROC curve
res.cat3$STATOS1 <- ifelse(res.cat3$OS < median(res.cat3$OS),1,0) #1=short/poor survival or high risk,0=long/good survival or low risk  
res.cat3$STATOS1

gfit1 <- roc(STATOS1 ~predict(res.cox3), data = res.cat3)

cutoff <- gfit1$thresholds[which.max(gfit1$sensitivities+gfit1$specificities)]
cutoff
plot(gfit1,
     print.auc=TRUE, #输出AUC值
     print.thres=TRUE, #输出cut-off值
     #main = "ROC Curve in Lung NAC tumor (n=77)", #设置图形的标题
     col= "red", #曲线颜色
     print.thres.col="black", #cut-off值字体的颜色
     identity.col="blue", #对角线颜色
     identity.lty=1,identity.lwd=1)

### cutoff-based cox model to stratify patients into different prognostic risk groups by survival curve
res.cat3[1:5,1:ncol(res.cat3)]
res.cat3$risk <- ifelse(predict(res.cox3) < cutoff,'low','high')
# fit survival curve
fit <- survfit(Surv(OS,STATOS)~risk,data =res.cat3)
##绘制生存曲线##
ggsurvplot(fit, data = res.cat3, conf.int = F, # 增加置信区间
           risk.table = TRUE, pval=T) #不同时间点风险人数表，增加P值


################# Nomogram绘制

library(car)
library(rms)
library(pROC)
library(timeROC)
library(ggDCA)

res.cat4 <- res.cat1[,c("OS","STATOS",rownames(summ))]
dim(res.cat4)
res.cat4[1:5,1:10]

res.cox4 <- cph(as.formula(paste0('Surv(OS,STATOS)~',paste(colnames(res.cat4)[-c(1:2)],sep = '',collapse = '+'))),
                x=T,y=T,surv=T,data =res.cat4)
print(res.cox4)

ddist <- datadist(res.cat4)
options(datadist='ddist')

surv <- Survival(res.cox4) # 建立生存函数
surv

surv1 <- function(x) surv(12,x) # 1年OS
surv2 <- function(x) surv(36,x) # 3年OS
surv3 <- function(x) surv(60,x) # 5年OS

nom <- nomogram(res.cox4,
                fun = list(surv1,surv2,surv3),
                funlabel = c('1-year survival Probability',
                             '3-year survival Probability',
                             '5-year survival Probability'))
plot(nom, 
     lplabel="Linear Predictor",
     xfrac = 0.2, # 左侧标签距离坐标轴的距离
     #varname.label = TRUE, 
     tcl = -0.2, # 刻度长短和方向 
     lmgp = 0.1, # 坐标轴标签距离坐标轴远近
     points.label ='Points', 
     total.points.label = 'Total Points',
     cap.labels = FALSE,
     cex.var = 0.8, # 左侧标签字体大小
     cex.axis = 0.7, # 坐标轴字体大小
     col.grid = gray(c(0.8, 0.95))) # 竖线颜色


############# ROC curve on the prediction of 1/3/5-year survival rate
library(timeROC)
res.cat4 <- res.cat1[,c("OS","STATOS",rownames(summ))]
dim(res.cat4)
res.cat4[1:5,1:5]

res.cox4 <- cph(as.formula(paste0('Surv(OS,STATOS)~',paste(colnames(res.cat4)[-c(1:2)],sep = '',collapse = '+'))),
                x=T,y=T,surv=T,data =res.cat4)

pred_f_training <- predict(res.cox4,res.cat4,type="lp")#!!!type="lp"
data_table <- data.frame(time=res.cat4[,"OS"],status=res.cat4[,"STATOS"],score=pred_f_training)


time_roc_res <- timeROC(
  T = data_table$time,
  delta = data_table$status,
  marker = data_table$score,
  cause = 1,
  weighting="marginal",
  times = c(12, 36, 60),
  ROC = TRUE,
  iid = TRUE
)


time_ROC_df <- data.frame(
  TP_1year = time_roc_res$TP[, 1],
  FP_1year = time_roc_res$FP[, 1],
  TP_3year = time_roc_res$TP[, 2],
  FP_3year = time_roc_res$FP[, 2],
  TP_5year = time_roc_res$TP[, 3],
  FP_5year = time_roc_res$FP[, 3]
)


ggplot(data = time_ROC_df) +
  geom_line(aes(x = FP_1year, y = TP_1year), linewidth = 1, color = "#BC3C29FF") +
  geom_line(aes(x = FP_3year, y = TP_3year), linewidth = 1, color = "#0072B5FF") +
  geom_line(aes(x = FP_5year, y = TP_5year), linewidth = 1, color = "#E18727FF") +
  geom_abline(slope = 1, intercept = 0, color = "grey", linewidth = 1, linetype = 2) +
  theme_bw() +
  annotate("text",
           x = 0.75, y = 0.25, size = 4.5,
           label = paste0("AUC at 1 year = ", sprintf("%.3f", time_roc_res$AUC[[1]])), color = "#BC3C29FF"
  ) +
  annotate("text",
           x = 0.75, y = 0.15, size = 4.5,
           label = paste0("AUC at 3 years = ", sprintf("%.3f", time_roc_res$AUC[[2]])), color = "#0072B5FF"
  ) +
  annotate("text",
           x = 0.75, y = 0.05, size = 4.5,
           label = paste0("AUC at 5 years = ", sprintf("%.3f", time_roc_res$AUC[[3]])), color = "#E18727FF"
  ) +
  labs(x = "Specificity", y = "Sensitivity") +
  theme(
    axis.text = element_text(size = 10, color = "black"),
    axis.title.x = element_text(size = 12, color = "black", margin = margin(c(15, 0, 0, 0))),
    axis.title.y = element_text(size = 12, color = "black", margin = margin(c(0, 15, 0, 0)))
  )

prognosis.n85 <- data.frame(Tumor=rep("Lung Primary Tumor",5),
                            Prognosis=c("Survival outcome","Risk stratification","1-year survival rate","3-year survival rate","5-year survival rate"),
                            AUC=c(gfit$auc,gfit1$auc,time_roc_res$AUC[[1]],time_roc_res$AUC[[2]],time_roc_res$AUC[[3]]))
prognosis.n85


############################## Lung NAC tumor (n=77) #####################################
library(survival)
library(survminer)

setwd("C:/MyRdata8/Lung NAC tumor (n=77)")

table <- read.csv("Lung NAC tumor (n=77).csv")
dim(table)
table[1:5,1:10]

table1 <- table[,c(2,3,4,10,22,42,44,50:10024)]
dim(table1)
table1[1:5,1:10]
rownames(table1) <- table1[,1]
table2 <- table1[,-1]
dim(table2)
table2[1:5,1:10]
colnames(table2)[2] <- 'Sex'
clin.n77 <- table2
table2[1:5,1:10]

table3 <- table2[,-c(1:6)]
dim(table3)
table3[1:5,1:7]
table3.n77 <- table3

# nucleotide profile data
table3 <- table3[,match(intersect(round(as.numeric(substr(colnames(table3),2,12)),4),nucleotide.annotation$query_mass),round(as.numeric(substr(colnames(table3),2,12)),4))]
dim(table3)

# median intensity on each mass
median1 <- apply(table3,MARGIN=2,FUN=median)
head(median1)
length(median1)
median1.n77 <- median1

# save nucleotide profile data
nucleotide.profile.lung.NAC.tumor.n77 <- table3
dim(nucleotide.profile.lung.NAC.tumor.n77)


# calculate zero number for each m/z value
len1 <- c()
for (i in 1:ncol(table3)) {
  len <- length(which(table3[,i]==0))
  len1 <- c(len1,len)
}
len1
table4 <- table3[,which(len1 <= 77)]
dim(table4)

table2[1:5,1:10]
table5 <- cbind(table2[,5:6],table4)
dim(table5)
table5[1:5,1:10]

res.cut <- surv_cutpoint(table5, time = "OS", event = "STATOS",
                         minprop = 0.1,
                         variables = colnames(table5)[3:ncol(table5)])
res.cat.1 <- surv_categorize(res.cut,labels = c(0, 1))
dim(res.cat.1)
res.cat.1[1:5,1:7]

table2[1:5,1:10]
res.cat <- cbind(table2[,1:4],res.cat.1)
dim(res.cat)
res.cat[1:5,1:10]
res.cat$Tumor.Regression <- ifelse(res.cat$Tumor.Regression=="Non-responder",0,1)
res.cat[1:5,1:10]

# save profile data with cutoff
res.cat.lung.NAC.tumor.n77 <- res.cat
dim(res.cat.lung.NAC.tumor.n77)
head(res.cat.lung.NAC.tumor.n77)


################ univariate cox analysis to calculate p value, HR
result <- data.frame()
for (i in c(1:4,7:ncol(res.cat))) {
  res.cox <- coxph(Surv(OS,STATOS)~res.cat[,i],data =res.cat)
  summ <- summary(res.cox)
  hypo <- cox.zph(res.cox)
  hypo1 <- hypo$table[-nrow(hypo$table),]
  result1 <- data.frame(mz.value=colnames(res.cat)[i],
                        HR = signif(summ$coef[2], digits=3),
                        HR.confint.lower = signif(summ$conf.int[,"lower .95"], 3),
                        HR.confint.upper = signif(summ$conf.int[,"upper .95"],3),
                        logtest.p=summ$logtest[3],
                        waldtest.p=summ$waldtest[3],
                        sctest.p=summ$sctest[3],
                        hypo.p=hypo1[3])
  result <- rbind(result,result1)
}
head(result)
dim(result)
univariate.lung.NAC.tumor.n77 <- result


# screen by sctest p value
result2 <- result[result$sctest.p<0.05&result$hypo.p>0.05,]
dim(result2)
head(result2)
result3 <- result2[-c(1:2),]
head(result3)
dim(result3)

# annotation on significant mass
result3$query_mass <- round(as.numeric(substr(result3$mz.value,2,12)),4)
result4 <- merge(result3,nucleotide.annotation)
result4$name.adduct <- paste(result4$adduct,result4$name,sep = '.')
result4$median.intensity <- median1[match(result4$mz.value,colnames(table3))]
dim(result4)
result5 <- result4[order(result4$ppm,decreasing = F),]
result6 <- result5[which(duplicated(result5$name.adduct)==F),]
dim(result6)
result7 <- result6[order(result6$median.intensity,decreasing = T),]
result8 <- result7[which(duplicated(result7$name)==F),]
result8 <- result8[order(result8$query_mass,decreasing = F),]
dim(result8)
head(result8)
univariate4.lung.NAC.tumor.n77 <- result8


############## validate by survival curve

# fit survival curve
fit <- survfit(Surv(OS,STATOS)~UICC,data =res.cat)
fit
survdiff(Surv(OS,STATOS)~UICC,data =res.cat) #采用log-rank 检验分析生存率差异

##绘制生存曲线##
ggsurvplot(fit, data = res.cat, conf.int = F, # 增加置信区间
           risk.table = TRUE, pval=T) #不同时间点风险人数表，增加P值


################# forest plot
library(forestplot)

dim(result8)
head(result8)
length(unique(result8$query_mass))
head(result2)

result10 <- data.frame(Feature=c(result2$mz.value[1:2],result8$query_mass),
                       Annotation=c("Responder","0/1/2/3/4",result8$name),
                       HR=c(result2$HR[1:2],result8$HR),
                       HR.lower=c(result2$HR.confint.lower[1:2],result8$HR.confint.lower),
                       HR.upper=c(result2$HR.confint.upper[1:2],result8$HR.confint.upper),
                       pvalue=round(c(result2$sctest.p[1:2],result8$sctest.p),3))
dim(result10)
head(result10)

result10$HR1 <- paste(result10$HR,"(",result10$HR.lower,"-",result10$HR.upper,")",sep = '')
dim(result10)
head(result10)

labeltext <- as.matrix(rbind(c("Feature","Annotation","HR","pvalue"),result10[,c(1,2,7,6)]))
labeltext

forestplot(labeltext, 
           mean = c(NA,result10$HR),
           lower = c(NA,result10$HR.lower),
           upper = c(NA,result10$HR.upper),
           zero = 1,lwd.zero = 2,
           clip = c(0,6),
           graph.pos = 3,
           boxsize = 0.5,
           align = "l",
           xlog = FALSE,
           col=fpColors(line = "black", zero = "black",
                        box="black"),
           lty.ci = 7,
           lwd.ci = 1,
           ci.vertices.height = 0.15,
           lwd.xaxis=2,
           xlab="Hazard Ratio",
           txt_gp = fpTxtGp(ticks = gpar(cex = 0.9), xlab = gpar(cex = 0.8), cex = 0.7),
           lineheight = unit(8,"mm")) 



################################################ multivariate cox analysis identify independent prognostic factor

dim(res.cat)
res.cat[1:5,1:10]

head(result2)
head(result8)

res.cat1 <- res.cat[,c("OS","STATOS",result2$mz.value[1:2],unique(result8$mz.value))]
dim(res.cat1)
res.cat1[1:5,1:7]

# save final profiling data after annotation
res.cat1.lung.primary.tumor.n85 <- res.cat1


################## feature selection from the beginning

### Step regression on Cox model
res.cox1 <- coxph(Surv(OS, STATOS) ~ .,data =res.cat1)
res.cox2 <- step(res.cox1)
res.cox3 <- coxph(Surv(OS, STATOS) ~ UICC + X226.0833 + X282.0844 + X287.0062 + 
                    X299.0647 + X304.0338 + X348.0341 + X365.0509 + X368.0347 + 
                    X447.991,data =res.cat1)
ggforest(model=res.cox3,data=res.cat1)


###等比例风险假定
hypo <- cox.zph(res.cox3) 
hypo1 <- hypo$table[-nrow(hypo$table),]
dim(hypo1)
hypo2 <- hypo1[hypo1[,3]>0.05,]
dim(hypo2)
hypo2

res.cat2 <- res.cat1[,c("OS","STATOS",rownames(hypo2))]
dim(res.cat2)
res.cat2[1:5,1:7]
res.cox4 <- coxph(Surv(OS, STATOS) ~ .,data = res.cat2)
ggforest(model=res.cox4,data=res.cat2)


summ <- as.data.frame(cbind(summary(res.cox4)$coefficients[,c(2,5)],summary(res.cox4)$conf.int[,c(3,4)]))
summ$mz.value <- rownames(summ)
head(summ)
dim(summ)
multivariate.lung.NAC.tumor.n77.sum <- summ

summ1 <- summ[-c(1),]
dim(summ1)
head(summ1)

summ2 <- merge(summ1,result8)
multivariate.lung.NAC.tumor.n77 <- summ2
dim(multivariate.lung.NAC.tumor.n77)
head(multivariate.lung.NAC.tumor.n77)
multivariate.lung.NAC.tumor.n77.sig <- summ2[summ2[,3]<0.05,]
dim(multivariate.lung.NAC.tumor.n77.sig)


#### forest plot
library(forestplot)

result12 <- summ2
dim(result12)
head(result12)

head(summ)
result13 <- summ[1,]
head(result13)

result14 <- data.frame(Feature=c(result13$mz.value,result12$query_mass),
                       Annotation=c("0/1/2/3/4",result12$name),
                       HR=c(result13[,1],result12[,2]),
                       HR.lower=c(result13[,3],result12[,4]),
                       HR.upper=c(result13[,4],result12[,5]),
                       pvalue=round(c(result13[,2],result12[,3]),3))
dim(result14)
head(result14)

result14$HR1 <- paste(round(result14$HR,2),"(",round(result14$HR.lower,2),"-",round(result14$HR.upper,2),")",sep = '')
result14$pvalue1 <- ifelse(result14$pvalue<0.001,"***",ifelse(result14$pvalue<0.01,"**",ifelse(result14$pvalue<0.05,"*"," ")))
result14$pvalue2 <- paste(result14$pvalue,result14$pvalue1,sep = '')
dim(result14)
head(result14)

labeltext <- as.matrix(rbind(c("Feature","Annotation","HR","pvalue"),result14[,c(1,2,7,9)]))
labeltext

forestplot(labeltext, 
           mean = c(NA,result14$HR),
           lower = c(NA,result14$HR.lower),
           upper = c(NA,result14$HR.upper),
           zero = 1,lwd.zero = 2,
           clip = c(0.01,40),
           graph.pos = 3,
           boxsize = 0.5,
           align = "l",
           xlog = T,
           col=fpColors(line = "black", zero = "black",
                        box="black"),
           lty.ci = 7,
           lwd.ci = 1,
           ci.vertices.height = 0.15,
           lwd.xaxis=2,
           xlab="Hazard Ratio",
           txt_gp = fpTxtGp(ticks = gpar(cex = 0.9), xlab = gpar(cex = 0.8), cex = 0.8),
           lineheight = unit(8,"mm")) 





############### use Cox model to predict survival outcome by ROC curve
res.cat3 <- res.cat[,c("OS","STATOS",rownames(summ))]
dim(res.cat3)
res.cat3[1:5,1:5]

res.cox3 <- coxph(Surv(OS,STATOS)~.,data =res.cat3)
sum2 <- summary(res.cox3)
# C index export
sum2$concordance

# forest plot
ggforest(model = res.cox3,data=res.cat3)

library(pROC)
gfit <- roc(STATOS ~predict(res.cox3), data = res.cat3)

plot(gfit,
     print.auc=TRUE, #输出AUC值
     print.thres=TRUE, #输出cut-off值
     #main = "ROC Curve in Lung NAC tumor (n=77)", #设置图形的标题
     col= "red", #曲线颜色
     print.thres.col="black", #cut-off值字体的颜色
     identity.col="blue", #对角线颜色
     identity.lty=1,identity.lwd=1)


################# use Cox model to predict good/poor prognosis or low/high risk groups by ROC curve
res.cat3$STATOS1 <- ifelse(res.cat3$OS < median(res.cat3$OS),0,1) #0=short/poor survival or high risk,1=long/good survival or low risk  
res.cat3$STATOS1

gfit1 <- roc(STATOS1 ~predict(res.cox3), data = res.cat3)
cutoff <- gfit1$thresholds[which.max(gfit1$sensitivities+gfit1$specificities)]
cutoff
plot(gfit1,
     print.auc=TRUE, #输出AUC值
     print.thres=TRUE, #输出cut-off值
     #main = "ROC Curve in Lung NAC tumor (n=77)", #设置图形的标题
     col= "red", #曲线颜色
     print.thres.col="black", #cut-off值字体的颜色
     identity.col="blue", #对角线颜色
     identity.lty=1,identity.lwd=1)

### cutoff-based cox model to stratify patients into different prognostic risk groups by survival curve
res.cat3[1:5,1:ncol(res.cat3)]
res.cat3$risk <- ifelse(predict(res.cox3) < cutoff,'low','high')
# fit survival curve
fit <- survfit(Surv(OS,STATOS)~risk,data =res.cat3)
##绘制生存曲线##
ggsurvplot(fit, data = res.cat3, conf.int = F, # 增加置信区间
           risk.table = TRUE, pval=T) #不同时间点风险人数表，增加P值



################# Nomogram绘制

library(car)
library(rms)
library(pROC)
library(timeROC)
library(ggDCA)

res.cat4 <- res.cat1[,c("OS","STATOS",rownames(summ))]
dim(res.cat4)
res.cat4[1:5,1:10]

res.cox4 <- cph(as.formula(paste0('Surv(OS,STATOS)~',paste(colnames(res.cat4)[-c(1:2)],sep = '',collapse = '+'))),
                x=T,y=T,surv=T,data =res.cat4)
print(res.cox4)

ddist <- datadist(res.cat4)
options(datadist='ddist')

surv <- Survival(res.cox4) # 建立生存函数
surv

surv1 <- function(x) surv(12,x) # 1年OS
surv2 <- function(x) surv(36,x) # 3年OS
surv3 <- function(x) surv(60,x) # 5年OS

nom <- nomogram(res.cox4,
                fun = list(surv1,surv2,surv3),
                funlabel = c('1-year survival Probability',
                             '3-year survival Probability',
                             '5-year survival Probability'))
plot(nom, 
     lplabel="Linear Predictor",
     xfrac = 0.2, # 左侧标签距离坐标轴的距离
     #varname.label = TRUE, 
     tcl = -0.2, # 刻度长短和方向 
     lmgp = 0.1, # 坐标轴标签距离坐标轴远近
     points.label ='Points', 
     total.points.label = 'Total Points',
     cap.labels = FALSE,
     cex.var = 0.8, # 左侧标签字体大小
     cex.axis = 0.7, # 坐标轴字体大小
     col.grid = gray(c(0.8, 0.95))) # 竖线颜色


############# ROC curve on the prediction of 1/3/5-year survival rate
library(timeROC)
res.cat4 <- res.cat1[,c("OS","STATOS",rownames(summ))]
dim(res.cat4)
res.cat4[1:5,1:5]

res.cox4 <- cph(as.formula(paste0('Surv(OS,STATOS)~',paste(colnames(res.cat4)[-c(1:2)],sep = '',collapse = '+'))),
                x=T,y=T,surv=T,data =res.cat4)
print(res.cox4)

pred_f_training <- predict(res.cox4,res.cat4,type="lp")#!!!type="lp",是他没错
data_table <- data.frame(time=res.cat4[,"OS"],status=res.cat4[,"STATOS"],score=pred_f_training)


time_roc_res <- timeROC(
  T = data_table$time,
  delta = data_table$status,
  marker = data_table$score,
  cause = 1,
  weighting="marginal",
  times = c(12, 36, 60),
  ROC = TRUE,
  iid = TRUE
)


time_ROC_df <- data.frame(
  TP_1year = time_roc_res$TP[, 1],
  FP_1year = time_roc_res$FP[, 1],
  TP_3year = time_roc_res$TP[, 2],
  FP_3year = time_roc_res$FP[, 2],
  TP_5year = time_roc_res$TP[, 3],
  FP_5year = time_roc_res$FP[, 3]
)


ggplot(data = time_ROC_df) +
  geom_line(aes(x = FP_1year, y = TP_1year), linewidth = 1, color = "#BC3C29FF") +
  geom_line(aes(x = FP_3year, y = TP_3year), linewidth = 1, color = "#0072B5FF") +
  geom_line(aes(x = FP_5year, y = TP_5year), linewidth = 1, color = "#E18727FF") +
  geom_abline(slope = 1, intercept = 0, color = "grey", linewidth = 1, linetype = 2) +
  theme_bw() +
  annotate("text",
           x = 0.75, y = 0.25, size = 4.5,
           label = paste0("AUC at 1 year = ", sprintf("%.3f", time_roc_res$AUC[[1]])), color = "#BC3C29FF"
  ) +
  annotate("text",
           x = 0.75, y = 0.15, size = 4.5,
           label = paste0("AUC at 3 years = ", sprintf("%.3f", time_roc_res$AUC[[2]])), color = "#0072B5FF"
  ) +
  annotate("text",
           x = 0.75, y = 0.05, size = 4.5,
           label = paste0("AUC at 5 years = ", sprintf("%.3f", time_roc_res$AUC[[3]])), color = "#E18727FF"
  ) +
  labs(x = "Specificity", y = "Sensitivity") +
  theme(
    axis.text = element_text(size = 10, color = "black"),
    axis.title.x = element_text(size = 12, color = "black", margin = margin(c(15, 0, 0, 0))),
    axis.title.y = element_text(size = 12, color = "black", margin = margin(c(0, 15, 0, 0)))
  )


prognosis.n77 <- data.frame(Tumor=rep("Lung NAC Tumor",5),
                            Prognosis=c("Survival outcome","Risk stratification","1-year survival rate","3-year survival rate","5-year survival rate"),
                            AUC=c(gfit$auc,gfit1$auc,time_roc_res$AUC[[1]],time_roc_res$AUC[[2]],time_roc_res$AUC[[3]]))
prognosis.n77


############################## Lunge Plattenepithel (n=238) #####################
library(survival)
library(survminer)

setwd("C:/MyRdata8/Lunge Plattenepithel (n=238)")
table <- read.csv("Lunge Plattenepithel (n=238).csv")
dim(table)
table[1:5,1:10]

table1 <- table[,c(1,2,3,5,6,21,53:11414)]
dim(table1)
table1[1:5,1:10]
rownames(table1) <- table1[,1]
table2 <- table1[,-1]
dim(table2)
table2[1:5,1:7]

table3 <- table2[,-c(1:5)]
dim(table3)
table3[1:5,1:7]
table3.n238 <- table3

# nucleotide profile data
table3 <- table3[,match(intersect(round(as.numeric(substr(colnames(table3),2,12)),4),nucleotide.annotation$query_mass),round(as.numeric(substr(colnames(table3),2,12)),4))]
dim(table3)

# median intensity on each mass
median1 <- apply(table3,MARGIN=2,FUN=median)
head(median1)
length(median1)
median1.n238 <- median1

# save nucleotide profile data
nucleotide.profile.lunge.plattenepithel.n238 <- table3
head(nucleotide.profile.lunge.plattenepithel.n238)


# calculate zero number for each m/z value
len1 <- c()
for (i in 1:ncol(table3)) {
  len <- length(which(table3[,i]==0))
  len1 <- c(len1,len)
}
len1
table4 <- table3[,which(len1 <= 238)]
dim(table4)

table5 <- cbind(table2[,1:2],table4)
dim(table5)
table5[1:5,1:7]

res.cut <- surv_cutpoint(table5, time = "OS", event = "STATOS",
                         minprop = 0.1,
                         variables = colnames(table5)[3:ncol(table5)])
res.cat.1 <- surv_categorize(res.cut,labels = c(0, 1))
dim(res.cat.1)
res.cat.1[1:5,1:7]


res.cat <- cbind(table2[,3:5],res.cat.1)
dim(res.cat)
res.cat[1:5,1:10]
colnames(res.cat)[c(1,3)] <- c("Gender","UICC") #1=men 2=women
res.cat$UICC <- ifelse(res.cat$UICC<5,1,ifelse(res.cat$UICC<7,2,ifelse(res.cat$UICC<10,3,4)))
res.cat$Age <- ifelse(res.cat$Age<70,0,1)
res.cat[1:5,1:10]


# save profile data with cutoff
res.cat.lunge.plattenepithel.n238 <- res.cat
dim(res.cat.lunge.plattenepithel.n238)
head(res.cat.lunge.plattenepithel.n238)


################ univariate cox analysis to calculate p value, HR
result <- data.frame()
for (i in c(1:3,6:165)) {
  res.cox <- coxph(Surv(OS,STATOS)~res.cat[,i],data =res.cat)
  summ <- summary(res.cox)
  hypo <- cox.zph(res.cox)
  hypo1 <- hypo$table[-nrow(hypo$table),]
  result1 <- data.frame(mz.value=colnames(res.cat)[i],
                        HR = signif(summ$coef[2], digits=3),
                        HR.confint.lower = signif(summ$conf.int[,"lower .95"], 3),
                        HR.confint.upper = signif(summ$conf.int[,"upper .95"],3),
                        logtest.p=summ$logtest[3],
                        waldtest.p=summ$waldtest[3],
                        sctest.p=summ$sctest[3],
                        hypo.p=hypo1[3])
  result <- rbind(result,result1)
}
head(result)
dim(result)
univariate.lunge.plattenepithel.n238 <- result


# screen by sctest p value
result2 <- result[result$sctest.p<0.05&result$hypo.p>0.05,]
dim(result2)
head(result2)
result3 <- result2[-c(1:2),]
head(result3)
dim(result3)

# annotation on significant mass
result3$query_mass <- round(as.numeric(substr(result3$mz.value,2,12)),4)
result4 <- merge(result3,nucleotide.annotation)
result4$name.adduct <- paste(result4$adduct,result4$name,sep = '.')
result4$median.intensity <- median1[match(result4$mz.value,colnames(table3))]
dim(result4)
result5 <- result4[order(result4$ppm,decreasing = F),]
result6 <- result5[which(duplicated(result5$name.adduct)==F),]
dim(result6)
result7 <- result6[order(result6$median.intensity,decreasing = T),]
result8 <- result7[which(duplicated(result7$name)==F),]
result8 <- result8[order(result8$query_mass,decreasing = F),]
dim(result8)
head(result8)
univariate4.lunge.plattenepithel.n238 <- result8


############## validate by survival curve

# fit survival curve
fit <- survfit(Surv(OS,STATOS)~Age,data =res.cat)
fit
survdiff(Surv(OS,STATOS)~Age ,data =res.cat) #采用log-rank 检验分析生存率差异

##绘制生存曲线##
ggsurvplot(fit, data = res.cat, conf.int = F, # 增加置信区间
           risk.table = TRUE, pval=T) #不同时间点风险人数表，增加P值



################# forest plot
library(forestplot)

dim(result8)
head(result8)
length(unique(result8$query_mass))
head(result2)

result10 <- data.frame(Feature=c(result2$mz.value[1:2],result8$query_mass),
                       Annotation=c("70","1/2/3/4",result8$name),
                       HR=c(result2$HR[1:2],result8$HR),
                       HR.lower=c(result2$HR.confint.lower[1:2],result8$HR.confint.lower),
                       HR.upper=c(result2$HR.confint.upper[1:2],result8$HR.confint.upper),
                       pvalue=round(c(result2$sctest.p[1:2],result8$sctest.p),3))
dim(result10)
head(result10)

result10$HR1 <- paste(result10$HR,"(",result10$HR.lower,"-",result10$HR.upper,")",sep = '')
dim(result10)
head(result10)

labeltext <- as.matrix(rbind(c("Feature","Annotation","HR","pvalue"),result10[,c(1,2,7,6)]))
labeltext

forestplot(labeltext, 
           mean = c(NA,result10$HR),
           lower = c(NA,result10$HR.lower),
           upper = c(NA,result10$HR.upper),
           zero = 1,lwd.zero = 2,
           clip = c(0,6),
           graph.pos = 3,
           boxsize = 0.5,
           align = "l",
           xlog = FALSE,
           col=fpColors(line = "black", zero = "black",
                        box="black"),
           lty.ci = 7,
           lwd.ci = 1,
           ci.vertices.height = 0.15,
           lwd.xaxis=2,
           xlab="Hazard Ratio",
           txt_gp = fpTxtGp(ticks = gpar(cex = 0.9), xlab = gpar(cex = 0.8), cex = 0.7),
           lineheight = unit(8,"mm")) 





################################################ multivariate cox analysis identify independent prognostic factor
dim(res.cat)
res.cat[1:5,1:10]

head(result2)
head(result8)

res.cat1 <- res.cat[,c("OS","STATOS",result2$mz.value[1:2],unique(result8$mz.value))]
dim(res.cat1)
res.cat1[1:5,1:7]

# save final profiling data after annotation
res.cat1.lung.plattenepithel.n238 <- res.cat1


################## method1 feature selection from the beginning

### Step regression on Cox model
res.cox1 <- coxph(Surv(OS, STATOS) ~ .,data =res.cat1)
res.cox2 <- step(res.cox1)
res.cox3 <- coxph(Surv(OS, STATOS) ~ Age + UICC + X111.0209659 + X192.964059 + 
                    X203.9939506 + X242.0791809 + X267.0735148 + X288.0391331 + 
                    X319.0448899 + X327.0000669 + X350.0162893 + X350.9979773 + 
                    X360.000614 + X368.0379518 + X384.9846068 + X401.027813 + 
                    X410.0257475 + X447.9737532 + X464.9630112 + X500.0189157 + 
                    X520.916907,data =res.cat1)
ggforest(model=res.cox3,data=res.cat1)


###等比例风险假定
hypo <- cox.zph(res.cox3) 
hypo1 <- hypo$table[-nrow(hypo$table),]
dim(hypo1)
hypo2 <- hypo1[hypo1[,3]>0.05,]
dim(hypo2)
hypo2

res.cat2 <- res.cat1[,c("OS","STATOS",rownames(hypo2))]
dim(res.cat2)
res.cat2[1:5,1:7]
res.cox4 <- coxph(Surv(OS, STATOS) ~ .,data = res.cat2)
ggforest(model=res.cox4,data=res.cat2)


summ <- as.data.frame(cbind(summary(res.cox4)$coefficients[,c(2,5)],summary(res.cox4)$conf.int[,c(3,4)]))
summ$mz.value <- rownames(summ)
head(summ)
dim(summ)
multivariate.lunge.plattenepithel.n238.sum <- summ

summ1 <- summ[-c(1:2),]
dim(summ1)
head(summ1)

summ2 <- merge(summ1,result8)
multivariate.lunge.plattenepithel.n238 <- summ2
dim(multivariate.lunge.plattenepithel.n238)
head(multivariate.lunge.plattenepithel.n238)
multivariate.lunge.plattenepithel.n238.sig <- summ2[summ2[,3]<0.05,]
dim(multivariate.lunge.plattenepithel.n238.sig)


#### forest plot
library(forestplot)

result12 <- summ2
dim(result12)
head(result12)

head(summ)
result13 <- summ[1:2,]
head(result13)

result14 <- data.frame(Feature=c(result13$mz.value,result12$query_mass),
                       Annotation=c("70","1/2/3/4",result12$name),
                       HR=c(result13[,1],result12[,2]),
                       HR.lower=c(result13[,3],result12[,4]),
                       HR.upper=c(result13[,4],result12[,5]),
                       pvalue=round(c(result13[,2],result12[,3]),3))
dim(result14)
head(result14)

result14$HR1 <- paste(round(result14$HR,2),"(",round(result14$HR.lower,2),"-",round(result14$HR.upper,2),")",sep = '')
result14$pvalue1 <- ifelse(result14$pvalue<0.001,"***",ifelse(result14$pvalue<0.01,"**",ifelse(result14$pvalue<0.05,"*"," ")))
result14$pvalue2 <- paste(result14$pvalue,result14$pvalue1,sep = '')
dim(result14)
head(result14)

labeltext <- as.matrix(rbind(c("Feature","Annotation","HR","pvalue"),result14[,c(1,2,7,9)]))
labeltext

forestplot(labeltext, 
           mean = c(NA,result14$HR),
           lower = c(NA,result14$HR.lower),
           upper = c(NA,result14$HR.upper),
           zero = 1,lwd.zero = 2,
           clip = c(0,10),
           graph.pos = 3,
           boxsize = 0.5,
           align = "l",
           xlog = F,
           col=fpColors(line = "black", zero = "black",
                        box="black"),
           lty.ci = 7,
           lwd.ci = 1,
           ci.vertices.height = 0.15,
           lwd.xaxis=2,
           xlab="Hazard Ratio",
           txt_gp = fpTxtGp(ticks = gpar(cex = 0.9), xlab = gpar(cex = 0.8), cex = 0.8),
           lineheight = unit(8,"mm")) 



############### use Cox model to predict survival outcome by ROC curve
res.cat3 <- res.cat1[,c("OS","STATOS",rownames(summ))]
dim(res.cat3)
res.cat3[1:5,1:5]

res.cox3 <- coxph(Surv(OS,STATOS)~.,data =res.cat3)
sum2 <- summary(res.cox3)
# C index export
sum2$concordance

# forest plot
ggforest(model = res.cox3,data=res.cat3)

library(pROC)
gfit <- roc(STATOS ~predict(res.cox3), data = res.cat3)

plot(gfit,
     print.auc=TRUE, #输出AUC值
     print.thres=TRUE, #输出cut-off值
     #main = "ROC Curve in Lung NAC tumor (n=77)", #设置图形的标题
     col= "red", #曲线颜色
     print.thres.col="black", #cut-off值字体的颜色
     identity.col="blue", #对角线颜色
     identity.lty=1,identity.lwd=1)


################# use Cox model to predict good/poor prognosis or low/high risk groups by ROC curve
res.cat3$STATOS1 <- ifelse(res.cat3$OS <= median(res.cat3$OS),1,0) #1=short/poor survival or high risk,0=long/good survival or low risk  
res.cat3$STATOS1

gfit1 <- roc(STATOS1 ~predict(res.cox3), data = res.cat3)
cutoff <- gfit1$thresholds[which.max(gfit1$sensitivities+gfit1$specificities)]
cutoff
plot(gfit1,
     print.auc=TRUE, #输出AUC值
     print.thres=TRUE, #输出cut-off值
     #main = "ROC Curve in Lung NAC tumor (n=77)", #设置图形的标题
     col= "red", #曲线颜色
     print.thres.col="black", #cut-off值字体的颜色
     identity.col="blue", #对角线颜色
     identity.lty=1,identity.lwd=1)

### cutoff-based cox model to stratify patients into different prognostic risk groups by survival curve
res.cat3[1:5,1:ncol(res.cat3)]
res.cat3$risk <- ifelse(predict(res.cox3) < cutoff,'low','high')
# fit survival curve
fit <- survfit(Surv(OS,STATOS)~risk,data =res.cat3)
fit
##绘制生存曲线##
ggsurvplot(fit, data = res.cat3, conf.int = F, # 增加置信区间
           risk.table = TRUE, pval=T) #不同时间点风险人数表，增加P值



################# Nomogram绘制

library(car)
library(rms)
library(pROC)
library(timeROC)
library(ggDCA)

res.cat4 <- res.cat1[,c("OS","STATOS",rownames(summ))]
dim(res.cat4)
res.cat4[1:5,]

res.cox4 <- cph(as.formula(paste0('Surv(OS,STATOS)~',paste(colnames(res.cat4)[-c(1:2)],sep = '',collapse = '+'))),
                x=T,y=T,surv=T,data =res.cat4)
print(res.cox4)

ddist <- datadist(res.cat4)
options(datadist='ddist')

surv <- Survival(res.cox4) # 建立生存函数
surv

surv1 <- function(x) surv(12,x) # 1年OS
surv2 <- function(x) surv(36,x) # 3年OS
surv3 <- function(x) surv(60,x) # 5年OS

nom <- nomogram(res.cox4,
                fun = list(surv1,surv2,surv3),
                funlabel = c('1-year survival Probability',
                             '3-year survival Probability',
                             '5-year survival Probability'))
plot(nom, 
     lplabel="Linear Predictor",
     xfrac = 0.2, # 左侧标签距离坐标轴的距离
     #varname.label = TRUE, 
     tcl = -0.2, # 刻度长短和方向 
     lmgp = 0.1, # 坐标轴标签距离坐标轴远近
     points.label ='Points', 
     total.points.label = 'Total Points',
     cap.labels = FALSE,
     cex.var = 0.7, # 左侧标签字体大小
     cex.axis = 0.6, # 坐标轴字体大小
     col.grid = gray(c(0.8, 0.95))) # 竖线颜色


############# ROC curve on the prediction of 1/3/5-year survival rate
res.cat4 <- res.cat1[,c("OS","STATOS",rownames(summ))]
dim(res.cat4)
res.cat4[1:5,1:5]

res.cox4 <- cph(as.formula(paste0('Surv(OS,STATOS)~',paste(colnames(res.cat4)[-c(1:2)],sep = '',collapse = '+'))),
                x=T,y=T,surv=T,data =res.cat4)
print(res.cox4)

pred_f_training <- predict(res.cox4,res.cat4,type="lp")#!!!type="lp",是他没错
data_table <- data.frame(time=res.cat4[,"OS"],status=res.cat4[,"STATOS"],score=pred_f_training)


time_roc_res <- timeROC(
  T = data_table$time,
  delta = data_table$status,
  marker = data_table$score,
  cause = 1,
  weighting="marginal",
  times = c(12, 36, 60),
  ROC = TRUE,
  iid = TRUE
)


time_ROC_df <- data.frame(
  TP_1year = time_roc_res$TP[, 1],
  FP_1year = time_roc_res$FP[, 1],
  TP_3year = time_roc_res$TP[, 2],
  FP_3year = time_roc_res$FP[, 2],
  TP_5year = time_roc_res$TP[, 3],
  FP_5year = time_roc_res$FP[, 3]
)


ggplot(data = time_ROC_df) +
  geom_line(aes(x = FP_1year, y = TP_1year), linewidth = 1, color = "#BC3C29FF") +
  geom_line(aes(x = FP_3year, y = TP_3year), linewidth = 1, color = "#0072B5FF") +
  geom_line(aes(x = FP_5year, y = TP_5year), linewidth = 1, color = "#E18727FF") +
  geom_abline(slope = 1, intercept = 0, color = "grey", linewidth = 1, linetype = 2) +
  theme_bw() +
  annotate("text",
           x = 0.75, y = 0.25, size = 4.5,
           label = paste0("AUC at 1 year = ", sprintf("%.3f", time_roc_res$AUC[[1]])), color = "#BC3C29FF"
  ) +
  annotate("text",
           x = 0.75, y = 0.15, size = 4.5,
           label = paste0("AUC at 3 years = ", sprintf("%.3f", time_roc_res$AUC[[2]])), color = "#0072B5FF"
  ) +
  annotate("text",
           x = 0.75, y = 0.05, size = 4.5,
           label = paste0("AUC at 5 years = ", sprintf("%.3f", time_roc_res$AUC[[3]])), color = "#E18727FF"
  ) +
  labs(x = "Specificity", y = "Sensitivity") +
  theme(
    axis.text = element_text(size = 10, color = "black"),
    axis.title.x = element_text(size = 12, color = "black", margin = margin(c(15, 0, 0, 0))),
    axis.title.y = element_text(size = 12, color = "black", margin = margin(c(0, 15, 0, 0)))
  )

prognosis.n238 <- data.frame(Tumor=rep("Lung squamous cell carcinoma",5),
                             Prognosis=c("Survival outcome","Risk stratification","1-year survival rate","3-year survival rate","5-year survival rate"),
                             AUC=c(gfit$auc,gfit1$auc,time_roc_res$AUC[[1]],time_roc_res$AUC[[2]],time_roc_res$AUC[[3]]))
prognosis.n238





############################## all primary EAC (n=102) #####################################
library(survival)
library(survminer)

setwd("C:/MyRdata8/All primary EAC (n=102)")
table <- read.csv("All primary EAC (n=102).csv")
dim(table)
table[1:5,1:15]

table1 <- table[,c(1,2,3,4,5,13:6959)]
dim(table1)
table1[1:5,1:7]

rownames(table1) <- table1[,1]
table2 <- table1[,-1]
dim(table2)
table2[1:5,1:8]

table3 <- table2[,-c(1:6)]
dim(table3)
table3[1:5,1:7]
table3.n102 <- table3

# nucleotide profile data
table3 <- table3[,match(intersect(round(as.numeric(substr(colnames(table3),2,12)),4),nucleotide.annotation$query_mass),round(as.numeric(substr(colnames(table3),2,12)),4))]
dim(table3)

# median intensity on each mass
median1 <- apply(table3,MARGIN=2,FUN=median)
head(median1)
length(median1)
median1.n102 <- median1

# save nucleotide profile data
nucleotide.profile.all.primary.EAC.n102 <- table3
head(nucleotide.profile.all.primary.EAC.n102)


# calculate zero number for each m/z value
len1 <- c()
for (i in 1:ncol(table3)) {
  len <- length(which(table3[,i]==0))
  len1 <- c(len1,len)
}
len1
table4 <- table3[,which(len1 <= 101)]
dim(table4)

table5 <- cbind(table2[,5:6],table4)
dim(table5)
table5[1:5,1:7]

res.cut <- surv_cutpoint(table5, time = "OS", event = "STATOS",
                         minprop = 0.1,
                         variables = colnames(table5)[3:ncol(table5)])
res.cat.1 <- surv_categorize(res.cut,labels = c(0, 1))
dim(res.cat.1)
res.cat.1[1:5,1:7]

res.cat <- cbind(table2[,1:4],res.cat.1)
dim(res.cat)
res.cat[1:5,1:10]

colnames(res.cat)[3] <- "UICC"
res.cat$UICC <- ifelse(res.cat$yt==0&res.cat$yn==0&res.cat$ym==0,0,
                       ifelse(res.cat$yt==1&res.cat$yn==0&res.cat$ym==0,1,
                              ifelse(res.cat$ym==1,4,
                                     ifelse(res.cat$ym==0&((res.cat$yt==0&res.cat$yn==1)|(res.cat$yt==1&res.cat$yn==1)|(res.cat$yt==2&res.cat$yn==0)|(res.cat$yt==2&res.cat$yn==1)|(res.cat$yt==3&res.cat$yn==0)),2,3))))
res.cat[1:5,1:12]
dim(res.cat)
res.cat <- res.cat[,-c(1,2,4)]


# save profile data with cutoff
res.cat.all.primary.EAC.n102 <- res.cat
dim(res.cat.all.primary.EAC.n102)
head(res.cat.all.primary.EAC.n102)


################ univariate cox analysis to calculate p value, HR
result <- data.frame()
for (i in c(1,4:ncol(res.cat))) {
  res.cox <- coxph(Surv(OS,STATOS)~res.cat[,i],data =res.cat)
  summ <- summary(res.cox)
  hypo <- cox.zph(res.cox)
  hypo1 <- hypo$table[-nrow(hypo$table),]
  result1 <- data.frame(mz.value=colnames(res.cat)[i],
                        HR = signif(summ$coef[2], digits=3),
                        HR.confint.lower = signif(summ$conf.int[,"lower .95"], 3),
                        HR.confint.upper = signif(summ$conf.int[,"upper .95"],3),
                        logtest.p=summ$logtest[3],
                        waldtest.p=summ$waldtest[3],
                        sctest.p=summ$sctest[3],
                        hypo.p=hypo1[3])
  result <- rbind(result,result1)
}
head(result)
dim(result)
univariate.all.primary.EAC.n102 <- result


# screen by sctest p value
result2 <- result[result$sctest.p<0.05&result$hypo.p>0.05,]
dim(result2)
head(result2)
result3 <- result2[-c(1),]
head(result3)
dim(result3)

# annotation on significant mass
result3$query_mass <- round(as.numeric(substr(result3$mz.value,2,12)),4)
result4 <- merge(result3,nucleotide.annotation)
result4$name.adduct <- paste(result4$adduct,result4$name,sep = '.')
result4$median.intensity <- median1[match(result4$mz.value,colnames(table3))]
dim(result4)
result5 <- result4[order(result4$ppm,decreasing = F),]
result6 <- result5[which(duplicated(result5$name.adduct)==F),]
dim(result6)
result7 <- result6[order(result6$median.intensity,decreasing = T),]
result8 <- result7[which(duplicated(result7$name)==F),]
result8 <- result8[order(result8$query_mass,decreasing = F),]
dim(result8)
head(result8)
univariate4.all.primary.EAC.n102 <- result8



############## validate by survival curve
# fit survival curve
fit <- survfit(Surv(OS,STATOS)~UICC ,data =res.cat)
fit
survdiff(Surv(OS,STATOS)~UICC ,data =res.cat) #采用log-rank 检验分析生存率差异

##绘制生存曲线##
ggsurvplot(fit, data = res.cat, conf.int = F, # 增加置信区间
           risk.table = TRUE, pval=T) #不同时间点风险人数表，增加P值



################# forest plot
library(forestplot)

dim(result8)
head(result8)
length(unique(result8$query_mass))
head(result2)

result10 <- data.frame(Feature=c(result2$mz.value[1],result8$query_mass),
                       Annotation=c("1/2/3/4",result8$name),
                       HR=c(result2$HR[1],result8$HR),
                       HR.lower=c(result2$HR.confint.lower[1],result8$HR.confint.lower),
                       HR.upper=c(result2$HR.confint.upper[1],result8$HR.confint.upper),
                       pvalue=round(c(result2$sctest.p[1],result8$sctest.p),3))
dim(result10)
head(result10)

result10$HR1 <- paste(result10$HR,"(",result10$HR.lower,"-",result10$HR.upper,")",sep = '')
dim(result10)
head(result10)

labeltext <- as.matrix(rbind(c("Feature","Annotation","HR","pvalue"),result10[,c(1,2,7,6)]))
labeltext

forestplot(labeltext, 
           mean = c(NA,result10$HR),
           lower = c(NA,result10$HR.lower),
           upper = c(NA,result10$HR.upper),
           zero = 1,lwd.zero = 2,
           clip = c(0,6),
           graph.pos = 3,
           boxsize = 0.5,
           align = "l",
           xlog = FALSE,
           col=fpColors(line = "black", zero = "black",
                        box="black"),
           lty.ci = 7,
           lwd.ci = 1,
           ci.vertices.height = 0.15,
           lwd.xaxis=2,
           xlab="Hazard Ratio",
           txt_gp = fpTxtGp(ticks = gpar(cex = 0.9), xlab = gpar(cex = 0.8), cex = 0.7),
           lineheight = unit(8,"mm")) 


################################################ multivariate cox analysis identify independent prognostic factor
dim(res.cat)
res.cat[1:5,1:10]

head(result2)
head(result8)

res.cat1 <- res.cat[,c("OS","STATOS",result2$mz.value[1],unique(result8$mz.value))]
dim(res.cat1)
res.cat1[1:5,1:7]

# save final profiling data after annotation
res.cat1.all.primary.EAC.n102 <- res.cat1


################## method1 feature selection from the beginning
############### lasso regression
res.cat2 <-res.cat1[res.cat1$OS!=0,]
dim(res.cat2)
res.cat2[1:5,1:10]

#调用glmnet包中的glmnet函数，注意family那里一定要制定是“cox”，如果是做logistic需要换成"binomial"。
library(glmnet)
x <- as.matrix(res.cat2[,-c(1,2)])
y <- data.matrix(Surv(res.cat2$OS,res.cat2$STATOS))

#主要在做交叉验证,lasso
fitcv <- cv.glmnet(x,y,family="cox", alpha=1,nfolds=10)
plot(fitcv)
coef <- coef(fitcv, s="lambda.min")
coef1 <- coef[which(coef!=0),]
names(coef1)
length(names(coef1)) #15

#res.cat1 <- res.cat[,c("OS","STATOS",names(coef1))]
dim(res.cat1)
res.cat1[1:5,1:10]


### Step regression on Cox model
res.cox1 <- coxph(Surv(OS, STATOS) ~ .,data =res.cat1)
res.cox2 <- step(res.cox1)
res.cox3 <- coxph(Surv(OS, STATOS) ~ UICC + X150.0424563 + X152.9959439 + X192.9646429 + 
                    X203.9944491 + X228.0272234 + X342.0254651 + X392.0179625 + 
                    X402.9943229 + X421.9957971,data =res.cat1)
ggforest(model=res.cox3,data=res.cat1)


###等比例风险假定
hypo <- cox.zph(res.cox3) 
#hypo <- cox.zph(res.cox4) 
hypo1 <- hypo$table[-nrow(hypo$table),]
dim(hypo1)
hypo2 <- hypo1[hypo1[,3]>0.05,]
dim(hypo2)
hypo2

res.cat2 <- res.cat1[,c("OS","STATOS",rownames(hypo2))]
dim(res.cat2)
res.cat2[1:5,1:7]
res.cox4 <- coxph(Surv(OS, STATOS) ~ .,data = res.cat2)
ggforest(model=res.cox4,data=res.cat2)


summ <- as.data.frame(cbind(summary(res.cox4)$coefficients[,c(2,5)],summary(res.cox4)$conf.int[,c(3,4)]))
summ$mz.value <- rownames(summ)
head(summ)
dim(summ)
multivariate.all.primary.EAC.n102.sum <- summ

summ1 <- summ[-c(1),]
dim(summ1)
head(summ1)

summ2 <- merge(summ1,result8)
multivariate.all.primary.EAC.n102 <- summ2
dim(multivariate.all.primary.EAC.n102)
head(multivariate.all.primary.EAC.n102)
multivariate.all.primary.EAC.n102.sig <- summ2[summ2[,3]<0.05,]
dim(multivariate.all.primary.EAC.n102.sig)


#### forest plot
library(forestplot)

result12 <- summ2
dim(result12)
head(result12)

head(summ)
result13 <- summ[1,]
head(result13)

result14 <- data.frame(Feature=c(result13$mz.value,result12$query_mass),
                       Annotation=c("1/2/3/4",result12$name),
                       HR=c(result13[,1],result12[,2]),
                       HR.lower=c(result13[,3],result12[,4]),
                       HR.upper=c(result13[,4],result12[,5]),
                       pvalue=round(c(result13[,2],result12[,3]),3))
dim(result14)
head(result14)

result14$HR1 <- paste(round(result14$HR,2),"(",round(result14$HR.lower,2),"-",round(result14$HR.upper,2),")",sep = '')
result14$pvalue1 <- ifelse(result14$pvalue<0.001,"***",ifelse(result14$pvalue<0.01,"**",ifelse(result14$pvalue<0.05,"*"," ")))
result14$pvalue2 <- paste(result14$pvalue,result14$pvalue1,sep = '')
dim(result14)
head(result14)

labeltext <- as.matrix(rbind(c("Feature","Annotation","HR","pvalue"),result14[,c(1,2,7,9)]))
labeltext

forestplot(labeltext, 
           mean = c(NA,result14$HR),
           lower = c(NA,result14$HR.lower),
           upper = c(NA,result14$HR.upper),
           zero = 1,lwd.zero = 2,
           clip = c(0.01,50),
           graph.pos = 3,
           boxsize = 0.5,
           align = "l",
           xlog = T,
           col=fpColors(line = "black", zero = "black",
                        box="black"),
           lty.ci = 7,
           lwd.ci = 1,
           ci.vertices.height = 0.15,
           lwd.xaxis=2,
           xlab="Hazard Ratio",
           txt_gp = fpTxtGp(ticks = gpar(cex = 0.9), xlab = gpar(cex = 0.8), cex = 0.8),
           lineheight = unit(8,"mm")) 




############### use Cox model to predict survival outcome by ROC curve
res.cat3 <- res.cat1[,c("OS","STATOS",rownames(summ))]
dim(res.cat3)
res.cat3[1:5,1:5]

res.cox3 <- coxph(Surv(OS,STATOS)~.,data =res.cat3)
sum2 <- summary(res.cox3)
# C index export
sum2$concordance

# forest plot
ggforest(model = res.cox3,data=res.cat3)

library(pROC)
gfit <- roc(STATOS ~predict(res.cox3), data = res.cat3)

plot(gfit,
     print.auc=TRUE, #输出AUC值
     print.thres=TRUE, #输出cut-off值
     #main = "ROC Curve in Lung NAC tumor (n=77)", #设置图形的标题
     col= "red", #曲线颜色
     print.thres.col="black", #cut-off值字体的颜色
     identity.col="blue", #对角线颜色
     identity.lty=1,identity.lwd=1)


################# use Cox model to predict good/poor prognosis or low/high risk groups by ROC curve
res.cat3$STATOS1 <- ifelse(res.cat3$OS < median(res.cat3$OS),1,0) #1=short/poor survival or high risk,0=long/good survival or low risk  
res.cat3$STATOS1

gfit1 <- roc(STATOS1 ~predict(res.cox3), data = res.cat3)
cutoff <- gfit1$thresholds[which.max(gfit1$sensitivities+gfit1$specificities)]
cutoff
plot(gfit1,
     print.auc=TRUE, #输出AUC值
     print.thres=T, #输出cut-off值
     #main = "ROC Curve in Lung NAC tumor (n=77)", #设置图形的标题
     col= "red", #曲线颜色
     print.thres.col="black", #cut-off值字体的颜色
     identity.col="blue", #对角线颜色
     identity.lty=1,identity.lwd=1)

### cutoff-based cox model to stratify patients into different prognostic risk groups by survival curve
res.cat3[1:5,1:ncol(res.cat3)]
res.cat3$risk <- ifelse(predict(res.cox3) < cutoff,'low','high')
# fit survival curve
fit <- survfit(Surv(OS,STATOS)~risk,data =res.cat3)
fit
##绘制生存曲线##
ggsurvplot(fit, data = res.cat3, conf.int = F, # 增加置信区间
           risk.table = TRUE, pval=T) #不同时间点风险人数表，增加P值



################# Nomogram绘制

library(car)
library(rms)
library(pROC)
library(timeROC)
library(ggDCA)

res.cat4 <- res.cat1[,c("OS","STATOS",rownames(summ))]
dim(res.cat4)
res.cat4[1:5,]

res.cox4 <- cph(as.formula(paste0('Surv(OS,STATOS)~',paste(colnames(res.cat4)[-c(1:2)],sep = '',collapse = '+'))),
                x=T,y=T,surv=T,data =res.cat4)
print(res.cox4)

ddist <- datadist(res.cat4)
options(datadist='ddist')

surv <- Survival(res.cox4) # 建立生存函数
surv

surv1 <- function(x) surv(12,x) # 1年OS
surv2 <- function(x) surv(36,x) # 3年OS
surv3 <- function(x) surv(60,x) # 5年OS

nom <- nomogram(res.cox4,
                fun = list(surv1,surv2,surv3),
                funlabel = c('1-year survival Probability',
                             '3-year survival Probability',
                             '5-year survival Probability'))
plot(nom, 
     lplabel="Linear Predictor",
     xfrac = 0.2, # 左侧标签距离坐标轴的距离
     #varname.label = TRUE, 
     tcl = -0.2, # 刻度长短和方向 
     lmgp = 0.1, # 坐标轴标签距离坐标轴远近
     points.label ='Points', 
     total.points.label = 'Total Points',
     cap.labels = FALSE,
     cex.var = 0.7, # 左侧标签字体大小
     cex.axis = 0.6, # 坐标轴字体大小
     col.grid = gray(c(0.8, 0.95))) # 竖线颜色


############# ROC curve on the prediction of 1/3/5-year survival rate
res.cat4 <- res.cat1[,c("OS","STATOS",rownames(summ))]
dim(res.cat4)
res.cat4[1:5,1:5]

res.cox4 <- cph(as.formula(paste0('Surv(OS,STATOS)~',paste(colnames(res.cat4)[-c(1:2)],sep = '',collapse = '+'))),
                x=T,y=T,surv=T,data =res.cat4)
print(res.cox4)

pred_f_training <- predict(res.cox4,res.cat4,type="lp")#!!!type="lp",是他没错
data_table <- data.frame(time=res.cat4[,"OS"],status=res.cat4[,"STATOS"],score=pred_f_training)


time_roc_res <- timeROC(
  T = data_table$time,
  delta = data_table$status,
  marker = data_table$score,
  cause = 1,
  weighting="marginal",
  times = c(12, 36, 60),
  ROC = TRUE,
  iid = TRUE
)


time_ROC_df <- data.frame(
  TP_1year = time_roc_res$TP[, 1],
  FP_1year = time_roc_res$FP[, 1],
  TP_3year = time_roc_res$TP[, 2],
  FP_3year = time_roc_res$FP[, 2],
  TP_5year = time_roc_res$TP[, 3],
  FP_5year = time_roc_res$FP[, 3]
)


ggplot(data = time_ROC_df) +
  geom_line(aes(x = FP_1year, y = TP_1year), linewidth = 1, color = "#BC3C29FF") +
  geom_line(aes(x = FP_3year, y = TP_3year), linewidth = 1, color = "#0072B5FF") +
  geom_line(aes(x = FP_5year, y = TP_5year), linewidth = 1, color = "#E18727FF") +
  geom_abline(slope = 1, intercept = 0, color = "grey", linewidth = 1, linetype = 2) +
  theme_bw() +
  annotate("text",
           x = 0.75, y = 0.25, size = 4.5,
           label = paste0("AUC at 1 year = ", sprintf("%.3f", time_roc_res$AUC[[1]])), color = "#BC3C29FF"
  ) +
  annotate("text",
           x = 0.75, y = 0.15, size = 4.5,
           label = paste0("AUC at 3 years = ", sprintf("%.3f", time_roc_res$AUC[[2]])), color = "#0072B5FF"
  ) +
  annotate("text",
           x = 0.75, y = 0.05, size = 4.5,
           label = paste0("AUC at 5 years = ", sprintf("%.3f", time_roc_res$AUC[[3]])), color = "#E18727FF"
  ) +
  labs(x = "Specificity", y = "Sensitivity") +
  theme(
    axis.text = element_text(size = 10, color = "black"),
    axis.title.x = element_text(size = 12, color = "black", margin = margin(c(15, 0, 0, 0))),
    axis.title.y = element_text(size = 12, color = "black", margin = margin(c(0, 15, 0, 0)))
  )


prognosis.n102 <- data.frame(Tumor=rep("Primary EAC",5),
                             Prognosis=c("Survival outcome","Risk stratification","1-year survival rate","3-year survival rate","5-year survival rate"),
                             AUC=c(gfit$auc,gfit1$auc,time_roc_res$AUC[[1]],time_roc_res$AUC[[2]],time_roc_res$AUC[[3]]))
prognosis.n102


############################## Neoadjuvant treated EACs (n=144) #####################################
library(survival)
library(survminer)

setwd("C:/MyRdata8/Neoadjuvant treated EACs (n=144)")

table <- read.csv("Neoadjuvant treated EACs (n=144).csv")
dim(table)
table[1:5,1:7]

table$Block.number[match(c("B07_26391", "B11_15981", "B11_8507"),table$Block.number)] <- c("B07_26391.1", "B11_15981.1", "B11_8507.1")

table1 <- table
dim(table1)
table1[1:5,1:7]

rownames(table1) <- table1[,1]
table2 <- table1[,-1]
dim(table2)
table2[1:5,1:12]
surv.n144 <- table2

table3 <- table2[,-c(1:11)]
dim(table3)
table3[1:5,1:7]
table3.n144 <- table3

# nucleotide profile data
table3 <- table3[,match(intersect(round(as.numeric(substr(colnames(table3),2,12)),4),nucleotide.annotation$query_mass),round(as.numeric(substr(colnames(table3),2,12)),4))]
dim(table3)

# median intensity on each mass
median1 <- apply(table3,MARGIN=2,FUN=median)
head(median1)
length(median1)
median1.n144 <- median1

# save nucleotide profile data
nucleotide.profile.Neoadjuvant.treated.EACs.n144 <- table3
head(nucleotide.profile.Neoadjuvant.treated.EACs.n144)


# calculate zero number for each m/z value
len1 <- c()
for (i in 1:ncol(table3)) {
  len <- length(which(table3[,i]==0))
  len1 <- c(len1,len)
}
len1
table4 <- table3[,which(len1 <= 134)]
dim(table4)

table5 <- cbind(table2[,10:11],table4)
dim(table5)
table5[1:5,1:7]

res.cut <- surv_cutpoint(table5, time = "OS", event = "STATOS",
                         minprop = 0.1,
                         variables = colnames(table5)[3:ncol(table5)])
res.cat.1 <- surv_categorize(res.cut,labels = c(0, 1))
dim(res.cat.1)
res.cat.1[1:5,1:7]

table2[1:5,1:12]
res.cat <- cbind(table2[,1:9],res.cat.1)
dim(res.cat)
res.cat[1:5,1:12]

colnames(res.cat)[9] <- "UICC"
res.cat$UICC <- ifelse(res.cat$yT.im.TNM==0&res.cat$N0.versus.N1==0&res.cat$yM.im.TNM==0,0,
                       ifelse(res.cat$yT.im.TNM==1&res.cat$yN.im.TNM==0&res.cat$yM.im.TNM==0,1,
                              ifelse(res.cat$yM.im.TNM==1,4,
                                     ifelse(res.cat$yM.im.TNM==0&((res.cat$yT.im.TNM==0&res.cat$yN.im.TNM==1)|(res.cat$yT.im.TNM==1&res.cat$yN.im.TNM==1)|(res.cat$yT.im.TNM==2&res.cat$yN.im.TNM==0)|(res.cat$yT.im.TNM==2&res.cat$yN.im.TNM==1)|(res.cat$yT.im.TNM==3&res.cat$yN.im.TNM==0)),2,3))))
res.cat[1:5,1:12]
dim(res.cat)

res.cat <- res.cat[,c(1,8:ncol(res.cat))]
dim(res.cat)
res.cat[1:5,1:10]
clin.n144 <- res.cat

# save profile data with cutoff
res.cat.Neoadjuvant.treated.EACs.n144 <- res.cat
dim(res.cat.Neoadjuvant.treated.EACs.n144)
head(res.cat.Neoadjuvant.treated.EACs.n144)


################ univariate cox analysis to calculate p value, HR
result <- data.frame()
for (i in c(1:3,6:ncol(res.cat))) {
  res.cox <- coxph(Surv(OS,STATOS)~res.cat[,i],data =res.cat)
  summ <- summary(res.cox)
  hypo <- cox.zph(res.cox)
  hypo1 <- hypo$table[-nrow(hypo$table),]
  result1 <- data.frame(mz.value=colnames(res.cat)[i],
                        HR = signif(summ$coef[2], digits=3),
                        HR.confint.lower = signif(summ$conf.int[,"lower .95"], 3),
                        HR.confint.upper = signif(summ$conf.int[,"upper .95"],3),
                        logtest.p=summ$logtest[3],
                        waldtest.p=summ$waldtest[3],
                        sctest.p=summ$sctest[3],
                        hypo.p=hypo1[3])
  result <- rbind(result,result1)
}
head(result)
dim(result)
univariate.Neoadjuvant.treated.EACs.n144 <- result


# screen by sctest p value
result2 <- result[result$sctest.p<0.05&result$hypo.p>0.05,]
dim(result2)
head(result2)
result3 <- result2[-c(1:2),]
head(result3)
dim(result3)

# annotation on significant mass
result3$query_mass <- round(as.numeric(substr(result3$mz.value,2,12)),4)
result4 <- merge(result3,nucleotide.annotation)
result4$name.adduct <- paste(result4$adduct,result4$name,sep = '.')
result4$median.intensity <- median1[match(result4$mz.value,colnames(table3))]
dim(result4)
result5 <- result4[order(result4$ppm,decreasing = F),]
result6 <- result5[which(duplicated(result5$name.adduct)==F),]
dim(result6)
result7 <- result6[order(result6$median.intensity,decreasing = T),]
result8 <- result7[which(duplicated(result7$name)==F),]
result8 <- result8[order(result8$query_mass,decreasing = F),]
dim(result8)
head(result8)
univariate4.Neoadjuvant.treated.EACs.n144 <- result8



############## validate by survival curve
res.cat[1:5,1:10]
# fit survival curve
fit <- survfit(Surv(OS,STATOS)~Becker.Resp.vs.Non.resp,data =res.cat)
fit
survdiff(Surv(OS,STATOS)~Becker.Resp.vs.Non.resp,data =res.cat) #采用log-rank 检验分析生存率差异

##绘制生存曲线##
ggsurvplot(fit, data = res.cat, conf.int = F, # 增加置信区间
           risk.table = TRUE, pval=T) #不同时间点风险人数表，增加P值


################# forest plot
library(forestplot)

dim(result8)
head(result8)
length(unique(result8$query_mass))
head(result2)

result10 <- data.frame(Feature=c(result2$mz.value[1:2],result8$query_mass),
                       Annotation=c("0/1","0/1/2/3/4",result8$name),
                       HR=c(result2$HR[1:2],result8$HR),
                       HR.lower=c(result2$HR.confint.lower[1:2],result8$HR.confint.lower),
                       HR.upper=c(result2$HR.confint.upper[1:2],result8$HR.confint.upper),
                       pvalue=round(c(result2$sctest.p[1:2],result8$sctest.p),3))
dim(result10)
head(result10)

result10$HR1 <- paste(result10$HR,"(",result10$HR.lower,"-",result10$HR.upper,")",sep = '')
dim(result10)
head(result10)

labeltext <- as.matrix(rbind(c("Feature","Annotation","HR","pvalue"),result10[,c(1,2,7,6)]))
labeltext

forestplot(labeltext, 
           mean = c(NA,result10$HR),
           lower = c(NA,result10$HR.lower),
           upper = c(NA,result10$HR.upper),
           zero = 1,lwd.zero = 2,
           clip = c(0,6),
           graph.pos = 3,
           boxsize = 0.5,
           align = "l",
           xlog = FALSE,
           col=fpColors(line = "black", zero = "black",
                        box="black"),
           lty.ci = 7,
           lwd.ci = 1,
           ci.vertices.height = 0.15,
           lwd.xaxis=2,
           xlab="Hazard Ratio",
           txt_gp = fpTxtGp(ticks = gpar(cex = 0.9), xlab = gpar(cex = 0.8), cex = 0.7),
           lineheight = unit(8,"mm")) 



################################################ multivariate cox analysis identify independent prognostic factor
dim(res.cat)
res.cat[1:5,1:10]

head(result2)
head(result8)

res.cat1 <- res.cat[,c("OS","STATOS",result2$mz.value[1:2],unique(result8$mz.value))]
dim(res.cat1)
res.cat1[1:5,1:7]

# save final profiling data after annotation
res.cat1.Neoadjuvant.treated.EACs.n144 <- res.cat1


################## method1 feature selection from the beginning

### Step regression on Cox model
res.cox1 <- coxph(Surv(OS, STATOS) ~ .,data =res.cat1)
res.cox2 <- step(res.cox1)
res.cox3 <- coxph(Surv(OS, STATOS) ~ UICC + X282.0853485 + X323.0275529 + X330.0602309 + 
                    X384.0313623 + X442.0187089 + X505.9845464,data =res.cat1)
ggforest(model=res.cox3,data=res.cat1)


###等比例风险假定
hypo <- cox.zph(res.cox3) 
hypo1 <- hypo$table[-nrow(hypo$table),]
dim(hypo1)
hypo2 <- hypo1[hypo1[,3]>0.05,]
dim(hypo2)
hypo2

res.cat2 <- res.cat1[,c("OS","STATOS",rownames(hypo2))]
dim(res.cat2)
res.cat2[1:5,1:7]
res.cox4 <- coxph(Surv(OS, STATOS) ~ .,data = res.cat2)
ggforest(model=res.cox4,data=res.cat2)


summ <- as.data.frame(cbind(summary(res.cox4)$coefficients[,c(2,5)],summary(res.cox4)$conf.int[,c(3,4)]))
summ$mz.value <- rownames(summ)
head(summ)
dim(summ)
multivariate.Neoadjuvant.treated.EACs.n144.sum <- summ

summ1 <- summ[-c(1),]
dim(summ1)
head(summ1)

summ2 <- merge(summ1,result8)
multivariate.Neoadjuvant.treated.EACs.n144 <- summ2
dim(multivariate.Neoadjuvant.treated.EACs.n144)
head(multivariate.Neoadjuvant.treated.EACs.n144)
multivariate.Neoadjuvant.treated.EACs.n144.sig <- summ2[summ2[,3]<0.05,]
dim(multivariate.Neoadjuvant.treated.EACs.n144.sig)


#### forest plot
library(forestplot)

result12 <- summ2
dim(result12)
head(result12)

head(summ)
result13 <- summ[1,]
head(result13)

result14 <- data.frame(Feature=c(result13$mz.value,result12$query_mass),
                       Annotation=c("0/1/2/3/4",result12$name),
                       HR=c(result13[,1],result12[,2]),
                       HR.lower=c(result13[,3],result12[,4]),
                       HR.upper=c(result13[,4],result12[,5]),
                       pvalue=round(c(result13[,2],result12[,3]),3))
dim(result14)
head(result14)

result14$HR1 <- paste(round(result14$HR,2),"(",round(result14$HR.lower,2),"-",round(result14$HR.upper,2),")",sep = '')
result14$pvalue1 <- ifelse(result14$pvalue<0.001,"***",ifelse(result14$pvalue<0.01,"**",ifelse(result14$pvalue<0.05,"*"," ")))
result14$pvalue2 <- paste(result14$pvalue,result14$pvalue1,sep = '')
dim(result14)
head(result14)

labeltext <- as.matrix(rbind(c("Feature","Annotation","HR","pvalue"),result14[,c(1,2,7,9)]))
labeltext

forestplot(labeltext, 
           mean = c(NA,result14$HR),
           lower = c(NA,result14$HR.lower),
           upper = c(NA,result14$HR.upper),
           zero = 1,lwd.zero = 2,
           clip = c(0.01,40),
           graph.pos = 3,
           boxsize = 0.5,
           align = "l",
           xlog = T,
           col=fpColors(line = "black", zero = "black",
                        box="black"),
           lty.ci = 7,
           lwd.ci = 1,
           ci.vertices.height = 0.15,
           lwd.xaxis=2,
           xlab="Hazard Ratio",
           txt_gp = fpTxtGp(ticks = gpar(cex = 0.9), xlab = gpar(cex = 0.8), cex = 0.8),
           lineheight = unit(8,"mm")) 




############### use Cox model to predict survival outcome by ROC curve
res.cat3 <- res.cat1[,c("OS","STATOS",rownames(summ))]
dim(res.cat3)
res.cat3[1:5,1:5]

res.cox3 <- coxph(Surv(OS,STATOS)~.,data =res.cat3)
sum2 <- summary(res.cox3)
# C index export
sum2$concordance

# forest plot
ggforest(model = res.cox3,data=res.cat3)

library(pROC)
gfit <- roc(STATOS ~predict(res.cox3), data = res.cat3)

plot(gfit,
     print.auc=TRUE, #输出AUC值
     print.thres=TRUE, #输出cut-off值
     #main = "ROC Curve in Lung NAC tumor (n=77)", #设置图形的标题
     col= "red", #曲线颜色
     print.thres.col="black", #cut-off值字体的颜色
     identity.col="blue", #对角线颜色
     identity.lty=1,identity.lwd=1)


################# use Cox model to predict good/poor prognosis or low/high risk groups by ROC curve
res.cat3$STATOS1 <- ifelse(res.cat3$OS < median(res.cat3$OS),1,0) #1=short/poor survival or high risk,0=long/good survival or low risk  
res.cat3$STATOS1

gfit1 <- roc(STATOS1 ~predict(res.cox3), data = res.cat3)
cutoff <- gfit1$thresholds[which.max(gfit1$sensitivities+gfit1$specificities)]
cutoff
plot(gfit1,
     print.auc=TRUE, #输出AUC值
     print.thres=TRUE, #输出cut-off值
     #main = "ROC Curve in Lung NAC tumor (n=77)", #设置图形的标题
     col= "red", #曲线颜色
     print.thres.col="black", #cut-off值字体的颜色
     identity.col="blue", #对角线颜色
     identity.lty=1,identity.lwd=1)

### cutoff-based cox model to stratify patients into different prognostic risk groups by survival curve
res.cat3[1:5,1:ncol(res.cat3)]
res.cat3$risk <- ifelse(predict(res.cox3) < cutoff,'low','high')
# fit survival curve
fit <- survfit(Surv(OS,STATOS)~risk,data =res.cat3)
fit
##绘制生存曲线##
ggsurvplot(fit, data = res.cat3, conf.int = F, # 增加置信区间
           risk.table = TRUE, pval=T) #不同时间点风险人数表，增加P值



################# Nomogram绘制

library(car)
library(rms)
library(pROC)
library(timeROC)
library(ggDCA)

res.cat4 <- res.cat1[,c("OS","STATOS",rownames(summ))]
dim(res.cat4)
res.cat4[1:5,]

res.cox4 <- cph(as.formula(paste0('Surv(OS,STATOS)~',paste(colnames(res.cat4)[-c(1:2)],sep = '',collapse = '+'))),
                x=T,y=T,surv=T,data =res.cat4)
print(res.cox4)

ddist <- datadist(res.cat4)
options(datadist='ddist')

surv <- Survival(res.cox4) # 建立生存函数
surv

surv1 <- function(x) surv(12,x) # 1年OS
surv2 <- function(x) surv(36,x) # 3年OS
surv3 <- function(x) surv(60,x) # 5年OS

nom <- nomogram(res.cox4,
                fun = list(surv1,surv2,surv3),
                funlabel = c('1-year survival Probability',
                             '3-year survival Probability',
                             '5-year survival Probability'))
plot(nom, 
     lplabel="Linear Predictor",
     xfrac = 0.2, # 左侧标签距离坐标轴的距离
     #varname.label = TRUE, 
     tcl = -0.2, # 刻度长短和方向 
     lmgp = 0.1, # 坐标轴标签距离坐标轴远近
     points.label ='Points', 
     total.points.label = 'Total Points',
     cap.labels = FALSE,
     cex.var = 0.7, # 左侧标签字体大小
     cex.axis = 0.6, # 坐标轴字体大小
     col.grid = gray(c(0.8, 0.95))) # 竖线颜色


############# ROC curve on the prediction of 1/3/5-year survival rate
library(timeROC)
res.cat4 <- res.cat1[,c("OS","STATOS",rownames(summ))]
dim(res.cat4)
res.cat4[1:5,1:5]

res.cox4 <- cph(as.formula(paste0('Surv(OS,STATOS)~',paste(colnames(res.cat4)[-c(1:2)],sep = '',collapse = '+'))),
                x=T,y=T,surv=T,data =res.cat4)
print(res.cox4)

pred_f_training <- predict(res.cox4,res.cat4,type="lp")#!!!type="lp",是他没错
data_table <- data.frame(time=res.cat4[,"OS"],status=res.cat4[,"STATOS"],score=pred_f_training)


time_roc_res <- timeROC(
  T = data_table$time,
  delta = data_table$status,
  marker = data_table$score,
  cause = 1,
  weighting="marginal",
  times = c(12, 36, 60),
  ROC = TRUE,
  iid = TRUE
)


time_ROC_df <- data.frame(
  TP_1year = time_roc_res$TP[, 1],
  FP_1year = time_roc_res$FP[, 1],
  TP_3year = time_roc_res$TP[, 2],
  FP_3year = time_roc_res$FP[, 2],
  TP_5year = time_roc_res$TP[, 3],
  FP_5year = time_roc_res$FP[, 3]
)


ggplot(data = time_ROC_df) +
  geom_line(aes(x = FP_1year, y = TP_1year), linewidth = 1, color = "#BC3C29FF") +
  geom_line(aes(x = FP_3year, y = TP_3year), linewidth = 1, color = "#0072B5FF") +
  geom_line(aes(x = FP_5year, y = TP_5year), linewidth = 1, color = "#E18727FF") +
  geom_abline(slope = 1, intercept = 0, color = "grey", linewidth = 1, linetype = 2) +
  theme_bw() +
  annotate("text",
           x = 0.75, y = 0.25, size = 4.5,
           label = paste0("AUC at 1 year = ", sprintf("%.3f", time_roc_res$AUC[[1]])), color = "#BC3C29FF"
  ) +
  annotate("text",
           x = 0.75, y = 0.15, size = 4.5,
           label = paste0("AUC at 3 years = ", sprintf("%.3f", time_roc_res$AUC[[2]])), color = "#0072B5FF"
  ) +
  annotate("text",
           x = 0.75, y = 0.05, size = 4.5,
           label = paste0("AUC at 5 years = ", sprintf("%.3f", time_roc_res$AUC[[3]])), color = "#E18727FF"
  ) +
  labs(x = "Specificity", y = "Sensitivity") +
  theme(
    axis.text = element_text(size = 10, color = "black"),
    axis.title.x = element_text(size = 12, color = "black", margin = margin(c(15, 0, 0, 0))),
    axis.title.y = element_text(size = 12, color = "black", margin = margin(c(0, 15, 0, 0)))
  )

prognosis.n144 <- data.frame(Tumor=rep("Neoadjuvant treated EAC",5),
                             Prognosis=c("Survival outcome","Risk stratification","1-year survival rate","3-year survival rate","5-year survival rate"),
                             AUC=c(gfit$auc,gfit1$auc,time_roc_res$AUC[[1]],time_roc_res$AUC[[2]],time_roc_res$AUC[[3]]))
prognosis.n144


############################## Chromophobe RCC (n=100) #####################
library(survival)
library(survminer)

setwd("C:/MyRdata8/Chromophobe RCC (n=100)")
table <- read.csv("Chromophobe RCC (n=100).csv")
dim(table)
table[1:5,1:10]

table.Chromophobe <- table[which(table$Subtype=="chromophobe"),]
dim(table.Chromophobe)
table.Chromophobe[1:5,1:15]

table1 <- table.Chromophobe[,c(1,8,9,13:2123)]
dim(table1)
table1[1:5,1:7]
rownames(table1) <- table1[,1]
table2 <- table1[,-1]
dim(table2)
table2[1:5,1:7]

table3 <- table2[,-c(1:2)]
dim(table3)
table3[1:5,1:7]
table3.n100 <- table3

# nucleotide profile data
table3 <- table3[,match(intersect(round(as.numeric(substr(colnames(table3),2,12)),4),nucleotide.annotation$query_mass),round(as.numeric(substr(colnames(table3),2,12)),4))]
dim(table3)

# median intensity on each mass
median1 <- apply(table3,MARGIN=2,FUN=median)
head(median1)
length(median1)
median1.n100 <- median1

# save nucleotide profile data
nucleotide.profile.Chromophobe.RCC.n100 <- table3
head(nucleotide.profile.Chromophobe.RCC.n100)


# calculate zero number for each m/z value
len1 <- c()
for (i in 1:ncol(table3)) {
  len <- length(which(table3[,i]==0))
  len1 <- c(len1,len)
}
len1
table4 <- table3[,which(len1 <= 90)]
dim(table4)

table5 <- cbind(table2[,1:2],table4)
dim(table5)
table5[1:5,1:7]

res.cut <- surv_cutpoint(table5, time = "OS", event = "STATOS",
                         minprop = 0.1,
                         variables = colnames(table5)[3:ncol(table5)])
res.cat <- surv_categorize(res.cut,labels = c(0, 1))
dim(res.cat)
res.cat$Gender <- table.Chromophobe$gender..0.female..1.male.
res.cat[1:5,1:7]

# save profile data with cutoff
res.cat.Chromophobe.RCC.n100 <- res.cat
dim(res.cat.Chromophobe.RCC.n100)
head(res.cat.Chromophobe.RCC.n100)


################ calculate p value, HR
result <- data.frame()
for (i in 3:ncol(res.cat)) {
  res.cox <- coxph(Surv(OS,STATOS)~res.cat[,i],data =res.cat)
  summ <- summary(res.cox)
  hypo <- cox.zph(res.cox)
  hypo1 <- hypo$table[-nrow(hypo$table),]
  result1 <- data.frame(mz.value=colnames(res.cat)[i],
                        #cutoff = res.cut[[i-2]][["estimate"]][["estimated cutpoint"]],
                        HR = signif(summ$coef[2], digits=3),
                        HR.confint.lower = signif(summ$conf.int[,"lower .95"], 3),
                        HR.confint.upper = signif(summ$conf.int[,"upper .95"],3),
                        logtest.p=summ$logtest[3],
                        waldtest.p=summ$waldtest[3],
                        sctest.p=summ$sctest[3],
                        hypo.p=hypo1[3])
  result <- rbind(result,result1)
}
head(result)
dim(result)
univariate.Chromophobe.RCC.n100 <- result

# screen by sctest p value
result2 <- result[result$sctest.p<0.05&result$hypo.p>0.05,]
dim(result2)
head(result2)
result3 <- result2
head(result3)
dim(result3)

# annotation on significant mass
result3$query_mass <- round(as.numeric(substr(result3$mz.value,2,12)),4)
result4 <- merge(result3,nucleotide.annotation)
result4$name.adduct <- paste(result4$adduct,result4$name,sep = '.')
result4$median.intensity <- median1[match(result4$mz.value,colnames(table3))]
dim(result4)
result5 <- result4[order(result4$ppm,decreasing = F),]
result6 <- result5[which(duplicated(result5$name.adduct)==F),]
dim(result6)
result7 <- result6[order(result6$median.intensity,decreasing = T),]
result8 <- result7[which(duplicated(result7$name)==F),]
result8 <- result8[order(result8$query_mass,decreasing = F),]
dim(result8)
head(result8)
univariate4.Chromophobe.RCC.n100 <- result8


############## validate by survival curve

# fit survival curve
fit <- survfit(Surv(OS,STATOS)~Gender,data =res.cat)
fit

survdiff(Surv(OS,STATOS)~Gender,data =res.cat) #采用log-rank 检验分析生存率差异

##绘制生存曲线##
ggsurvplot(fit, data = res.cat, conf.int = F, # 增加置信区间
           risk.table = TRUE, pval=T) #不同时间点风险人数表，增加P值




################# forest plot
library(forestplot)

dim(result8)
head(result8)
length(unique(result8$query_mass))
head(result2)

result10 <- data.frame(Feature=c(result8$query_mass),
                       Annotation=c(result8$name),
                       HR=c(result8$HR),
                       HR.lower=c(result8$HR.confint.lower),
                       HR.upper=c(result8$HR.confint.upper),
                       pvalue=round(c(result8$sctest.p),3))
dim(result10)
head(result10)

result10$HR1 <- paste(result10$HR,"(",result10$HR.lower,"-",result10$HR.upper,")",sep = '')
dim(result10)
head(result10)

labeltext <- as.matrix(rbind(c("Feature","Annotation","HR","pvalue"),result10[,c(1,2,7,6)]))
labeltext

forestplot(labeltext, 
           mean = c(NA,result10$HR),
           lower = c(NA,result10$HR.lower),
           upper = c(NA,result10$HR.upper),
           zero = 1,lwd.zero = 2,
           clip = c(0,10),
           graph.pos = 3,
           boxsize = 0.5,
           align = "l",
           xlog = FALSE,
           col=fpColors(line = "black", zero = "black",
                        box="black"),
           lty.ci = 7,
           lwd.ci = 1,
           ci.vertices.height = 0.15,
           lwd.xaxis=2,
           xlab="Hazard Ratio",
           txt_gp = fpTxtGp(ticks = gpar(cex = 0.9), xlab = gpar(cex = 0.8), cex = 0.8),
           lineheight = unit(8,"mm")) 





################################### multivariate cox analysis identify independent prognostic factor
dim(res.cat)
res.cat[1:5,1:10]

head(result2)
head(result8)

res.cat1 <- res.cat[,c("OS","STATOS",unique(result8$mz.value))]
dim(res.cat1)
res.cat1[1:5,1:5]

# save final profiling data after annotation
res.cat1.Chromophobe.RCC.n100 <- res.cat1


################## method1 feature selection from the beginning

### Step regression on Cox model
res.cox1 <- coxph(Surv(OS, STATOS) ~ .,data =res.cat1)
res.cox2 <- step(res.cox1)
ggforest(model=res.cox1,data=res.cat1)


###等比例风险假定
hypo <- cox.zph(res.cox1) 
hypo1 <- hypo$table[-nrow(hypo$table),]
dim(hypo1)
hypo2 <- hypo1[hypo1[,3]>0.05,]
dim(hypo2)
hypo2

res.cat2 <- res.cat1[,c("OS","STATOS",rownames(hypo2))]
dim(res.cat2)
res.cat2[1:5,1:5]
res.cox4 <- coxph(Surv(OS, STATOS) ~ .,data = res.cat2)
ggforest(model=res.cox4,data=res.cat2)


summ <- as.data.frame(cbind(summary(res.cox4)$coefficients[,c(2,5)],summary(res.cox4)$conf.int[,c(3,4)]))
summ$mz.value <- rownames(summ)
head(summ)
dim(summ)
multivariate.Chromophobe.RCC.n100.sum <- summ

summ1 <- summ
dim(summ1)
head(summ1)

summ2 <- merge(summ1,result8)
multivariate.Chromophobe.RCC.n100 <- summ2
dim(multivariate.Chromophobe.RCC.n100)
head(multivariate.Chromophobe.RCC.n100)
multivariate.Chromophobe.RCC.n100.sig <- summ2[summ2[,3]<0.05,]
dim(multivariate.Chromophobe.RCC.n100.sig)


#### forest plot
library(forestplot)

result12 <- summ2
dim(result12)
head(result12)


result14 <- data.frame(Feature=c(result12$query_mass),
                       Annotation=c(result12$name),
                       HR=c(result12[,2]),
                       HR.lower=c(result12[,4]),
                       HR.upper=c(result12[,5]),
                       pvalue=round(c(result12[,3]),3))
dim(result14)
head(result14)

result14$HR1 <- paste(round(result14$HR,2),"(",round(result14$HR.lower,2),"-",round(result14$HR.upper,2),")",sep = '')
result14$pvalue1 <- ifelse(result14$pvalue<0.001,"***",ifelse(result14$pvalue<0.01,"**",ifelse(result14$pvalue<0.05,"*"," ")))
result14$pvalue2 <- paste(result14$pvalue,result14$pvalue1,sep = '')
dim(result14)
head(result14)

labeltext <- as.matrix(rbind(c("Feature","Annotation","HR","pvalue"),result14[,c(1,2,7,9)]))
labeltext

forestplot(labeltext, 
           mean = c(NA,result14$HR),
           lower = c(NA,result14$HR.lower),
           upper = c(NA,result14$HR.upper),
           zero = 1,lwd.zero = 2,
           clip = c(0.001,40),
           graph.pos = 3,
           boxsize = 0.5,
           align = "l",
           xlog = T,
           col=fpColors(line = "black", zero = "black",
                        box="black"),
           lty.ci = 7,
           lwd.ci = 1,
           ci.vertices.height = 0.15,
           lwd.xaxis=2,
           xlab="Hazard Ratio",
           txt_gp = fpTxtGp(ticks = gpar(cex = 0.9), xlab = gpar(cex = 0.8), cex = 0.8),
           lineheight = unit(8,"mm")) 




############### use Cox model to predict survival outcome by ROC curve
res.cat3 <- res.cat1[,c("OS","STATOS",rownames(summ))]
dim(res.cat3)
res.cat3[1:5,1:5]

res.cox3 <- coxph(Surv(OS,STATOS)~.,data =res.cat3)
sum2 <- summary(res.cox3)
# C index export
sum2$concordance

# forest plot
ggforest(model = res.cox3,data=res.cat3)

library(pROC)
gfit <- roc(STATOS ~predict(res.cox3), data = res.cat3)

plot(gfit,
     print.auc=TRUE, #输出AUC值
     print.thres=TRUE, #输出cut-off值
     #main = "ROC Curve in Lung NAC tumor (n=77)", #设置图形的标题
     col= "red", #曲线颜色
     print.thres.col="black", #cut-off值字体的颜色
     identity.col="blue", #对角线颜色
     identity.lty=1,identity.lwd=1)


################# use Cox model to predict good/poor prognosis or low/high risk groups by ROC curve
dim(res.cat3)
res.cat3$STATOS1 <- ifelse(res.cat3$OS < median(res.cat3$OS),1,0) #1=short/poor survival or high risk,0=long/good survival or low risk  
res.cat3$STATOS1

gfit1 <- roc(STATOS1 ~predict(res.cox3), data = res.cat3)
cutoff <- gfit1$thresholds[which.max(gfit1$sensitivities+gfit1$specificities)]
cutoff
plot(gfit1,
     print.auc=TRUE, #输出AUC值
     print.thres=TRUE, #输出cut-off值
     #main = "ROC Curve in Lung NAC tumor (n=77)", #设置图形的标题
     col= "red", #曲线颜色
     print.thres.col="black", #cut-off值字体的颜色
     identity.col="blue", #对角线颜色
     identity.lty=1,identity.lwd=1)


### cutoff-based cox model to stratify patients into different prognostic risk groups by survival curve
res.cat3[1:5,1:ncol(res.cat3)]
res.cat3$risk <- ifelse(predict(res.cox3) <= median(predict(res.cox3)),'low','high')
# fit survival curve
fit <- survfit(Surv(OS,STATOS)~risk,data =res.cat3)
fit
##绘制生存曲线##
ggsurvplot(fit, data = res.cat3, conf.int = F, # 增加置信区间
           risk.table = TRUE, pval=T) #不同时间点风险人数表，增加P值



################# Nomogram绘制
library(car)
library(rms)
library(pROC)
library(timeROC)
library(ggDCA)

res.cat4 <- res.cat1[,c("OS","STATOS",rownames(summ))]
dim(res.cat4)
res.cat4[1:5,]

res.cox4 <- cph(Surv(OS,STATOS)~X211.00145+X344.03995+X384.03295,x=T,y=T,surv=T,data =res.cat4)
print(res.cox4)

ddist <- datadist(res.cat4)
options(datadist='ddist')

surv <- Survival(res.cox4) # 建立生存函数
surv

surv1 <- function(x) surv(12,x) # 1年OS
surv2 <- function(x) surv(36,x) # 3年OS
surv3 <- function(x) surv(60,x) # 5年OS

nom <- nomogram(res.cox4,
                fun = list(surv1,surv2,surv3),
                funlabel = c('1-year survival Probability',
                             '3-year survival Probability',
                             '5-year survival Probability'))
plot(nom, 
     lplabel="Linear Predictor",
     xfrac = 0.2, # 左侧标签距离坐标轴的距离
     #varname.label = TRUE, 
     tcl = -0.2, # 刻度长短和方向 
     lmgp = 0.1, # 坐标轴标签距离坐标轴远近
     points.label ='Points', 
     total.points.label = 'Total Points',
     cap.labels = FALSE,
     cex.var = 1, # 左侧标签字体大小
     cex.axis = 0.8, # 坐标轴字体大小
     col.grid = gray(c(0.8, 0.95))) # 竖线颜色


############# ROC curve on the prediction of 1/3/5-year survival rate
res.cat4 <- res.cat1[,c("OS","STATOS",rownames(summ))]
dim(res.cat4)
res.cat4[1:5,1:5]

res.cox4 <- cph(Surv(OS,STATOS)~X211.00145+X344.03995+X384.03295,x=T,y=T,surv=T,data =res.cat4)
print(res.cox4)

pred_f_training <- predict(res.cox4,res.cat4,type="lp")#!!!type="lp",是他没错
data_table <- data.frame(time=res.cat4[,"OS"],status=res.cat4[,"STATOS"],score=pred_f_training)


time_roc_res <- timeROC(
  T = data_table$time,
  delta = data_table$status,
  marker = data_table$score,
  cause = 1,
  weighting="marginal",
  times = c(12, 36, 60),
  ROC = TRUE,
  iid = TRUE
)


time_ROC_df <- data.frame(
  TP_1year = time_roc_res$TP[, 1],
  FP_1year = time_roc_res$FP[, 1],
  TP_3year = time_roc_res$TP[, 2],
  FP_3year = time_roc_res$FP[, 2],
  TP_5year = time_roc_res$TP[, 3],
  FP_5year = time_roc_res$FP[, 3]
)


ggplot(data = time_ROC_df) +
  geom_line(aes(x = FP_1year, y = TP_1year), linewidth = 1, color = "#BC3C29FF") +
  geom_line(aes(x = FP_3year, y = TP_3year), linewidth = 1, color = "#0072B5FF") +
  geom_line(aes(x = FP_5year, y = TP_5year), linewidth = 1, color = "#E18727FF") +
  geom_abline(slope = 1, intercept = 0, color = "grey", linewidth = 1, linetype = 2) +
  theme_bw() +
  annotate("text",
           x = 0.75, y = 0.25, size = 4.5,
           label = paste0("AUC at 1 year = ", sprintf("%.3f", time_roc_res$AUC[[1]])), color = "#BC3C29FF"
  ) +
  annotate("text",
           x = 0.75, y = 0.15, size = 4.5,
           label = paste0("AUC at 3 years = ", sprintf("%.3f", time_roc_res$AUC[[2]])), color = "#0072B5FF"
  ) +
  annotate("text",
           x = 0.75, y = 0.05, size = 4.5,
           label = paste0("AUC at 5 years = ", sprintf("%.3f", time_roc_res$AUC[[3]])), color = "#E18727FF"
  ) +
  labs(x = "Specificity", y = "Sensitivity") +
  theme(
    axis.text = element_text(size = 10, color = "black"),
    axis.title.x = element_text(size = 12, color = "black", margin = margin(c(15, 0, 0, 0))),
    axis.title.y = element_text(size = 12, color = "black", margin = margin(c(0, 15, 0, 0)))
  )

prognosis.n100 <- data.frame(Tumor=rep("Chromophobe RCC",5),
                             Prognosis=c("Survival outcome","Risk stratification","1-year survival rate","3-year survival rate","5-year survival rate"),
                             AUC=c(gfit$auc,gfit1$auc,time_roc_res$AUC[[1]],time_roc_res$AUC[[2]],time_roc_res$AUC[[3]]))
prognosis.n100


############################## Clear cell RCC (n=476) #####################
library(survival)
library(survminer)

setwd("C:/MyRdata8/Clear cell RCC (n=476)")
table <- read.csv("Clear cell RCC (n=476).csv")
dim(table)
table[1:5,1:10]

table.Clear.cell <- table[which(table$Subtype=="clear cell"),]
dim(table.Clear.cell)

table1 <- table.Clear.cell[,c(1,2,8,9,13:2123)]
dim(table1)
table1[1:5,1:7]
rownames(table1) <- table1[,1]
table2 <- table1[,-1]
dim(table2)
table2[1:5,1:7]

table3 <- table2[,-c(1:3)]
dim(table3)
table3[1:5,1:7]
table3.n476 <- table3

# nucleotide profile data
table3 <- table3[,match(intersect(round(as.numeric(substr(colnames(table3),2,12)),4),nucleotide.annotation$query_mass),round(as.numeric(substr(colnames(table3),2,12)),4))]
dim(table3)

# median intensity on each mass
median1 <- apply(table3,MARGIN=2,FUN=median)
head(median1)
length(median1)
median1.n476 <- median1

# save nucleotide profile data
nucleotide.profile.Clear.cell.RCC.n476 <- table3
dim(nucleotide.profile.Clear.cell.RCC.n476)


# calculate zero number for each m/z value
len1 <- c()
for (i in 1:ncol(table3)) {
  len <- length(which(table3[,i]==0))
  len1 <- c(len1,len)
}
len1
table4 <- table3[,which(len1 <= 435)]
dim(table4)

table2[1:5,1:7]
table5 <- cbind(table2[,2:3],table4)
dim(table5)
table5[1:5,1:7]

res.cut <- surv_cutpoint(table5, time = "OS", event = "STATOS",
                         minprop = 0.1,
                         variables = colnames(table5)[3:ncol(table5)])


res.cat.1 <- surv_categorize(res.cut,labels = c(0, 1))
dim(res.cat.1)
res.cat.1[1:5,1:7]

res.cat <- cbind(table2$gender..0.female..1.male.,res.cat.1)
dim(res.cat)
res.cat[1:5,1:10]
colnames(res.cat)[1] <- "Gender"
res.cat[1:5,1:10]

# save profile data with cutoff
res.cat.Clear.cell.RCC.n476 <- res.cat
dim(res.cat.Clear.cell.RCC.n476)
head(res.cat.Clear.cell.RCC.n476)


################ univariate cox analysis to calculate p value, HR
result <- data.frame()
for (i in c(1,4:ncol(res.cat))) {
  res.cox <- coxph(Surv(OS,STATOS)~res.cat[,i],data =res.cat)
  summ <- summary(res.cox)
  hypo <- cox.zph(res.cox)
  hypo1 <- hypo$table[-nrow(hypo$table),]
  result1 <- data.frame(mz.value=colnames(res.cat)[i],
                        HR = signif(summ$coef[2], digits=3),
                        HR.confint.lower = signif(summ$conf.int[,"lower .95"], 3),
                        HR.confint.upper = signif(summ$conf.int[,"upper .95"],3),
                        logtest.p=summ$logtest[3],
                        waldtest.p=summ$waldtest[3],
                        sctest.p=summ$sctest[3],
                        hypo.p=hypo1[3])
  result <- rbind(result,result1)
}
head(result)
dim(result)
univariate.Clear.cell.RCC.n476 <- result

# screen by sctest p value
result2 <- result[result$sctest.p<0.05&result$hypo.p>0.05,]
dim(result2)
head(result2)
result3 <- result2[-c(1),]
head(result3)
dim(result3)

# annotation on significant mass
result3$query_mass <- round(as.numeric(substr(result3$mz.value,2,12)),4)
result4 <- merge(result3,nucleotide.annotation)
result4$name.adduct <- paste(result4$adduct,result4$name,sep = '.')
result4$median.intensity <- median1[match(result4$mz.value,colnames(table3))]
dim(result4)
result5 <- result4[order(result4$ppm,decreasing = F),]
result6 <- result5[which(duplicated(result5$name.adduct)==F),]
dim(result6)
result7 <- result6[order(result6$median.intensity,decreasing = T),]
result8 <- result7[which(duplicated(result7$name)==F),]
result8 <- result8[order(result8$query_mass,decreasing = F),]
dim(result8)
head(result8)
univariate4.Clear.cell.RCC.n476 <- result8


############## validate by survival curve
# fit survival curve
fit <- survfit(Surv(OS,STATOS)~Gender ,data =res.cat)
fit
summary(fit)

survdiff(Surv(OS,STATOS)~Gender ,data =res.cat) #采用log-rank 检验分析生存率差异

##绘制生存曲线##
ggsurvplot(fit, data = res.cat, conf.int = F, # 增加置信区间
           risk.table = TRUE, pval=T) #不同时间点风险人数表，增加P值


################# forest plot
library(forestplot)

dim(result8)
head(result8)
length(unique(result8$query_mass))
head(result2)

result10 <- data.frame(Feature=c(result2$mz.value[1],result8$query_mass),
                       Annotation=c("Male:1",result8$name),
                       HR=c(result2$HR[1],result8$HR),
                       HR.lower=c(result2$HR.confint.lower[1],result8$HR.confint.lower),
                       HR.upper=c(result2$HR.confint.upper[1],result8$HR.confint.upper),
                       pvalue=round(c(result2$sctest.p[1],result8$sctest.p),3))
dim(result10)
head(result10)

result10$HR1 <- paste(result10$HR,"(",result10$HR.lower,"-",result10$HR.upper,")",sep = '')
dim(result10)
head(result10)

labeltext <- as.matrix(rbind(c("Feature","Annotation","HR","pvalue"),result10[,c(1,2,7,6)]))
labeltext

forestplot(labeltext, 
           mean = c(NA,result10$HR),
           lower = c(NA,result10$HR.lower),
           upper = c(NA,result10$HR.upper),
           zero = 1,lwd.zero = 2,
           clip = c(0.01,300),
           graph.pos = 3,
           boxsize = 0.5,
           align = "l",
           xlog = T,
           col=fpColors(line = "black", zero = "black",
                        box="black"),
           lty.ci = 7,
           lwd.ci = 1,
           ci.vertices.height = 0.15,
           lwd.xaxis=2,
           xlab="Hazard Ratio",
           txt_gp = fpTxtGp(ticks = gpar(cex = 0.9), xlab = gpar(cex = 0.8), cex = 0.7),
           lineheight = unit(8,"mm")) 



################################################ multivariate cox analysis identify independent prognostic factor
dim(res.cat)
res.cat[1:5,1:10]

head(result2)
head(result8)

res.cat1 <- res.cat[,c("OS","STATOS",result2$mz.value[1],unique(result8$mz.value))]
dim(res.cat1)
res.cat1[1:5,1:7]

# save final profiling data after annotation
res.cat1.Clear.cell.RCC.n476 <- res.cat1


################## method1 feature selection from the beginning

### Step regression on Cox model
res.cox1 <- coxph(Surv(OS, STATOS) ~ .,data =res.cat1)
res.cox2 <- step(res.cox1)
res.cox3 <- coxph(Surv(OS, STATOS) ~ Gender + X152.9957 + X211.00145 + X346.0547 + 
                    X362.0508667 + X426.0220333,data =res.cat1)
ggforest(model=res.cox3,data=res.cat1)


###等比例风险假定
hypo <- cox.zph(res.cox3) 
hypo1 <- hypo$table[-nrow(hypo$table),]
dim(hypo1)
hypo2 <- hypo1[hypo1[,3]>0.05,]
dim(hypo2)
hypo2

res.cat2 <- res.cat1[,c("OS","STATOS",rownames(hypo2))]
dim(res.cat2)
res.cat2[1:5,1:5]
res.cox4 <- coxph(Surv(OS, STATOS) ~ .,data = res.cat2)
ggforest(model=res.cox4,data=res.cat2)


summ <- as.data.frame(cbind(summary(res.cox4)$coefficients[,c(2,5)],summary(res.cox4)$conf.int[,c(3,4)]))
summ$mz.value <- rownames(summ)
head(summ)
dim(summ)
multivariate.Clear.cell.RCC.n476.sum <- summ

summ1 <- summ[-c(1),]
dim(summ1)
head(summ1)

summ2 <- merge(summ1,result8)
multivariate.Clear.cell.RCC.n476 <- summ2
dim(multivariate.Clear.cell.RCC.n476)
head(multivariate.Clear.cell.RCC.n476)
multivariate.Clear.cell.RCC.n476.sig <- summ2[summ2[,3]<0.05,]
dim(multivariate.Clear.cell.RCC.n476.sig)


#### forest plot
library(forestplot)

result12 <- summ2
dim(result12)
head(result12)

head(summ)
result13 <- summ[1,]
head(result13)

result14 <- data.frame(Feature=c(result13$mz.value,result12$query_mass),
                       Annotation=c("Male:1",result12$name),
                       HR=c(result13[,1],result12[,2]),
                       HR.lower=c(result13[,3],result12[,4]),
                       HR.upper=c(result13[,4],result12[,5]),
                       pvalue=round(c(result13[,2],result12[,3]),3))
dim(result14)
head(result14)

result14$HR1 <- paste(round(result14$HR,2),"(",round(result14$HR.lower,2),"-",round(result14$HR.upper,2),")",sep = '')
result14$pvalue1 <- ifelse(result14$pvalue<0.001,"***",ifelse(result14$pvalue<0.01,"**",ifelse(result14$pvalue<0.05,"*"," ")))
result14$pvalue2 <- paste(result14$pvalue,result14$pvalue1,sep = '')
dim(result14)
head(result14)

labeltext <- as.matrix(rbind(c("Feature","Annotation","HR","pvalue"),result14[,c(1,2,7,9)]))
labeltext

forestplot(labeltext, 
           mean = c(NA,result14$HR),
           lower = c(NA,result14$HR.lower),
           upper = c(NA,result14$HR.upper),
           zero = 1,lwd.zero = 2,
           clip = c(0.01,300),
           graph.pos = 3,
           boxsize = 0.5,
           align = "l",
           xlog = T,
           col=fpColors(line = "black", zero = "black",
                        box="black"),
           lty.ci = 7,
           lwd.ci = 1,
           ci.vertices.height = 0.15,
           lwd.xaxis=2,
           xlab="Hazard Ratio",
           txt_gp = fpTxtGp(ticks = gpar(cex = 0.9), xlab = gpar(cex = 0.8), cex = 0.8),
           lineheight = unit(8,"mm")) 


############### use Cox model to predict survival outcome by ROC curve
res.cat3 <- res.cat[,c("OS","STATOS",rownames(summ))]
dim(res.cat3)
res.cat3[1:5,1:5]

res.cox3 <- coxph(Surv(OS,STATOS)~.,data =res.cat3)
sum2 <- summary(res.cox3)
# C index export
sum2$concordance

# forest plot
ggforest(model = res.cox3,data=res.cat3)

library(pROC)
gfit <- roc(STATOS ~predict(res.cox3), data = res.cat3)

plot(gfit,
     print.auc=TRUE, #输出AUC值
     print.thres=TRUE, #输出cut-off值
     #main = "ROC Curve in Lung NAC tumor (n=77)", #设置图形的标题
     col= "red", #曲线颜色
     print.thres.col="black", #cut-off值字体的颜色
     identity.col="blue", #对角线颜色
     identity.lty=1,identity.lwd=1)


################# use Cox model to predict good/poor prognosis or low/high risk groups by ROC curve
res.cat3$STATOS1 <- ifelse(res.cat3$OS < median(res.cat3$OS),1,0) #1=short/poor survival or high risk,0=long/good survival or low risk  
res.cat3$STATOS1

gfit1 <- roc(STATOS1 ~ predict(res.cox3), data = res.cat3)
cutoff <- gfit1$thresholds[which.max(gfit1$sensitivities+gfit1$specificities)]
cutoff
plot(gfit1,
     print.auc=TRUE, #输出AUC值
     print.thres=TRUE, #输出cut-off值
     #main = "ROC Curve in Lung NAC tumor (n=77)", #设置图形的标题
     col= "red", #曲线颜色
     print.thres.col="black", #cut-off值字体的颜色
     identity.col="blue", #对角线颜色
     identity.lty=1,identity.lwd=1)

### cutoff-based cox model to stratify patients into different prognostic risk groups by survival curve
res.cat3[1:5,1:ncol(res.cat3)]
res.cat3$risk <- ifelse(predict(res.cox3) < cutoff,'low','high')
# fit survival curve
fit <- survfit(Surv(OS,STATOS)~risk,data =res.cat3)
fit
##绘制生存曲线##
ggsurvplot(fit, data = res.cat3, conf.int = F, # 增加置信区间
           risk.table = TRUE, pval=T) #不同时间点风险人数表，增加P值


################# Nomogram绘制
library(car)
library(rms)
library(pROC)
library(timeROC)
library(ggDCA)

res.cat4 <- res.cat1[,c("OS","STATOS",rownames(summ))]
dim(res.cat4)
res.cat4[1:5,]

res.cox4 <- cph(as.formula(paste0('Surv(OS,STATOS)~',paste(colnames(res.cat4)[-c(1:2)],sep = '',collapse = '+'))),
                x=T,y=T,surv=T,data =res.cat4)
print(res.cox4)

ddist <- datadist(res.cat4)
options(datadist='ddist')

surv <- Survival(res.cox4) # 建立生存函数
surv

surv1 <- function(x) surv(12,x) # 1年OS
surv2 <- function(x) surv(36,x) # 3年OS
surv3 <- function(x) surv(60,x) # 5年OS

nom <- nomogram(res.cox4,
                fun = list(surv1,surv2,surv3),
                funlabel = c('1-year survival Probability',
                             '3-year survival Probability',
                             '5-year survival Probability'))
plot(nom, 
     lplabel="Linear Predictor",
     xfrac = 0.2, # 左侧标签距离坐标轴的距离
     #varname.label = TRUE, 
     tcl = -0.2, # 刻度长短和方向 
     lmgp = 0.1, # 坐标轴标签距离坐标轴远近
     points.label ='Points', 
     total.points.label = 'Total Points',
     cap.labels = FALSE,
     cex.var = 0.7, # 左侧标签字体大小
     cex.axis = 0.6, # 坐标轴字体大小
     col.grid = gray(c(0.8, 0.95))) # 竖线颜色


############# ROC curve on the prediction of 1/3/5-year survival rate
res.cat4 <- res.cat1[,c("OS","STATOS",rownames(summ))]
dim(res.cat4)
res.cat4[1:5,1:5]

res.cox4 <- cph(as.formula(paste0('Surv(OS,STATOS)~',paste(colnames(res.cat4)[-c(1:2)],sep = '',collapse = '+'))),
                x=T,y=T,surv=T,data =res.cat4)
print(res.cox4)

pred_f_training <- predict(res.cox4,res.cat4,type="lp")#!!!type="lp",是他没错
data_table <- data.frame(time=res.cat4[,"OS"],status=res.cat4[,"STATOS"],score=pred_f_training)


time_roc_res <- timeROC(
  T = data_table$time,
  delta = data_table$status,
  marker = data_table$score,
  cause = 1,
  weighting="marginal",
  times = c(12, 36, 60),
  ROC = TRUE,
  iid = TRUE
)


time_ROC_df <- data.frame(
  TP_1year = time_roc_res$TP[, 1],
  FP_1year = time_roc_res$FP[, 1],
  TP_3year = time_roc_res$TP[, 2],
  FP_3year = time_roc_res$FP[, 2],
  TP_5year = time_roc_res$TP[, 3],
  FP_5year = time_roc_res$FP[, 3]
)


ggplot(data = time_ROC_df) +
  geom_line(aes(x = FP_1year, y = TP_1year), linewidth = 1, color = "#BC3C29FF") +
  geom_line(aes(x = FP_3year, y = TP_3year), linewidth = 1, color = "#0072B5FF") +
  geom_line(aes(x = FP_5year, y = TP_5year), linewidth = 1, color = "#E18727FF") +
  geom_abline(slope = 1, intercept = 0, color = "grey", linewidth = 1, linetype = 2) +
  theme_bw() +
  annotate("text",
           x = 0.75, y = 0.25, size = 4.5,
           label = paste0("AUC at 1 year = ", sprintf("%.3f", time_roc_res$AUC[[1]])), color = "#BC3C29FF"
  ) +
  annotate("text",
           x = 0.75, y = 0.15, size = 4.5,
           label = paste0("AUC at 3 years = ", sprintf("%.3f", time_roc_res$AUC[[2]])), color = "#0072B5FF"
  ) +
  annotate("text",
           x = 0.75, y = 0.05, size = 4.5,
           label = paste0("AUC at 5 years = ", sprintf("%.3f", time_roc_res$AUC[[3]])), color = "#E18727FF"
  ) +
  labs(x = "Specificity", y = "Sensitivity") +
  theme(
    axis.text = element_text(size = 10, color = "black"),
    axis.title.x = element_text(size = 12, color = "black", margin = margin(c(15, 0, 0, 0))),
    axis.title.y = element_text(size = 12, color = "black", margin = margin(c(0, 15, 0, 0)))
  )

prognosis.n476 <- data.frame(Tumor=rep("Clear cell RCC",5),
                             Prognosis=c("Survival outcome","Risk stratification","1-year survival rate","3-year survival rate","5-year survival rate"),
                             AUC=c(gfit$auc,gfit1$auc,time_roc_res$AUC[[1]],time_roc_res$AUC[[2]],time_roc_res$AUC[[3]]))
prognosis.n476



############################## Papillary RCC (n=101) #####################
library(survival)
library(survminer)

setwd("C:/MyRdata8/Papillary RCC (n=101)")
table <- read.csv("Papillary RCC (n=101).csv")
dim(table)
table[1:5,1:10]

table.Papillary <- table[which(table$Subtype=="papillary RCC"),]
dim(table.Papillary)
table.Papillary[1:5,1:10]

table1 <- table.Papillary[,c(1,8,9,13:2123)]
dim(table1)
table1[1:5,1:7]
rownames(table1) <- table1[,1]
table2 <- table1[,-1]
dim(table2)
table2[1:5,1:7]

table3 <- table2[,-c(1:2)]
dim(table3)
table3[1:5,1:7]
table3.n101 <- table3

# nucleotide profile data
table3 <- table3[,match(intersect(round(as.numeric(substr(colnames(table3),2,12)),4),nucleotide.annotation$query_mass),round(as.numeric(substr(colnames(table3),2,12)),4))]
dim(table3)

# median intensity on each mass
median1 <- apply(table3,MARGIN=2,FUN=median)
head(median1)
length(median1)
median1.n101 <- median1

# save nucleotide profile data
nucleotide.profile.Papillary.RCC.n101 <- table3
head(nucleotide.profile.Papillary.RCC.n101)


# calculate zero number for each m/z value
len1 <- c()
for (i in 1:ncol(table3)) {
  len <- length(which(table3[,i]==0))
  len1 <- c(len1,len)
}
len1
table4 <- table3[,which(len1 <= 91)]
dim(table4)

table5 <- cbind(table2[,1:2],table4)
dim(table5)
table5[1:5,1:7]

res.cut <- surv_cutpoint(table5, time = "OS", event = "STATOS",
                         minprop = 0.1,
                         variables = colnames(table5)[3:ncol(table5)])
res.cat <- surv_categorize(res.cut,labels = c(0, 1))
dim(res.cat)
res.cat[1:5,1:7]
res.cat$Gender <- table.Papillary$gender..0.female..1.male.
head(res.cat)

# save profile data with cutoff
res.cat.Papillary.RCC.n101 <- res.cat
dim(res.cat.Papillary.RCC.n101)
head(res.cat.Papillary.RCC.n101)


################ calculate p value, HR
result <- data.frame()
for (i in 3:ncol(res.cat)) {
  res.cox <- coxph(Surv(OS,STATOS)~res.cat[,i],data =res.cat)
  summ <- summary(res.cox)
  hypo <- cox.zph(res.cox)
  hypo1 <- hypo$table[-nrow(hypo$table),]
  result1 <- data.frame(mz.value=colnames(res.cat)[i],
                        #cutoff = res.cut[[i-2]][["estimate"]][["estimated cutpoint"]],
                        HR = signif(summ$coef[2], digits=3),
                        HR.confint.lower = signif(summ$conf.int[,"lower .95"], 3),
                        HR.confint.upper = signif(summ$conf.int[,"upper .95"],3),
                        logtest.p=summ$logtest[3],
                        waldtest.p=summ$waldtest[3],
                        sctest.p=summ$sctest[3],
                        hypo.p=hypo1[3])
  result <- rbind(result,result1)
}
head(result)
dim(result)
univariate.Papillary.RCC.n101 <- result


# screen by sctest p value
result2 <- result[result$sctest.p<0.05&result$hypo.p>0.05,]
dim(result2)
head(result2)
result3 <- result2[-c(5),]
head(result3)
dim(result3)

# annotation on significant mass
result3$query_mass <- round(as.numeric(substr(result3$mz.value,2,12)),4)
result4 <- merge(result3,nucleotide.annotation)
result4$name.adduct <- paste(result4$adduct,result4$name,sep = '.')
result4$median.intensity <- median1[match(result4$mz.value,colnames(table3))]
dim(result4)
result5 <- result4[order(result4$ppm,decreasing = F),]
result6 <- result5[which(duplicated(result5$name.adduct)==F),]
dim(result6)
result7 <- result6[order(result6$median.intensity,decreasing = T),]
result8 <- result7[which(duplicated(result7$name)==F),]
result8 <- result8[order(result8$query_mass,decreasing = F),]
dim(result8)
head(result8)
univariate4.Papillary.RCC.n101 <- result8



############## validate by survival curve

# fit survival curve
fit <- survfit(Surv(OS,STATOS)~Gender,data =res.cat)
survdiff(Surv(OS,STATOS)~Gender,data =res.cat) #采用log-rank 检验分析生存率差异

##绘制生存曲线##
ggsurvplot(fit, data = res.cat, conf.int = F, # 增加置信区间
           risk.table = TRUE, pval=T) #不同时间点风险人数表，增加P值



################# forest plot
library(forestplot)

dim(result8)
head(result8)
length(unique(result8$query_mass))
head(result2)

result10 <- data.frame(Feature=c(result8$query_mass),
                       Annotation=c(result8$name),
                       HR=c(result8$HR),
                       HR.lower=c(result8$HR.confint.lower),
                       HR.upper=c(result8$HR.confint.upper),
                       pvalue=round(c(result8$sctest.p),3))
dim(result10)
head(result10)

result10$HR1 <- paste(result10$HR,"(",result10$HR.lower,"-",result10$HR.upper,")",sep = '')
dim(result10)
head(result10)

labeltext <- as.matrix(rbind(c("Feature","Annotation","HR","pvalue"),result10[,c(1,2,7,6)]))
labeltext

forestplot(labeltext, 
           mean = c(NA,result10$HR),
           lower = c(NA,result10$HR.lower),
           upper = c(NA,result10$HR.upper),
           zero = 1,lwd.zero = 2,
           clip = c(0,6),
           graph.pos = 3,
           boxsize = 0.5,
           align = "l",
           xlog = FALSE,
           col=fpColors(line = "black", zero = "black",
                        box="black"),
           lty.ci = 7,
           lwd.ci = 1,
           ci.vertices.height = 0.15,
           lwd.xaxis=2,
           xlab="Hazard Ratio",
           txt_gp = fpTxtGp(ticks = gpar(cex = 0.9), xlab = gpar(cex = 0.8), cex = 0.7),
           lineheight = unit(8,"mm")) 



################################### multivariate cox analysis identify independent prognostic factor
dim(res.cat)
res.cat[1:5,1:10]

head(result2)
head(result8)

res.cat1 <- res.cat[,c("OS","STATOS",unique(result8$mz.value))]
dim(res.cat1)
res.cat1[1:5,1:5]

# save final profiling data after annotation
res.cat1.Papillary.RCC.n101 <- res.cat1


################## method1 feature selection from the beginning

### Step regression on Cox model
res.cox1 <- coxph(Surv(OS, STATOS) ~ .,data =res.cat1)
res.cox2 <- step(res.cox1)
res.cox3 <- coxph(Surv(OS, STATOS) ~ X346.0547 + X362.0508667 + X384.9842,data =res.cat1)
ggforest(model=res.cox3,data=res.cat1)


###等比例风险假定
hypo <- cox.zph(res.cox3) 
hypo1 <- hypo$table[-nrow(hypo$table),]
dim(hypo1)
hypo2 <- hypo1[hypo1[,3]>0.05,]
dim(hypo2)
hypo2

res.cat2 <- res.cat1[,c("OS","STATOS",rownames(hypo2))]
dim(res.cat2)
res.cat2[1:5,1:5]
res.cox4 <- coxph(Surv(OS, STATOS) ~ .,data = res.cat2)
ggforest(model=res.cox4,data=res.cat2)


summ <- as.data.frame(cbind(summary(res.cox4)$coefficients[,c(2,5)],summary(res.cox4)$conf.int[,c(3,4)]))
summ$mz.value <- rownames(summ)
head(summ)
dim(summ)
multivariate.Papillary.RCC.n101.sum <- summ

summ1 <- summ
dim(summ1)
head(summ1)

summ2 <- merge(summ1,result8)
multivariate.Papillary.RCC.n101 <- summ2
dim(multivariate.Papillary.RCC.n101)
head(multivariate.Papillary.RCC.n101)
multivariate.Papillary.RCC.n101.sig <- summ2[summ2[,3]<0.05,]
dim(multivariate.Papillary.RCC.n101.sig)


#### forest plot
library(forestplot)

result12 <- summ2
dim(result12)
head(result12)


result14 <- data.frame(Feature=c(result12$query_mass),
                       Annotation=c(result12$name),
                       HR=c(result12[,2]),
                       HR.lower=c(result12[,4]),
                       HR.upper=c(result12[,5]),
                       pvalue=round(c(result12[,3]),3))
dim(result14)
head(result14)

result14$HR1 <- paste(round(result14$HR,2),"(",round(result14$HR.lower,2),"-",round(result14$HR.upper,2),")",sep = '')
result14$pvalue1 <- ifelse(result14$pvalue<0.001,"***",ifelse(result14$pvalue<0.01,"**",ifelse(result14$pvalue<0.05,"*"," ")))
result14$pvalue2 <- paste(result14$pvalue,result14$pvalue1,sep = '')
dim(result14)
head(result14)

labeltext <- as.matrix(rbind(c("Feature","Annotation","HR","pvalue"),result14[,c(1,2,7,9)]))
labeltext

forestplot(labeltext, 
           mean = c(NA,result14$HR),
           lower = c(NA,result14$HR.lower),
           upper = c(NA,result14$HR.upper),
           zero = 1,lwd.zero = 2,
           clip = c(0.01,40),
           graph.pos = 3,
           boxsize = 0.5,
           align = "l",
           xlog = T,
           col=fpColors(line = "black", zero = "black",
                        box="black"),
           lty.ci = 7,
           lwd.ci = 1,
           ci.vertices.height = 0.15,
           lwd.xaxis=2,
           xlab="Hazard Ratio",
           txt_gp = fpTxtGp(ticks = gpar(cex = 0.9), xlab = gpar(cex = 0.8), cex = 0.8),
           lineheight = unit(8,"mm")) 



############### use Cox model to predict survival outcome by ROC curve
res.cat3 <- res.cat1[,c("OS","STATOS",rownames(summ))]
dim(res.cat3)
res.cat3[1:5,1:5]

res.cox3 <- coxph(Surv(OS,STATOS)~.,data =res.cat3)
sum2 <- summary(res.cox3)
# C index export
sum2$concordance

# forest plot
ggforest(model = res.cox3,data=res.cat3)

library(pROC)
gfit <- roc(STATOS ~predict(res.cox3), data = res.cat3)

plot(gfit,
     print.auc=TRUE, #输出AUC值
     print.thres=TRUE, #输出cut-off值
     #main = "ROC Curve in Lung NAC tumor (n=77)", #设置图形的标题
     col= "red", #曲线颜色
     print.thres.col="black", #cut-off值字体的颜色
     identity.col="blue", #对角线颜色
     identity.lty=1,identity.lwd=1)


################# use Cox model to predict good/poor prognosis or low/high risk groups by ROC curve
res.cat3$STATOS1 <- ifelse(res.cat3$OS < median(res.cat3$OS),1,0) #1=short/poor survival or high risk,0=long/good survival or low risk  
res.cat3$STATOS1

gfit1 <- roc(STATOS1 ~predict(res.cox3), data = res.cat3)
cutoff <- gfit1$thresholds[which.max(gfit1$sensitivities+gfit1$specificities)]
cutoff
plot(gfit1,
     print.auc=TRUE, #输出AUC值
     print.thres=TRUE, #输出cut-off值
     #main = "ROC Curve in Lung NAC tumor (n=77)", #设置图形的标题
     col= "red", #曲线颜色
     print.thres.col="black", #cut-off值字体的颜色
     identity.col="blue", #对角线颜色
     identity.lty=1,identity.lwd=1)

### cutoff-based cox model to stratify patients into different prognostic risk groups by survival curve
res.cat3[1:5,1:ncol(res.cat3)]
res.cat3$risk <- ifelse(predict(res.cox3) < cutoff,'low','high')
# fit survival curve
fit <- survfit(Surv(OS,STATOS)~risk,data =res.cat3)
fit
##绘制生存曲线##
ggsurvplot(fit, data = res.cat3, conf.int = F, # 增加置信区间
           risk.table = TRUE, pval=T) #不同时间点风险人数表，增加P值




################# Nomogram绘制

library(car)
library(rms)
library(pROC)
library(timeROC)
library(ggDCA)

res.cat4 <- res.cat1[,c("OS","STATOS",rownames(summ))]
dim(res.cat4)
res.cat4[1:5,]

res.cox4 <- cph(as.formula(paste0('Surv(OS,STATOS)~',paste(colnames(res.cat4)[-c(1:2)],sep = '',collapse = '+'))),
                x=T,y=T,surv=T,data =res.cat4)
print(res.cox4)

ddist <- datadist(res.cat4)
options(datadist='ddist')

surv <- Survival(res.cox4) # 建立生存函数
surv

surv1 <- function(x) surv(12,x) # 1年OS
surv2 <- function(x) surv(36,x) # 3年OS
surv3 <- function(x) surv(60,x) # 5年OS

nom <- nomogram(res.cox4,
                fun = list(surv1,surv2,surv3),
                funlabel = c('1-year survival Probability',
                             '3-year survival Probability',
                             '5-year survival Probability'))
plot(nom, 
     lplabel="Linear Predictor",
     xfrac = 0.2, # 左侧标签距离坐标轴的距离
     #varname.label = TRUE, 
     tcl = -0.2, # 刻度长短和方向 
     lmgp = 0.1, # 坐标轴标签距离坐标轴远近
     points.label ='Points', 
     total.points.label = 'Total Points',
     cap.labels = FALSE,
     cex.var = 1, # 左侧标签字体大小
     cex.axis = 0.8, # 坐标轴字体大小
     col.grid = gray(c(0.8, 0.95))) # 竖线颜色


############# ROC curve on the prediction of 1/3/5-year survival rate
res.cat4 <- res.cat1[,c("OS","STATOS",rownames(summ))]
dim(res.cat4)
res.cat4[1:5,1:5]

res.cox4 <- cph(as.formula(paste0('Surv(OS,STATOS)~',paste(colnames(res.cat4)[-c(1:2)],sep = '',collapse = '+'))),
                x=T,y=T,surv=T,data =res.cat4)
print(res.cox4)

pred_f_training <- predict(res.cox4,res.cat4,type="lp")#!!!type="lp",是他没错
data_table <- data.frame(time=res.cat4[,"OS"],status=res.cat4[,"STATOS"],score=pred_f_training)


time_roc_res <- timeROC(
  T = data_table$time,
  delta = data_table$status,
  marker = data_table$score,
  cause = 1,
  weighting="marginal",
  times = c(12, 36, 60),
  ROC = TRUE,
  iid = TRUE
)


time_ROC_df <- data.frame(
  TP_1year = time_roc_res$TP[, 1],
  FP_1year = time_roc_res$FP[, 1],
  TP_3year = time_roc_res$TP[, 2],
  FP_3year = time_roc_res$FP[, 2],
  TP_5year = time_roc_res$TP[, 3],
  FP_5year = time_roc_res$FP[, 3]
)


ggplot(data = time_ROC_df) +
  geom_line(aes(x = FP_1year, y = TP_1year), linewidth = 1, color = "#BC3C29FF") +
  geom_line(aes(x = FP_3year, y = TP_3year), linewidth = 1, color = "#0072B5FF") +
  geom_line(aes(x = FP_5year, y = TP_5year), linewidth = 1, color = "#E18727FF") +
  geom_abline(slope = 1, intercept = 0, color = "grey", linewidth = 1, linetype = 2) +
  theme_bw() +
  annotate("text",
           x = 0.75, y = 0.25, size = 4.5,
           label = paste0("AUC at 1 year = ", sprintf("%.3f", time_roc_res$AUC[[1]])), color = "#BC3C29FF"
  ) +
  annotate("text",
           x = 0.75, y = 0.15, size = 4.5,
           label = paste0("AUC at 3 years = ", sprintf("%.3f", time_roc_res$AUC[[2]])), color = "#0072B5FF"
  ) +
  annotate("text",
           x = 0.75, y = 0.05, size = 4.5,
           label = paste0("AUC at 5 years = ", sprintf("%.3f", time_roc_res$AUC[[3]])), color = "#E18727FF"
  ) +
  labs(x = "Specificity", y = "Sensitivity") +
  theme(
    axis.text = element_text(size = 10, color = "black"),
    axis.title.x = element_text(size = 12, color = "black", margin = margin(c(15, 0, 0, 0))),
    axis.title.y = element_text(size = 12, color = "black", margin = margin(c(0, 15, 0, 0)))
  )


prognosis.n101 <- data.frame(Tumor=rep("Papillary RCC",5),
                             Prognosis=c("Survival outcome","Risk stratification","1-year survival rate","3-year survival rate","5-year survival rate"),
                             AUC=c(gfit$auc,gfit1$auc,time_roc_res$AUC[[1]],time_roc_res$AUC[[2]],time_roc_res$AUC[[3]]))
prognosis.n101




############################## Adrenocortical Carcinoma (ACC) (n=72) #####################
library(survival)
library(survminer)

setwd("C:/MyRdata8/Adrenocortical Carcinoma (ACC) (n=72)")
table1 <- read.csv("Adrenocortical Carcinoma (ACC) (n=72).csv")
dim(table1)
table1[1:5,1:7]
rownames(table1) <- table1[,1]
table2 <- table1[,-1]
dim(table2)
table2[1:5,1:7]

table3 <- table2[,-c(1:3)]
dim(table3)
table3[1:5,1:7]
table3.n72 <- table3

# nucleotide profile data
table3 <- table3[,match(intersect(round(as.numeric(substr(colnames(table3),2,12)),4),nucleotide.annotation$query_mass),round(as.numeric(substr(colnames(table3),2,12)),4))]
dim(table3)

# median intensity on each mass
median1 <- apply(table3,MARGIN=2,FUN=median)
head(median1)
length(median1)
median1.n72 <- median1

# save nucleotide profile data
nucleotide.profile.ACC.n72 <- table3
head(nucleotide.profile.ACC.n72)


# calculate zero number for each m/z value
len1 <- c()
for (i in 1:ncol(table3)) {
  len <- length(which(table3[,i]==0))
  len1 <- c(len1,len)
}
len1
table4 <- table3[,which(len1 <= 72)]
dim(table4)

table5 <- cbind(table2[,c(2:3)],table4)
dim(table5)
table5[1:5,1:7]

res.cut <- surv_cutpoint(table5, time = "OS", event = "STATOS",
                         minprop = 0.1,
                         variables = colnames(table5)[3:ncol(table5)])
res.cat.1 <- surv_categorize(res.cut,labels = c(0, 1))
dim(res.cat.1)
res.cat.1[1:5,1:7]

res.cat <- cbind(table2[,1],res.cat.1)
dim(res.cat)
res.cat[1:5,1:10]
colnames(res.cat)[1] <- "UICC"
res.cat <- res.cat[,c(2,3,1,4:ncol(res.cat))]


# save profile data with cutoff
res.cat.ACC.n72 <- res.cat
dim(res.cat.ACC.n72)
head(res.cat.ACC.n72)


################ univariate cox analysis to calculate p value, HR
result <- data.frame()
for (i in 3:ncol(res.cat)) {
  res.cox <- coxph(Surv(OS,STATOS)~res.cat[,i],data =res.cat)
  summ <- summary(res.cox)
  hypo <- cox.zph(res.cox)
  hypo1 <- hypo$table[-nrow(hypo$table),]
  result1 <- data.frame(mz.value=colnames(res.cat)[i],
                        #cutoff = res.cut[[i-2]][["estimate"]][["estimated cutpoint"]],
                        HR = signif(summ$coef[2], digits=3),
                        HR.confint.lower = signif(summ$conf.int[,"lower .95"], 3),
                        HR.confint.upper = signif(summ$conf.int[,"upper .95"],3),
                        logtest.p=summ$logtest[3],
                        waldtest.p=summ$waldtest[3],
                        sctest.p=summ$sctest[3],
                        hypo.p=hypo1[3])
  result <- rbind(result,result1)
}
head(result)
dim(result)
univariate.ACC.n72 <- result


# screen by sctest p value
result2 <- result[result$sctest.p<0.05&result$hypo.p>0.05,]
dim(result2)
head(result2)
result3 <- result2[-1,]
head(result3)
dim(result3)

# annotation on significant mass
result3$query_mass <- round(as.numeric(substr(result3$mz.value,2,12)),4)
result4 <- merge(result3,nucleotide.annotation)
result4$name.adduct <- paste(result4$adduct,result4$name,sep = '.')
result4$median.intensity <- median1[match(result4$mz.value,colnames(table3))]
dim(result4)
result5 <- result4[order(result4$ppm,decreasing = F),]
result6 <- result5[which(duplicated(result5$name.adduct)==F),]
dim(result6)
result7 <- result6[order(result6$median.intensity,decreasing = T),]
result8 <- result7[which(duplicated(result7$name)==F),]
result8 <- result8[order(result8$query_mass,decreasing = F),]
dim(result8)
head(result8)
univariate4.ACC.n72 <- result8


############## validate by survival curve
# fit survival curve
fit <- survfit(Surv(OS,STATOS)~res.cat[,3],data =res.cat)
survdiff(Surv(OS,STATOS)~res.cat[,3],data =res.cat) #采用log-rank 检验分析生存率差异

##绘制生存曲线##
ggsurvplot(fit, data = res.cat, conf.int = TRUE, # 增加置信区间
           risk.table = TRUE, pval=T) #不同时间点风险人数表，增加P值



################# forest plot
library(forestplot)

dim(result8)
head(result8)
length(unique(result8$query_mass))
head(result2)

result10 <- data.frame(Feature=c(result2$mz.value[1],result8$query_mass),
                       Annotation=c("1/2/3/4",result8$name),
                       HR=c(result2$HR[1],result8$HR),
                       HR.lower=c(result2$HR.confint.lower[1],result8$HR.confint.lower),
                       HR.upper=c(result2$HR.confint.upper[1],result8$HR.confint.upper),
                       pvalue=round(c(result2$sctest.p[1],result8$sctest.p),3))
dim(result10)
head(result10)

result10$HR1 <- paste(result10$HR,"(",result10$HR.lower,"-",result10$HR.upper,")",sep = '')
dim(result10)
head(result10)

labeltext <- as.matrix(rbind(c("Feature","Annotation","HR","pvalue"),result10[,c(1,2,7,6)]))
labeltext

forestplot(labeltext, 
           mean = c(NA,result10$HR),
           lower = c(NA,result10$HR.lower),
           upper = c(NA,result10$HR.upper),
           zero = 1,lwd.zero = 2,
           clip = c(0,6),
           graph.pos = 3,
           boxsize = 0.5,
           align = "l",
           xlog = FALSE,
           col=fpColors(line = "black", zero = "black",
                        box="black"),
           lty.ci = 7,
           lwd.ci = 1,
           ci.vertices.height = 0.15,
           lwd.xaxis=2,
           xlab="Hazard Ratio",
           txt_gp = fpTxtGp(ticks = gpar(cex = 0.9), xlab = gpar(cex = 0.8), cex = 0.8),
           lineheight = unit(8,"mm")) 





################################### multivariate cox analysis identify independent prognostic factor
dim(res.cat)
res.cat[1:5,1:10]

head(result2)
head(result8)

res.cat1 <- res.cat[,c("OS","STATOS","UICC",unique(result8$mz.value))]
dim(res.cat1)
res.cat1[1:5,1:5]

# save final profiling data after annotation
res.cat1.ACC.n72 <- res.cat1


################## method1 feature selection from the beginning

### Step regression on Cox model
res.cox1 <- coxph(Surv(OS, STATOS) ~ .,data =res.cat1)
res.cox2 <- step(res.cox1)
res.cox3 <- coxph(Surv(OS, STATOS) ~ UICC + X150.042 + X195.0065 + X273.0595 + 
                    X348.0355 + X408.0127 + X421.9955 + X424.9794,data =res.cat1)
ggforest(model=res.cox3,data=res.cat1)


###等比例风险假定
hypo <- cox.zph(res.cox3) 
hypo1 <- hypo$table[-nrow(hypo$table),]
dim(hypo1)
hypo2 <- hypo1[hypo1[,3]>0.05,]
dim(hypo2)
hypo2

res.cat2 <- res.cat1[,c("OS","STATOS",rownames(hypo2))]
dim(res.cat2)
res.cat2[1:5,1:5]
res.cox4 <- coxph(Surv(OS, STATOS) ~ .,data = res.cat2)
ggforest(model=res.cox4,data=res.cat2)


summ <- as.data.frame(cbind(summary(res.cox4)$coefficients[,c(2,5)],summary(res.cox4)$conf.int[,c(3,4)]))
summ$mz.value <- rownames(summ)
head(summ)
dim(summ)
multivariate.ACC.n72.sum <- summ

summ1 <- summ[-1,]
dim(summ1)
head(summ1)

summ2 <- merge(summ1,result8)
multivariate.ACC.n72 <- summ2
dim(multivariate.ACC.n72)
head(multivariate.ACC.n72)
multivariate.ACC.n72.sig <- summ2[summ2[,3]<0.05,]
dim(multivariate.ACC.n72.sig)


#### forest plot
library(forestplot)

result12 <- summ2
dim(result12)
head(result12)


result14 <- data.frame(Feature=c(rownames(summ)[1],result12$query_mass),
                       Annotation=c("1/2/3/4",result12$name),
                       HR=c(summ$`exp(coef)`[1],result12[,2]),
                       HR.lower=c(summ$`lower .95`[1],result12[,4]),
                       HR.upper=c(summ$`upper .95`[1],result12[,5]),
                       pvalue=round(c(summ$`Pr(>|z|)`[1],result12[,3]),3))
dim(result14)
head(result14)

result14$HR1 <- paste(round(result14$HR,2),"(",round(result14$HR.lower,2),"-",round(result14$HR.upper,2),")",sep = '')
result14$pvalue1 <- ifelse(result14$pvalue<0.001,"***",ifelse(result14$pvalue<0.01,"**",ifelse(result14$pvalue<0.05,"*"," ")))
result14$pvalue2 <- paste(result14$pvalue,result14$pvalue1,sep = '')
dim(result14)
head(result14)

labeltext <- as.matrix(rbind(c("Feature","Annotation","HR","pvalue"),result14[,c(1,2,7,9)]))
labeltext

forestplot(labeltext, 
           mean = c(NA,result14$HR),
           lower = c(NA,result14$HR.lower),
           upper = c(NA,result14$HR.upper),
           zero = 1,lwd.zero = 2,
           clip = c(0.01,50),
           graph.pos = 3,
           boxsize = 0.5,
           align = "l",
           xlog = T,
           col=fpColors(line = "black", zero = "black",
                        box="black"),
           lty.ci = 7,
           lwd.ci = 1,
           ci.vertices.height = 0.15,
           lwd.xaxis=2,
           xlab="Hazard Ratio",
           txt_gp = fpTxtGp(ticks = gpar(cex = 0.9), xlab = gpar(cex = 0.8), cex = 0.8),
           lineheight = unit(8,"mm")) 


############### use Cox model to predict survival outcome by ROC curve
res.cat3 <- res.cat[,c("OS","STATOS",rownames(summ))]
dim(res.cat3)
res.cat3[1:5,1:5]

res.cox3 <- coxph(Surv(OS,STATOS)~.,data =res.cat3)
sum2 <- summary(res.cox3)
# C index export
sum2$concordance

# forest plot
ggforest(model = res.cox3,data=res.cat3)

library(pROC)
gfit <- roc(STATOS ~predict(res.cox3), data = res.cat3)

plot(gfit,
     print.auc=TRUE, #输出AUC值
     print.thres=TRUE, #输出cut-off值
     #main = "ROC Curve in Lung NAC tumor (n=77)", #设置图形的标题
     col= "red", #曲线颜色
     print.thres.col="black", #cut-off值字体的颜色
     identity.col="blue", #对角线颜色
     identity.lty=1,identity.lwd=1)


################# use Cox model to predict good/poor prognosis or low/high risk groups by ROC curve
res.cat3$STATOS1 <- ifelse(res.cat3$OS < median(res.cat3$OS),1,0) #1=short/poor survival or high risk,0=long/good survival or low risk  
res.cat3$STATOS1

gfit1 <- roc(STATOS1 ~predict(res.cox3), data = res.cat3)
cutoff <- gfit1$thresholds[which.max(gfit1$sensitivities+gfit1$specificities)]
cutoff
plot(gfit1,
     print.auc=TRUE, #输出AUC值
     print.thres=TRUE, #输出cut-off值
     #main = "ROC Curve in Lung NAC tumor (n=77)", #设置图形的标题
     col= "red", #曲线颜色
     print.thres.col="black", #cut-off值字体的颜色
     identity.col="blue", #对角线颜色
     identity.lty=1,identity.lwd=1)


### cutoff-based cox model to stratify patients into different prognostic risk groups by survival curve
res.cat3[1:5,1:ncol(res.cat3)]
res.cat3$risk <- ifelse(predict(res.cox3) < cutoff,'low','high')
# fit survival curve
fit <- survfit(Surv(OS,STATOS)~risk,data =res.cat3)
fit
##绘制生存曲线##
ggsurvplot(fit, data = res.cat3, conf.int = F, # 增加置信区间
           risk.table = TRUE, pval=T) #不同时间点风险人数表，增加P值




################# Nomogram绘制

library(car)
library(rms)
library(pROC)
library(timeROC)
library(ggDCA)

res.cat4 <- res.cat1[,c("OS","STATOS",rownames(summ))]
dim(res.cat4)
res.cat4[1:5,]

res.cox4 <- cph(as.formula(paste0('Surv(OS,STATOS)~',paste(colnames(res.cat4)[-c(1:2)],sep = '',collapse = '+'))),
                x=T,y=T,surv=T,data =res.cat4)
print(res.cox4)

ddist <- datadist(res.cat4)
options(datadist='ddist')

surv <- Survival(res.cox4) # 建立生存函数
surv

surv1 <- function(x) surv(12,x) # 1年OS
surv2 <- function(x) surv(36,x) # 3年OS
surv3 <- function(x) surv(60,x) # 5年OS

nom <- nomogram(res.cox4,
                fun = list(surv1,surv2,surv3),
                funlabel = c('1-year survival Probability',
                             '3-year survival Probability',
                             '5-year survival Probability'))
plot(nom, 
     lplabel="Linear Predictor",
     xfrac = 0.2, # 左侧标签距离坐标轴的距离
     #varname.label = TRUE, 
     tcl = -0.2, # 刻度长短和方向 
     lmgp = 0.1, # 坐标轴标签距离坐标轴远近
     points.label ='Points', 
     total.points.label = 'Total Points',
     cap.labels = FALSE,
     cex.var = 0.8, # 左侧标签字体大小
     cex.axis = 0.4, # 坐标轴字体大小
     col.grid = gray(c(0.8, 0.95))) # 竖线颜色


############# ROC curve on the prediction of 1/3/5-year survival rate
res.cat4 <- res.cat1[,c("OS","STATOS",rownames(summ))]
dim(res.cat4)
res.cat4[1:5,1:5]

res.cox4 <- cph(as.formula(paste0('Surv(OS,STATOS)~',paste(colnames(res.cat4)[-c(1:2)],sep = '',collapse = '+'))),
                x=T,y=T,surv=T,data =res.cat4)
print(res.cox4)

pred_f_training <- predict(res.cox4,res.cat4,type="lp")#!!!type="lp",是他没错
data_table <- data.frame(time=res.cat4[,"OS"],status=res.cat4[,"STATOS"],score=pred_f_training)


time_roc_res <- timeROC(
  T = data_table$time,
  delta = data_table$status,
  marker = data_table$score,
  cause = 1,
  weighting="marginal",
  times = c(12, 36, 60),
  ROC = TRUE,
  iid = TRUE
)


time_ROC_df <- data.frame(
  TP_1year = time_roc_res$TP[, 1],
  FP_1year = time_roc_res$FP[, 1],
  TP_3year = time_roc_res$TP[, 2],
  FP_3year = time_roc_res$FP[, 2],
  TP_5year = time_roc_res$TP[, 3],
  FP_5year = time_roc_res$FP[, 3]
)


ggplot(data = time_ROC_df) +
  geom_line(aes(x = FP_1year, y = TP_1year), linewidth = 1, color = "#BC3C29FF") +
  geom_line(aes(x = FP_3year, y = TP_3year), linewidth = 1, color = "#0072B5FF") +
  geom_line(aes(x = FP_5year, y = TP_5year), linewidth = 1, color = "#E18727FF") +
  geom_abline(slope = 1, intercept = 0, color = "grey", linewidth = 1, linetype = 2) +
  theme_bw() +
  annotate("text",
           x = 0.75, y = 0.25, size = 4.5,
           label = paste0("AUC at 1 year = ", sprintf("%.3f", time_roc_res$AUC[[1]])), color = "#BC3C29FF"
  ) +
  annotate("text",
           x = 0.75, y = 0.15, size = 4.5,
           label = paste0("AUC at 3 years = ", sprintf("%.3f", time_roc_res$AUC[[2]])), color = "#0072B5FF"
  ) +
  annotate("text",
           x = 0.75, y = 0.05, size = 4.5,
           label = paste0("AUC at 5 years = ", sprintf("%.3f", time_roc_res$AUC[[3]])), color = "#E18727FF"
  ) +
  labs(x = "Specificity", y = "Sensitivity") +
  theme(
    axis.text = element_text(size = 10, color = "black"),
    axis.title.x = element_text(size = 12, color = "black", margin = margin(c(15, 0, 0, 0))),
    axis.title.y = element_text(size = 12, color = "black", margin = margin(c(0, 15, 0, 0)))
  )


prognosis.n72 <- data.frame(Tumor=rep("Adrenocortical Carcinoma",5),
                            Prognosis=c("Survival outcome","Risk stratification","1-year survival rate","3-year survival rate","5-year survival rate"),
                            AUC=c(gfit$auc,gfit1$auc,time_roc_res$AUC[[1]],time_roc_res$AUC[[2]],time_roc_res$AUC[[3]]))
prognosis.n72





############################## Pancreas tumor (n=107) #####################
library(survival)
library(survminer)

setwd("C:/MyRdata8/Pancreas tumor (n=107)")
table <- read.csv("Pancreas tumor (n=107).csv")
dim(table)
table[1:5,1:10]

table1 <- table[,c(1,2,3,8,9,14,15,17:5253)]
dim(table1)
table1[1:5,1:10]
rownames(table1) <- table1[,1]
table2 <- table1[,-1]
dim(table2)
table2[1:5,1:10]

table3 <- table2[,-c(1:6)]
dim(table3)
table3[1:5,1:7]
table3.n107 <- table3

# nucleotide profile data
table3 <- table3[,match(intersect(round(as.numeric(substr(colnames(table3),2,12)),4),nucleotide.annotation$query_mass),round(as.numeric(substr(colnames(table3),2,12)),4))]
dim(table3)


# median intensity on each mass
median1 <- apply(table3,MARGIN=2,FUN=median)
head(median1)
length(median1)
median1.n107 <- median1

# save nucleotide profile data
nucleotide.profile.Pancreas.tumor.n107 <- table3
head(nucleotide.profile.Pancreas.tumor.n107)


# calculate zero number for each m/z value
len1 <- c()
for (i in 1:ncol(table3)) {
  len <- length(which(table3[,i]==0))
  len1 <- c(len1,len)
}
len1
table4 <- table3[,which(len1 <= 96)]
dim(table4)

table2[1:5,1:10]
table5 <- cbind(table2[,5:6],table4)
dim(table5)
table5[1:5,1:7]

res.cut <- surv_cutpoint(table5, time = "OS", event = "STATOS",
                         minprop = 0.1,
                         variables = colnames(table5)[3:ncol(table5)])


res.cat.1 <- surv_categorize(res.cut,labels = c(0, 1))
dim(res.cat.1)
res.cat.1[1:5,1:7]

table2[1:5,1:10]
res.cat <- cbind(table2[,c(1:3)],res.cat.1)
dim(res.cat)
res.cat[1:5,1:10]

colnames(res.cat)[c(1,2)] <- c("Gender","Age")
res.cat[1:5,1:10]


# save profile data with cutoff
res.cat.Pancreas.tumor.n107 <- res.cat
dim(res.cat.Pancreas.tumor.n107)
head(res.cat.Pancreas.tumor.n107)


################ univariate cox analysis to calculate p value, HR
result <- data.frame()
for (i in c(1:3,6:ncol(res.cat))) {
  res.cox <- coxph(Surv(OS,STATOS)~res.cat[,i],data =res.cat)
  summ <- summary(res.cox)
  hypo <- cox.zph(res.cox)
  hypo1 <- hypo$table[-nrow(hypo$table),]
  result1 <- data.frame(mz.value=colnames(res.cat)[i],
                        HR = signif(summ$coef[2], digits=3),
                        HR.confint.lower = signif(summ$conf.int[,"lower .95"], 3),
                        HR.confint.upper = signif(summ$conf.int[,"upper .95"],3),
                        logtest.p=summ$logtest[3],
                        waldtest.p=summ$waldtest[3],
                        sctest.p=summ$sctest[3],
                        hypo.p=hypo1[3])
  result <- rbind(result,result1)
}
head(result)
dim(result)
univariate.Pancreas.tumor.n107 <- result


# screen by sctest p value
result2 <- result[result$sctest.p<0.05&result$hypo.p>0.05,]
dim(result2)
head(result2)
result3 <- result2[-c(1),]
head(result3)
dim(result3)

# annotation on significant mass
result3$query_mass <- round(as.numeric(substr(result3$mz.value,2,12)),4)
result4 <- merge(result3,nucleotide.annotation)
result4$name.adduct <- paste(result4$adduct,result4$name,sep = '.')
result4$median.intensity <- median1[match(result4$mz.value,colnames(table3))]
dim(result4)
result5 <- result4[order(result4$ppm,decreasing = F),]
result6 <- result5[which(duplicated(result5$name.adduct)==F),]
dim(result6)
result7 <- result6[order(result6$median.intensity,decreasing = T),]
result8 <- result7[which(duplicated(result7$name)==F),]
result8 <- result8[order(result8$query_mass,decreasing = F),]
dim(result8)
head(result8)
univariate4.Pancreas.tumor.n107 <- result8


############## validate by survival curve
# fit survival curve
fit <- survfit(Surv(OS,STATOS)~UICC,data =res.cat)
fit
survdiff(Surv(OS,STATOS)~UICC,data =res.cat) #采用log-rank 检验分析生存率差异

##绘制生存曲线##
ggsurvplot(fit, data = res.cat, conf.int = F, # 增加置信区间
           risk.table = TRUE, pval=T) #不同时间点风险人数表，增加P值





################# forest plot
library(forestplot)

dim(result8)
head(result8)
length(unique(result8$query_mass))
head(result2)

result10 <- data.frame(Feature=c(result2$mz.value[1],result8$query_mass),
                       Annotation=c("1/2/3/4",result8$name),
                       HR=c(result2$HR[1],result8$HR),
                       HR.lower=c(result2$HR.confint.lower[1],result8$HR.confint.lower),
                       HR.upper=c(result2$HR.confint.upper[1],result8$HR.confint.upper),
                       pvalue=round(c(result2$sctest.p[1],result8$sctest.p),3))
dim(result10)
head(result10)

result10$HR1 <- paste(result10$HR,"(",result10$HR.lower,"-",result10$HR.upper,")",sep = '')
dim(result10)
head(result10)

labeltext <- as.matrix(rbind(c("Feature","Annotation","HR","pvalue"),result10[,c(1,2,7,6)]))
labeltext

forestplot(labeltext, 
           mean = c(NA,result10$HR),
           lower = c(NA,result10$HR.lower),
           upper = c(NA,result10$HR.upper),
           zero = 1,lwd.zero = 2,
           clip = c(0,6),
           graph.pos = 3,
           boxsize = 0.5,
           align = "l",
           xlog = FALSE,
           col=fpColors(line = "black", zero = "black",
                        box="black"),
           lty.ci = 7,
           lwd.ci = 1,
           ci.vertices.height = 0.15,
           lwd.xaxis=2,
           xlab="Hazard Ratio",
           txt_gp = fpTxtGp(ticks = gpar(cex = 0.9), xlab = gpar(cex = 0.8), cex = 0.7),
           lineheight = unit(8,"mm")) 





################################################ multivariate cox analysis identify independent prognostic factor
dim(res.cat)
res.cat[1:5,1:10]

head(result2)
head(result8)

res.cat1 <- res.cat[,c("OS","STATOS",result2$mz.value[1],unique(result8$mz.value))]
dim(res.cat1)
res.cat1[1:5,1:7]

# save final profiling data after annotation
res.cat1.Pancreas.tumor.n107 <- res.cat1


################## method1 feature selection from the beginning

### Step regression on Cox model
res.cox1 <- coxph(Surv(OS, STATOS) ~ .,data =res.cat1)
res.cox2 <- step(res.cox1)
res.cox3 <- coxph(Surv(OS, STATOS) ~ UICC + X152.996 ,data =res.cat1)
ggforest(model=res.cox3,data=res.cat1)


###等比例风险假定
hypo <- cox.zph(res.cox3) 
hypo1 <- hypo$table[-nrow(hypo$table),]
dim(hypo1)
hypo2 <- hypo1[hypo1[,3]>0.05,]
dim(hypo2)
hypo2

res.cat2 <- res.cat1[,c("OS","STATOS",rownames(hypo2))]
dim(res.cat2)
res.cat2[1:5,1:2]
res.cox4 <- coxph(Surv(OS, STATOS) ~ .,data = res.cat2)
ggforest(model=res.cox4,data=res.cat2)


summ <- as.data.frame(cbind(summary(res.cox4)$coefficients[,c(2,5)],summary(res.cox4)$conf.int[,c(3,4)]))
summ$mz.value <- rownames(summ)
head(summ)
dim(summ)
multivariate.Pancreas.tumor.n107.sum <- summ

summ1 <- summ[-c(1),]
dim(summ1)
head(summ1)

summ2 <- merge(summ1,result8)
multivariate.Pancreas.tumor.n107 <- summ2
dim(multivariate.Pancreas.tumor.n107)
head(multivariate.Pancreas.tumor.n107)
multivariate.Pancreas.tumor.n107.sig <- summ2[summ2[,3]<0.05,]
dim(multivariate.Pancreas.tumor.n107.sig)


#### forest plot
library(forestplot)

result12 <- summ2
dim(result12)
head(result12)

head(summ)
result13 <- summ[1,]
head(result13)

result14 <- data.frame(Feature=c(result13$mz.value,result12$query_mass),
                       Annotation=c("1/2/3/4",result12$name),
                       HR=c(result13[,1],result12[,2]),
                       HR.lower=c(result13[,3],result12[,4]),
                       HR.upper=c(result13[,4],result12[,5]),
                       pvalue=round(c(result13[,2],result12[,3]),3))
dim(result14)
head(result14)

result14$HR1 <- paste(round(result14$HR,2),"(",round(result14$HR.lower,2),"-",round(result14$HR.upper,2),")",sep = '')
result14$pvalue1 <- ifelse(result14$pvalue<0.001,"***",ifelse(result14$pvalue<0.01,"**",ifelse(result14$pvalue<0.05,"*"," ")))
result14$pvalue2 <- paste(result14$pvalue,result14$pvalue1,sep = '')
dim(result14)
head(result14)

labeltext <- as.matrix(rbind(c("Feature","Annotation","HR","pvalue"),result14[,c(1,2,7,9)]))
labeltext

forestplot(labeltext, 
           mean = c(NA,result14$HR),
           lower = c(NA,result14$HR.lower),
           upper = c(NA,result14$HR.upper),
           zero = 1,lwd.zero = 2,
           clip = c(0,10),
           graph.pos = 3,
           boxsize = 0.5,
           align = "l",
           xlog = F,
           col=fpColors(line = "black", zero = "black",
                        box="black"),
           lty.ci = 7,
           lwd.ci = 1,
           ci.vertices.height = 0.15,
           lwd.xaxis=2,
           xlab="Hazard Ratio",
           txt_gp = fpTxtGp(ticks = gpar(cex = 0.9), xlab = gpar(cex = 0.8), cex = 0.8),
           lineheight = unit(8,"mm")) 



############### use Cox model to predict survival outcome by ROC curve
res.cat3 <- res.cat1[,c("OS","STATOS",rownames(summ))]
dim(res.cat3)
res.cat3[1:5,1:4]

res.cox3 <- coxph(Surv(OS,STATOS)~.,data =res.cat3)
sum2 <- summary(res.cox3)
# C index export
sum2$concordance

# forest plot
ggforest(model = res.cox3,data=res.cat3)

library(pROC)
gfit <- roc(STATOS ~predict(res.cox3), data = res.cat3)

plot(gfit,
     print.auc=TRUE, #输出AUC值
     print.thres=TRUE, #输出cut-off值
     #main = "ROC Curve in Lung NAC tumor (n=77)", #设置图形的标题
     col= "red", #曲线颜色
     print.thres.col="black", #cut-off值字体的颜色
     identity.col="blue", #对角线颜色
     identity.lty=1,identity.lwd=1)


################# use Cox model to predict good/poor prognosis or low/high risk groups by ROC curve
res.cat3$STATOS1 <- ifelse(res.cat3$OS < median(res.cat3$OS),1,0) #1=short/poor survival or high risk,0=long/good survival or low risk  
res.cat3$STATOS1

gfit1 <- roc(STATOS1 ~predict(res.cox3), data = res.cat3)
cutoff <- gfit1$thresholds[which.max(gfit1$sensitivities+gfit1$specificities)]
cutoff
plot(gfit1,
     print.auc=TRUE, #输出AUC值
     print.thres=TRUE, #输出cut-off值
     #main = "ROC Curve in Lung NAC tumor (n=77)", #设置图形的标题
     col= "red", #曲线颜色
     print.thres.col="black", #cut-off值字体的颜色
     identity.col="blue", #对角线颜色
     identity.lty=1,identity.lwd=1)

### cutoff-based cox model to stratify patients into different prognostic risk groups by survival curve
res.cat3[1:5,1:ncol(res.cat3)]
res.cat3$risk <- ifelse(predict(res.cox3) < cutoff,'low','high')
# fit survival curve
fit <- survfit(Surv(OS,STATOS)~risk,data =res.cat3)
fit
##绘制生存曲线##
ggsurvplot(fit, data = res.cat3, conf.int = F, # 增加置信区间
           risk.table = TRUE, pval=T) #不同时间点风险人数表，增加P值



################# Nomogram绘制
library(car)
library(rms)
library(pROC)
library(timeROC)
library(ggDCA)

res.cat4 <- res.cat1[,c("OS","STATOS",rownames(summ))]
dim(res.cat4)
res.cat4[1:5,]

res.cox4 <- cph(as.formula(paste0('Surv(OS,STATOS)~',paste(colnames(res.cat4)[-c(1:2)],sep = '',collapse = '+'))),
                x=T,y=T,surv=T,data =res.cat4)
print(res.cox4)

ddist <- datadist(res.cat4)
options(datadist='ddist')

surv <- Survival(res.cox4) # 建立生存函数
surv

surv1 <- function(x) surv(12,x) # 1年OS
surv2 <- function(x) surv(36,x) # 3年OS
surv3 <- function(x) surv(60,x) # 5年OS

nom <- nomogram(res.cox4,
                fun = list(surv1,surv2,surv3),
                funlabel = c('1-year survival Probability',
                             '3-year survival Probability',
                             '5-year survival Probability'))
plot(nom, 
     lplabel="Linear Predictor",
     xfrac = 0.2, # 左侧标签距离坐标轴的距离
     #varname.label = TRUE, 
     tcl = -0.2, # 刻度长短和方向 
     lmgp = 0.1, # 坐标轴标签距离坐标轴远近
     points.label ='Points', 
     total.points.label = 'Total Points',
     cap.labels = FALSE,
     cex.var = 0.7, # 左侧标签字体大小
     cex.axis = 0.6, # 坐标轴字体大小
     col.grid = gray(c(0.8, 0.95))) # 竖线颜色


############# ROC curve on the prediction of 1/3/5-year survival rate
res.cat4 <- res.cat1[,c("OS","STATOS",rownames(summ))]
dim(res.cat4)
res.cat4[1:5,1:3]

res.cox4 <- cph(as.formula(paste0('Surv(OS,STATOS)~',paste(colnames(res.cat4)[-c(1:2)],sep = '',collapse = '+'))),
                x=T,y=T,surv=T,data =res.cat4)
print(res.cox4)

pred_f_training <- predict(res.cox4,res.cat4,type="lp")#!!!type="lp",是他没错
data_table <- data.frame(time=res.cat4[,"OS"],status=res.cat4[,"STATOS"],score=pred_f_training)


time_roc_res <- timeROC(
  T = data_table$time,
  delta = data_table$status,
  marker = data_table$score,
  cause = 1,
  weighting="marginal",
  times = c(12, 36, 60),
  ROC = TRUE,
  iid = TRUE
)


time_ROC_df <- data.frame(
  TP_1year = time_roc_res$TP[, 1],
  FP_1year = time_roc_res$FP[, 1],
  TP_3year = time_roc_res$TP[, 2],
  FP_3year = time_roc_res$FP[, 2],
  TP_5year = time_roc_res$TP[, 3],
  FP_5year = time_roc_res$FP[, 3]
)


ggplot(data = time_ROC_df) +
  geom_line(aes(x = FP_1year, y = TP_1year), linewidth = 1, color = "#BC3C29FF") +
  geom_line(aes(x = FP_3year, y = TP_3year), linewidth = 1, color = "#0072B5FF") +
  geom_line(aes(x = FP_5year, y = TP_5year), linewidth = 1, color = "#E18727FF") +
  geom_abline(slope = 1, intercept = 0, color = "grey", linewidth = 1, linetype = 2) +
  theme_bw() +
  annotate("text",
           x = 0.75, y = 0.25, size = 4.5,
           label = paste0("AUC at 1 year = ", sprintf("%.3f", time_roc_res$AUC[[1]])), color = "#BC3C29FF"
  ) +
  annotate("text",
           x = 0.75, y = 0.15, size = 4.5,
           label = paste0("AUC at 3 years = ", sprintf("%.3f", time_roc_res$AUC[[2]])), color = "#0072B5FF"
  ) +
  annotate("text",
           x = 0.75, y = 0.05, size = 4.5,
           label = paste0("AUC at 5 years = ", sprintf("%.3f", time_roc_res$AUC[[3]])), color = "#E18727FF"
  ) +
  labs(x = "Specificity", y = "Sensitivity") +
  theme(
    axis.text = element_text(size = 10, color = "black"),
    axis.title.x = element_text(size = 12, color = "black", margin = margin(c(15, 0, 0, 0))),
    axis.title.y = element_text(size = 12, color = "black", margin = margin(c(0, 15, 0, 0)))
  )


prognosis.n107 <- data.frame(Tumor=rep("Pancreatic cancer",5),
                             Prognosis=c("Survival outcome","Risk stratification","1-year survival rate","3-year survival rate","5-year survival rate"),
                             AUC=c(gfit$auc,gfit1$auc,time_roc_res$AUC[[1]],time_roc_res$AUC[[2]],time_roc_res$AUC[[3]]))
prognosis.n107


############################## Gastric Cancer (n=246) #####################
library(survival)
library(survminer)

setwd("C:/MyRdata8/Gastric Cancer (n=246)")
table <- read.csv("Gastric Cancer (n=246).csv")
dim(table)
table[1:5,1:15]

table1 <- table[,c(1:5,10,11,13:9290)]
dim(table1)
table1[1:5,1:10]
rownames(table1) <- table1[,1]
table2 <- table1[,-1]
dim(table2)
table2[1:5,1:10]

table3 <- table2[,-c(1:6)]
dim(table3)
table3[1:5,1:7]
table3.n246 <- table3

# nucleotide profile data
table3 <- table3[,match(intersect(round(as.numeric(substr(colnames(table3),2,12)),4),nucleotide.annotation$query_mass),round(as.numeric(substr(colnames(table3),2,12)),4))]
dim(table3)

# median intensity on each mass
median1 <- apply(table3,MARGIN=2,FUN=median)
head(median1)
length(median1)
median1.n246 <- median1

# save nucleotide profile data
nucleotide.profile.Gastric.Cancer.n246 <- table3
head(nucleotide.profile.Gastric.Cancer.n246)


# calculate zero number for each m/z value
len1 <- c()
for (i in 1:ncol(table3)) {
  len <- length(which(table3[,i]==0))
  len1 <- c(len1,len)
}
len1
table4 <- table3[,which(len1 <= 246)]
dim(table4)

table5 <- cbind(table2[,3:4],table4)
dim(table5)
table5[1:5,1:7]

res.cut <- surv_cutpoint(table5, time = "OS", event = "STATOS",
                         minprop = 0.1,
                         variables = colnames(table5)[3:ncol(table5)])


res.cat.1 <- surv_categorize(res.cut,labels = c(0, 1))
dim(res.cat.1)
res.cat.1[1:5,1:7]

table2[1:5,1:10]
res.cat <- cbind(table2[,c(1,2,6)],res.cat.1)
dim(res.cat)
res.cat[1:5,1:10]

colnames(res.cat)[c(1:3)] <- c("Age","Gender","UICC")
res.cat[1:5,1:10]
res.cat$Age <- ifelse(res.cat$Age<60,0,1)
res.cat[1:5,1:10]

# save profile data with cutoff
res.cat.Gastric.Cancer.n246 <- res.cat
dim(res.cat.Gastric.Cancer.n246)
head(res.cat.Gastric.Cancer.n246)



################ univariate cox analysis to calculate p value, HR
result <- data.frame()
for (i in c(1:3,6:ncol(res.cat))) {
  res.cox <- coxph(Surv(OS,STATOS)~res.cat[,i],data =res.cat)
  summ <- summary(res.cox)
  hypo <- cox.zph(res.cox)
  hypo1 <- hypo$table[-nrow(hypo$table),]
  result1 <- data.frame(mz.value=colnames(res.cat)[i],
                        HR = signif(summ$coef[2], digits=3),
                        HR.confint.lower = signif(summ$conf.int[,"lower .95"], 3),
                        HR.confint.upper = signif(summ$conf.int[,"upper .95"],3),
                        logtest.p=summ$logtest[3],
                        waldtest.p=summ$waldtest[3],
                        sctest.p=summ$sctest[3],
                        hypo.p=hypo1[3])
  result <- rbind(result,result1)
}
head(result)
dim(result)
univariate.Gastric.Cancer.n246 <- result


# screen by sctest p value
result2 <- result[result$sctest.p<0.05&result$hypo.p>0.05,]
dim(result2)
head(result2)
result3 <- result2[-c(1:2),]
head(result3)
dim(result3)

# annotation on significant mass
result3$query_mass <- round(as.numeric(substr(result3$mz.value,2,12)),4)
result4 <- merge(result3,nucleotide.annotation)
result4$name.adduct <- paste(result4$adduct,result4$name,sep = '.')
result4$median.intensity <- median1[match(result4$mz.value,colnames(table3))]
dim(result4)
result5 <- result4[order(result4$ppm,decreasing = F),]
result6 <- result5[which(duplicated(result5$name.adduct)==F),]
dim(result6)
result7 <- result6[order(result6$median.intensity,decreasing = T),]
result8 <- result7[which(duplicated(result7$name)==F),]
result8 <- result8[order(result8$query_mass,decreasing = F),]
dim(result8)
head(result8)
univariate4.Gastric.Cancer.n246 <- result8


############## validate by survival curve
# fit survival curve
fit <- survfit(Surv(OS,STATOS)~Age,data =res.cat)
fit
survdiff(Surv(OS,STATOS)~Age,data =res.cat) #采用log-rank 检验分析生存率差异

##绘制生存曲线##
ggsurvplot(fit, data = res.cat, conf.int = F, # 增加置信区间
           risk.table = TRUE, pval=T) #不同时间点风险人数表，增加P值



################# forest plot
library(forestplot)

dim(result8)
head(result8)
length(unique(result8$query_mass))
head(result2)

result10 <- data.frame(Feature=c(result2$mz.value[1:2],result8$query_mass),
                       Annotation=c("60","0/1/2/3/4",result8$name),
                       HR=c(result2$HR[1:2],result8$HR),
                       HR.lower=c(result2$HR.confint.lower[1:2],result8$HR.confint.lower),
                       HR.upper=c(result2$HR.confint.upper[1:2],result8$HR.confint.upper),
                       pvalue=round(c(result2$sctest.p[1:2],result8$sctest.p),3))
dim(result10)
head(result10)

result10$HR1 <- paste(result10$HR,"(",result10$HR.lower,"-",result10$HR.upper,")",sep = '')
dim(result10)
head(result10)

labeltext <- as.matrix(rbind(c("Feature","Annotation","HR","pvalue"),result10[,c(1,2,7,6)]))
labeltext

forestplot(labeltext, 
           mean = c(NA,result10$HR),
           lower = c(NA,result10$HR.lower),
           upper = c(NA,result10$HR.upper),
           zero = 1,lwd.zero = 2,
           clip = c(0,6),
           graph.pos = 3,
           boxsize = 0.5,
           align = "l",
           xlog = FALSE,
           col=fpColors(line = "black", zero = "black",
                        box="black"),
           lty.ci = 7,
           lwd.ci = 1,
           ci.vertices.height = 0.15,
           lwd.xaxis=2,
           xlab="Hazard Ratio",
           txt_gp = fpTxtGp(ticks = gpar(cex = 0.9), xlab = gpar(cex = 0.8), cex = 0.7),
           lineheight = unit(8,"mm")) 



################################################ multivariate cox analysis identify independent prognostic factor
dim(res.cat)
res.cat[1:5,1:10]

head(result2)
head(result8)

res.cat1 <- res.cat[,c("OS","STATOS",result2$mz.value[c(1,2)],unique(result8$mz.value))]
dim(res.cat1)
res.cat1[1:5,1:7]

# save final profiling data after annotation
res.cat1.Gastric.Cancer.n246 <- res.cat1


################## method1 feature selection from the beginning

### Step regression on Cox model
res.cox1 <- coxph(Surv(OS, STATOS) ~ .,data =res.cat1)
res.cox2 <- step(res.cox1)
res.cox3 <- coxph(Surv(OS, STATOS) ~ Age + UICC + X152.9958 + X171.0065 + X192.964 + 
                    X233.0667 + X265.0257 + X272.0775 + X279.0387 + X304.034 + 
                    X313.033 + X322.0446 + X328.0346 + X343.0342 + X375.0119 + 
                    X401.0272 + X424.9184,data =res.cat1)
ggforest(model=res.cox3,data=res.cat1)


###等比例风险假定
hypo <- cox.zph(res.cox3)
#hypo <- cox.zph(res.cox4) 
hypo1 <- hypo$table[-nrow(hypo$table),]
dim(hypo1)
hypo2 <- hypo1[hypo1[,3]>0.05,]
dim(hypo2)
hypo2

res.cat2 <- res.cat1[,c("OS","STATOS",rownames(hypo2))]
dim(res.cat2)
res.cat2[1:5,1:5]
res.cox4 <- coxph(Surv(OS, STATOS) ~ .,data = res.cat2)
ggforest(model=res.cox4,data=res.cat2)


summ <- as.data.frame(cbind(summary(res.cox4)$coefficients[,c(2,5)],summary(res.cox4)$conf.int[,c(3,4)]))
summ$mz.value <- rownames(summ)
head(summ)
dim(summ)
multivariate.Gastric.Cancer.n246.sum <- summ

summ1 <- summ[-c(1),]
dim(summ1)
head(summ1)

summ2 <- merge(summ1,result8)
multivariate.Gastric.Cancer.n246 <- summ2
dim(multivariate.Gastric.Cancer.n246)
head(multivariate.Gastric.Cancer.n246)
multivariate.Gastric.Cancer.n246.sig <- summ2[summ2[,3]<0.05,]
dim(multivariate.Gastric.Cancer.n246.sig)


#### forest plot
library(forestplot)

result12 <- summ2
dim(result12)
head(result12)

head(summ)
result13 <- summ[1,]
head(result13)

result14 <- data.frame(Feature=c(result13$mz.value,result12$query_mass),
                       Annotation=c("60",result12$name),
                       HR=c(result13[,1],result12[,2]),
                       HR.lower=c(result13[,3],result12[,4]),
                       HR.upper=c(result13[,4],result12[,5]),
                       pvalue=round(c(result13[,2],result12[,3]),3))
dim(result14)
head(result14)

result14$HR1 <- paste(round(result14$HR,2),"(",round(result14$HR.lower,2),"-",round(result14$HR.upper,2),")",sep = '')
result14$pvalue1 <- ifelse(result14$pvalue<0.001,"***",ifelse(result14$pvalue<0.01,"**",ifelse(result14$pvalue<0.05,"*"," ")))
result14$pvalue2 <- paste(result14$pvalue,result14$pvalue1,sep = '')
dim(result14)
head(result14)

labeltext <- as.matrix(rbind(c("Feature","Annotation","HR","pvalue"),result14[,c(1,2,7,9)]))
labeltext

forestplot(labeltext, 
           mean = c(NA,result14$HR),
           lower = c(NA,result14$HR.lower),
           upper = c(NA,result14$HR.upper),
           zero = 1,lwd.zero = 2,
           clip = c(0.01,40),
           graph.pos = 3,
           boxsize = 0.5,
           align = "l",
           xlog = T,
           col=fpColors(line = "black", zero = "black",
                        box="black"),
           lty.ci = 7,
           lwd.ci = 1,
           ci.vertices.height = 0.15,
           lwd.xaxis=2,
           xlab="Hazard Ratio",
           txt_gp = fpTxtGp(ticks = gpar(cex = 0.9), xlab = gpar(cex = 0.8), cex = 0.8),
           lineheight = unit(8,"mm")) 




############### use Cox model to predict survival outcome by ROC curve
res.cat3 <- res.cat1[,c("OS","STATOS",rownames(summ))]
dim(res.cat3)
res.cat3[1:5,1:5]

res.cox3 <- coxph(Surv(OS,STATOS)~.,data =res.cat3)
sum2 <- summary(res.cox3)
# C index export
sum2$concordance

# forest plot
ggforest(model = res.cox3,data=res.cat3)

library(pROC)
gfit <- roc(STATOS ~predict(res.cox3), data = res.cat3)

plot(gfit,
     print.auc=TRUE, #输出AUC值
     print.thres=TRUE, #输出cut-off值
     #main = "ROC Curve in Lung NAC tumor (n=77)", #设置图形的标题
     col= "red", #曲线颜色
     print.thres.col="black", #cut-off值字体的颜色
     identity.col="blue", #对角线颜色
     identity.lty=1,identity.lwd=1)


################# use Cox model to predict good/poor prognosis or low/high risk groups by ROC curve
res.cat3$STATOS1 <- ifelse(res.cat3$OS < median(res.cat3$OS),1,0) #1=short/poor survival or high risk,0=long/good survival or low risk  
res.cat3$STATOS1

gfit1 <- roc(STATOS1 ~predict(res.cox3), data = res.cat3)
cutoff <- gfit1$thresholds[which.max(gfit1$sensitivities+gfit1$specificities)]
cutoff
plot(gfit1,
     print.auc=TRUE, #输出AUC值
     print.thres=TRUE, #输出cut-off值
     #main = "ROC Curve in Lung NAC tumor (n=77)", #设置图形的标题
     col= "red", #曲线颜色
     print.thres.col="black", #cut-off值字体的颜色
     identity.col="blue", #对角线颜色
     identity.lty=1,identity.lwd=1)

### cutoff-based cox model to stratify patients into different prognostic risk groups by survival curve
res.cat3[1:5,1:ncol(res.cat3)]
res.cat3$risk <- ifelse(predict(res.cox3) < cutoff,'low','high')
# fit survival curve
fit <- survfit(Surv(OS,STATOS)~risk,data =res.cat3)
fit
##绘制生存曲线##
ggsurvplot(fit, data = res.cat3, conf.int = F, # 增加置信区间
           risk.table = TRUE, pval=T) #不同时间点风险人数表，增加P值



################# Nomogram绘制

library(car)
library(rms)
library(pROC)
library(timeROC)
library(ggDCA)

res.cat4 <- res.cat1[,c("OS","STATOS",rownames(summ))]
dim(res.cat4)
res.cat4[1:5,]

res.cox4 <- cph(as.formula(paste0('Surv(OS,STATOS)~',paste(colnames(res.cat4)[-c(1:2)],sep = '',collapse = '+'))),
                x=T,y=T,surv=T,data =res.cat4)
print(res.cox4)

ddist <- datadist(res.cat4)
options(datadist='ddist')

surv <- Survival(res.cox4) # 建立生存函数
surv

surv1 <- function(x) surv(12,x) # 1年OS
surv2 <- function(x) surv(36,x) # 3年OS
surv3 <- function(x) surv(60,x) # 5年OS

nom <- nomogram(res.cox4,
                fun = list(surv1,surv2,surv3),
                funlabel = c('1-year survival Probability',
                             '3-year survival Probability',
                             '5-year survival Probability'))
plot(nom, 
     lplabel="Linear Predictor",
     xfrac = 0.2, # 左侧标签距离坐标轴的距离
     #varname.label = TRUE, 
     tcl = -0.2, # 刻度长短和方向 
     lmgp = 0.1, # 坐标轴标签距离坐标轴远近
     points.label ='Points', 
     total.points.label = 'Total Points',
     cap.labels = FALSE,
     cex.var = 0.7, # 左侧标签字体大小
     cex.axis = 0.6, # 坐标轴字体大小
     col.grid = gray(c(0.8, 0.95))) # 竖线颜色


############# ROC curve on the prediction of 1/3/5-year survival rate
res.cat4 <- res.cat1[,c("OS","STATOS",rownames(summ))]
dim(res.cat4)
res.cat4[1:5,1:5]

res.cox4 <- cph(as.formula(paste0('Surv(OS,STATOS)~',paste(colnames(res.cat4)[-c(1:2)],sep = '',collapse = '+'))),
                x=T,y=T,surv=T,data =res.cat4)
print(res.cox4)

pred_f_training <- predict(res.cox4,res.cat4,type="lp")#!!!type="lp",是他没错
data_table <- data.frame(time=res.cat4[,"OS"],status=res.cat4[,"STATOS"],score=pred_f_training)


time_roc_res <- timeROC(
  T = data_table$time,
  delta = data_table$status,
  marker = data_table$score,
  cause = 1,
  weighting="marginal",
  times = c(12, 36, 60),
  ROC = TRUE,
  iid = TRUE
)


time_ROC_df <- data.frame(
  TP_1year = time_roc_res$TP[, 1],
  FP_1year = time_roc_res$FP[, 1],
  TP_3year = time_roc_res$TP[, 2],
  FP_3year = time_roc_res$FP[, 2],
  TP_5year = time_roc_res$TP[, 3],
  FP_5year = time_roc_res$FP[, 3]
)


ggplot(data = time_ROC_df) +
  geom_line(aes(x = FP_1year, y = TP_1year), linewidth = 1, color = "#BC3C29FF") +
  geom_line(aes(x = FP_3year, y = TP_3year), linewidth = 1, color = "#0072B5FF") +
  geom_line(aes(x = FP_5year, y = TP_5year), linewidth = 1, color = "#E18727FF") +
  geom_abline(slope = 1, intercept = 0, color = "grey", linewidth = 1, linetype = 2) +
  theme_bw() +
  annotate("text",
           x = 0.75, y = 0.25, size = 4.5,
           label = paste0("AUC at 1 year = ", sprintf("%.3f", time_roc_res$AUC[[1]])), color = "#BC3C29FF"
  ) +
  annotate("text",
           x = 0.75, y = 0.15, size = 4.5,
           label = paste0("AUC at 3 years = ", sprintf("%.3f", time_roc_res$AUC[[2]])), color = "#0072B5FF"
  ) +
  annotate("text",
           x = 0.75, y = 0.05, size = 4.5,
           label = paste0("AUC at 5 years = ", sprintf("%.3f", time_roc_res$AUC[[3]])), color = "#E18727FF"
  ) +
  labs(x = "Specificity", y = "Sensitivity") +
  theme(
    axis.text = element_text(size = 10, color = "black"),
    axis.title.x = element_text(size = 12, color = "black", margin = margin(c(15, 0, 0, 0))),
    axis.title.y = element_text(size = 12, color = "black", margin = margin(c(0, 15, 0, 0)))
  )


prognosis.n246 <- data.frame(Tumor=rep("Gastric cancer",5),
                             Prognosis=c("Survival outcome","Risk stratification","1-year survival rate","3-year survival rate","5-year survival rate"),
                             AUC=c(gfit$auc,gfit1$auc,time_roc_res$AUC[[1]],time_roc_res$AUC[[2]],time_roc_res$AUC[[3]]))
prognosis.n246



############################## Advanced gastric cancer (n=82) #####################
library(survival)
library(survminer)

setwd("C:/MyRdata8/HER2 Test Deviations (n=82)")
table <- read.csv("HER2 Test Deviations (n=82).csv")
dim(table)
table[1:5,1:20]

table1 <- table
dim(table1)
table1[1:5,1:7]
rownames(table1) <- table1[,1]
table2 <- table1[,-1]
dim(table2)
table2[1:5,1:15]
surv.n82 <- table2
surv.n82[1:5,1:15]

table3 <- table2[,-c(1:12)]
dim(table3)
table3[1:5,1:7]
table3.n82 <- table3
table3.n82.1 <- table3[which(table2$classifier=="resistant"|table2$classifier=="sensitive"),]

# nucleotide profile data
table3 <- table3[,match(intersect(round(as.numeric(substr(colnames(table3),2,12)),4),nucleotide.annotation$query_mass),
                        round(as.numeric(substr(colnames(table3),2,12)),4))]
dim(table3)


# median intensity on each mass
median1 <- apply(table3,MARGIN=2,FUN=median)
head(median1)
length(median1)
median1.n82 <- median1

# save nucleotide profile data
nucleotide.profile.HER2.n82 <- table3
head(nucleotide.profile.HER2.n82)


# calculate zero number for each m/z value
len1 <- c()
for (i in 1:ncol(table3)) {
  len <- length(which(table3[,i]==0))
  len1 <- c(len1,len)
}
len1
table4 <- table3[,which(len1 <= 74)]
dim(table4)

table5 <- cbind(table2[,1:2],table4)
dim(table5)
table5[1:5,1:7]

res.cut <- surv_cutpoint(table5, time = "OS", event = "STATOS",
                         minprop = 0.1,
                         variables = colnames(table5)[3:ncol(table5)])
res.cat.1 <- surv_categorize(res.cut,labels = c(0, 1))
dim(res.cat.1)
res.cat.1[1:5,1:7]

res.cat <- cbind(table2[,c(3:12)],res.cat.1)
dim(res.cat)
res.cat[1:5,1:15]

res.cat$age_years_calc <- ifelse(res.cat$age_years_calc<70,1,2)
res.cat$Visite1.Geschlecht <- ifelse(res.cat$Visite1.Geschlecht=="m",1,2)
res.cat$Tumorprobe.SP.state <- ifelse(res.cat$Tumorprobe.SP.state=="positiv",1,0)
#res.cat$calc.Trastuzumab.Visite.1.5.Trastuzumab <- ifelse((res.cat$classifier=="sensitive")|(res.cat$classifier=="resistant"),1,0)
res.cat$calc.Trastuzumab.Visite.1.5.Trastuzumab <- ifelse(res.cat$calc.Trastuzumab.Visite.1.5.Trastuzumab=="Ja",1,0)
res.cat$Tumorprobe.Typ <- ifelse(res.cat$Tumorprobe.Typ=="Prim\xe4rbiopsie",1,ifelse(res.cat$Tumorprobe.Typ=="Metastase",3,2))
res.cat[1:5,1:15]
dim(res.cat)

# save profile data with cutoff
res.cat.HER2.n82 <- res.cat
dim(res.cat.HER2.n82)
head(res.cat.HER2.n82)


################ calculate p value, HR
result <- data.frame()
for (i in c(1:8,10,13:ncol(res.cat))) {
  res.cox <- coxph(Surv(OS,STATOS)~res.cat[,i],data =res.cat)
  summ <- summary(res.cox)
  hypo <- cox.zph(res.cox)
  hypo1 <- hypo$table[-nrow(hypo$table),]
  result1 <- data.frame(mz.value=colnames(res.cat)[i],
                        #cutoff = res.cut[[i-2]][["estimate"]][["estimated cutpoint"]],
                        HR = signif(summ$coef[2], digits=3),
                        HR.confint.lower = signif(summ$conf.int[,"lower .95"], 3),
                        HR.confint.upper = signif(summ$conf.int[,"upper .95"],3),
                        logtest.p=summ$logtest[3],
                        waldtest.p=summ$waldtest[3],
                        sctest.p=summ$sctest[3],
                        hypo.p=hypo1[3])
  result <- rbind(result,result1)
}
head(result)
dim(result)
univariate.HER2.n82 <- result

# screen by sctest p value
result2 <- result[result$sctest.p<0.05&result$hypo.p>0.05,]
dim(result2)
head(result2)
result3 <- result2[-c(1),]
head(result3)
dim(result3)

# annotation on significant mass
result3$query_mass <- round(as.numeric(substr(result3$mz.value,2,12)),4)
result4 <- merge(result3,nucleotide.annotation)
result4$name.adduct <- paste(result4$adduct,result4$name,sep = '.')
result4$median.intensity <- median1[match(result4$mz.value,colnames(table3))]
dim(result4)
result5 <- result4[order(result4$ppm,decreasing = F),]
result6 <- result5[which(duplicated(result5$name.adduct)==F),]
dim(result6)
result7 <- result6[order(result6$median.intensity,decreasing = T),]
result8 <- result7[which(duplicated(result7$name)==F),]
result8 <- result8[order(result8$query_mass,decreasing = F),]
dim(result8)
head(result8)
univariate4.HER2.n82 <- result8


############## validate by survival curve
res.cat[1:5,1:15]
res.cat1 <- res.cat[res.cat$classifier!="",]
# fit survival curve
fit <- survfit(Surv(OS,STATOS)~classifier,data =res.cat)
fit
survdiff(Surv(OS,STATOS)~classifier,data =res.cat) #采用log-rank 检验分析生存率差异

##绘制生存曲线##
ggsurvplot(fit, data = res.cat, conf.int = F, # 增加置信区间
           risk.table = TRUE, pval=T) #不同时间点风险人数表，增加P值



################# forest plot
library(forestplot)

dim(result8)
head(result8)
length(unique(result8$query_mass))
head(result2)
#result2$mz.value[1] <- 'Treatment'

result10 <- data.frame(Feature=c(result2$mz.value[1],result8$query_mass),
                       Annotation=c("0/1",result8$name),
                       HR=c(result2$HR[1],result8$HR),
                       HR.lower=c(result2$HR.confint.lower[1],result8$HR.confint.lower),
                       HR.upper=c(result2$HR.confint.upper[1],result8$HR.confint.upper),
                       pvalue=round(c(result2$sctest.p[1],result8$sctest.p),3))
dim(result10)
head(result10)

result10$HR1 <- paste(result10$HR,"(",result10$HR.lower,"-",result10$HR.upper,")",sep = '')
dim(result10)
head(result10)

labeltext <- as.matrix(rbind(c("Feature","Annotation","HR","pvalue"),result10[,c(1,2,7,6)]))
labeltext

forestplot(labeltext, 
           mean = c(NA,result10$HR),
           lower = c(NA,result10$HR.lower),
           upper = c(NA,result10$HR.upper),
           zero = 1,lwd.zero = 2,
           clip = c(0,6),
           graph.pos = 3,
           boxsize = 0.5,
           align = "l",
           xlog = FALSE,
           col=fpColors(line = "black", zero = "black",
                        box="black"),
           lty.ci = 7,
           lwd.ci = 1,
           ci.vertices.height = 0.15,
           lwd.xaxis=2,
           xlab="Hazard Ratio",
           txt_gp = fpTxtGp(ticks = gpar(cex = 0.9), xlab = gpar(cex = 0.8), cex = 0.7),
           lineheight = unit(8,"mm")) 



################################### multivariate cox analysis identify independent prognostic factor
dim(res.cat)
res.cat[1:5,1:10]

head(result2)
head(result8)

res.cat1 <- res.cat[,c("OS","STATOS",result2$mz.value[1],unique(result8$mz.value))]
dim(res.cat1)
res.cat1[1:5,1:7]

# save final profiling data after annotation
res.cat1.HER2.n82 <- res.cat1


################## method1 feature selection from the beginning

### Step regression on Cox model
res.cox1 <- coxph(Surv(OS, STATOS) ~ .,data =res.cat1)
res.cox2 <- step(res.cox1)
res.cox3 <- coxph(Surv(OS, STATOS) ~ X150.0419 + X171.0065 + X304.0325 + X321.0499 + 
                    X322.0445 + X408.0117 + X424.0053,data =res.cat1)
ggforest(model=res.cox3,data=res.cat1)


###等比例风险假定
hypo <- cox.zph(res.cox3) 
#hypo <- cox.zph(res.cox4) 
hypo1 <- hypo$table[-nrow(hypo$table),]
dim(hypo1)
hypo2 <- hypo1[hypo1[,3]>0.05,]
dim(hypo2)
hypo2

res.cat2 <- res.cat1[,c("OS","STATOS",rownames(hypo2))]
dim(res.cat2)
res.cat2[1:5,1:7]
res.cox4 <- coxph(Surv(OS, STATOS) ~ .,data = res.cat2)
ggforest(model=res.cox4,data=res.cat2)


summ <- as.data.frame(cbind(summary(res.cox4)$coefficients[,c(2,5)],summary(res.cox4)$conf.int[,c(3,4)]))
summ$mz.value <- rownames(summ)
head(summ)
dim(summ)
multivariate.HER2.n82.sum <- summ

summ1 <- summ
dim(summ1)
head(summ1)

summ2 <- merge(summ1,result8)
multivariate.HER2.n82 <- summ2
dim(multivariate.HER2.n82)
head(multivariate.HER2.n82)
multivariate.HER2.n82.sig <- summ2[summ2[,3]<0.05,]
dim(multivariate.HER2.n82.sig)


#### forest plot
library(forestplot)

result12 <- summ2
dim(result12)
head(result12)


result14 <- data.frame(Feature=c(result12$query_mass),
                       Annotation=c(result12$name),
                       HR=c(result12[,2]),
                       HR.lower=c(result12[,4]),
                       HR.upper=c(result12[,5]),
                       pvalue=round(c(result12[,3]),3))
dim(result14)
head(result14)

result14$HR1 <- paste(round(result14$HR,2),"(",round(result14$HR.lower,2),"-",round(result14$HR.upper,2),")",sep = '')
result14$pvalue1 <- ifelse(result14$pvalue<0.001,"***",ifelse(result14$pvalue<0.01,"**",ifelse(result14$pvalue<0.05,"*"," ")))
result14$pvalue2 <- paste(result14$pvalue,result14$pvalue1,sep = '')
dim(result14)
head(result14)

labeltext <- as.matrix(rbind(c("Feature","Annotation","HR","pvalue"),result14[,c(1,2,7,9)]))
labeltext

forestplot(labeltext, 
           mean = c(NA,result14$HR),
           lower = c(NA,result14$HR.lower),
           upper = c(NA,result14$HR.upper),
           zero = 1,lwd.zero = 2,
           clip = c(0.01,40),
           graph.pos = 3,
           boxsize = 0.5,
           align = "l",
           xlog = T,
           col=fpColors(line = "black", zero = "black",
                        box="black"),
           lty.ci = 7,
           lwd.ci = 1,
           ci.vertices.height = 0.15,
           lwd.xaxis=2,
           xlab="Hazard Ratio",
           txt_gp = fpTxtGp(ticks = gpar(cex = 0.9), xlab = gpar(cex = 0.8), cex = 0.8),
           lineheight = unit(8,"mm")) 



############### use Cox model to predict survival outcome by ROC curve
res.cat3 <- res.cat1[,c("OS","STATOS",rownames(summ))]
dim(res.cat3)
res.cat3[1:5,1:5]

res.cox3 <- coxph(Surv(OS,STATOS)~.,data =res.cat3)
sum2 <- summary(res.cox3)
# C index export
sum2$concordance

# forest plot
ggforest(model = res.cox3,data=res.cat3)

library(pROC)
gfit <- roc(STATOS ~predict(res.cox3), data = res.cat3)

plot(gfit,
     print.auc=TRUE, #输出AUC值
     print.thres=TRUE, #输出cut-off值
     #main = "ROC Curve in Lung NAC tumor (n=77)", #设置图形的标题
     col= "red", #曲线颜色
     print.thres.col="black", #cut-off值字体的颜色
     identity.col="blue", #对角线颜色
     identity.lty=1,identity.lwd=1)


################# use Cox model to predict good/poor prognosis or low/high risk groups or chemotherapy-sensitive/resistant by ROC curve
res.cat3$STATOS1 <- ifelse(res.cat3$OS < median(res.cat3$OS),1,0) #1=short/poor survival or high risk,0=long/good survival or low risk  
res.cat3$STATOS1

gfit1 <- roc(STATOS1 ~predict(res.cox3), data = res.cat3)
cutoff <- gfit1$thresholds[which.max(gfit1$sensitivities+gfit1$specificities)]
cutoff
plot(gfit1,
     print.auc=TRUE, #输出AUC值
     print.thres=TRUE, #输出cut-off值
     #main = "ROC Curve in Lung NAC tumor (n=77)", #设置图形的标题
     col= "red", #曲线颜色
     print.thres.col="black", #cut-off值字体的颜色
     identity.col="blue", #对角线颜色
     identity.lty=1,identity.lwd=1)

### cutoff-based cox model to stratify patients into different prognostic risk groups by survival curve
res.cat3[1:5,1:ncol(res.cat3)]
res.cat3$risk <- ifelse(predict(res.cox3) < cutoff,'low','high')
# fit survival curve
fit <- survfit(Surv(OS,STATOS)~risk,data =res.cat3)
fit
##绘制生存曲线##
ggsurvplot(fit, data = res.cat3, conf.int = F, # 增加置信区间
           risk.table = TRUE, pval=T) #不同时间点风险人数表，增加P值





################# Nomogram绘制
library(car)
library(rms)
library(pROC)
library(timeROC)
library(ggDCA)

res.cat4 <- res.cat1[,c("OS","STATOS",rownames(summ))]
dim(res.cat4)
res.cat4[1:5,]

res.cox4 <- cph(as.formula(paste0('Surv(OS,STATOS)~',paste(colnames(res.cat4)[-c(1:2)],sep = '',collapse = '+'))),
                x=T,y=T,surv=T,data =res.cat4)
print(res.cox4)

ddist <- datadist(res.cat4)
options(datadist='ddist')

surv <- Survival(res.cox4) # 建立生存函数
surv

surv1 <- function(x) surv(12,x) # 1年OS
surv2 <- function(x) surv(36,x) # 3年OS
surv3 <- function(x) surv(60,x) # 5年OS

nom <- nomogram(res.cox4,
                fun = list(surv1,surv2,surv3),
                funlabel = c('1-year survival Probability',
                             '3-year survival Probability',
                             '5-year survival Probability'))
plot(nom, 
     lplabel="Linear Predictor",
     xfrac = 0.2, # 左侧标签距离坐标轴的距离
     #varname.label = TRUE, 
     tcl = -0.2, # 刻度长短和方向 
     lmgp = 0.1, # 坐标轴标签距离坐标轴远近
     points.label ='Points', 
     total.points.label = 'Total Points',
     cap.labels = FALSE,
     cex.var = 0.8, # 左侧标签字体大小
     cex.axis = 0.8, # 坐标轴字体大小
     col.grid = gray(c(0.8, 0.95))) # 竖线颜色


############# ROC curve on the prediction of 1/3/5-year survival rate
res.cat4 <- res.cat1[,c("OS","STATOS",rownames(summ))]
dim(res.cat4)
res.cat4[1:5,1:5]

res.cox4 <- cph(as.formula(paste0('Surv(OS,STATOS)~',paste(colnames(res.cat4)[-c(1:2)],sep = '',collapse = '+'))),
                x=T,y=T,surv=T,data =res.cat4)
print(res.cox4)

pred_f_training <- predict(res.cox4,res.cat4,type="lp")#!!!type="lp",是他没错
data_table <- data.frame(time=res.cat4[,"OS"],status=res.cat4[,"STATOS"],score=pred_f_training)


time_roc_res <- timeROC(
  T = data_table$time,
  delta = data_table$status,
  marker = data_table$score,
  cause = 1,
  weighting="marginal",
  times = c(12, 36, 60),
  ROC = TRUE,
  iid = TRUE
)


time_ROC_df <- data.frame(
  TP_1year = time_roc_res$TP[, 1],
  FP_1year = time_roc_res$FP[, 1],
  TP_3year = time_roc_res$TP[, 2],
  FP_3year = time_roc_res$FP[, 2],
  TP_5year = time_roc_res$TP[, 3],
  FP_5year = time_roc_res$FP[, 3]
)


ggplot(data = time_ROC_df) +
  geom_line(aes(x = FP_1year, y = TP_1year), linewidth = 1, color = "#BC3C29FF") +
  geom_line(aes(x = FP_3year, y = TP_3year), linewidth = 1, color = "#0072B5FF") +
  geom_line(aes(x = FP_5year, y = TP_5year), linewidth = 1, color = "#E18727FF") +
  geom_abline(slope = 1, intercept = 0, color = "grey", linewidth = 1, linetype = 2) +
  theme_bw() +
  annotate("text",
           x = 0.75, y = 0.25, size = 4.5,
           label = paste0("AUC at 1 year = ", sprintf("%.3f", time_roc_res$AUC[[1]])), color = "#BC3C29FF"
  ) +
  annotate("text",
           x = 0.75, y = 0.15, size = 4.5,
           label = paste0("AUC at 3 years = ", sprintf("%.3f", time_roc_res$AUC[[2]])), color = "#0072B5FF"
  ) +
  annotate("text",
           x = 0.75, y = 0.05, size = 4.5,
           label = paste0("AUC at 5 years = ", sprintf("%.3f", time_roc_res$AUC[[3]])), color = "#E18727FF"
  ) +
  labs(x = "Specificity", y = "Sensitivity") +
  theme(
    axis.text = element_text(size = 10, color = "black"),
    axis.title.x = element_text(size = 12, color = "black", margin = margin(c(15, 0, 0, 0))),
    axis.title.y = element_text(size = 12, color = "black", margin = margin(c(0, 15, 0, 0)))
  )



prognosis.n82 <- data.frame(Tumor=rep("HER2-Gastric cancer",5),
                            Prognosis=c("Survival outcome","Risk stratification","1-year survival rate","3-year survival rate","5-year survival rate"),
                            AUC=c(gfit$auc,gfit1$auc,time_roc_res$AUC[[1]],time_roc_res$AUC[[2]],time_roc_res$AUC[[3]]))
prognosis.n82



##################################### RF predict chemotherapy-sensitive/resistant
dim(res.cat.1)
head(res.cat.1)
res.cat.22 <- res.cat.1
res.cat.2 <- res.cat.22[,-1]
res.cat.2$STATOS <- as.factor(res.cat.2$STATOS) #survival outcome
#res.cat.2$STATOS <- as.factor(ifelse(res.cat.1$OS<median(res.cat.1$OS),1,0)) # survival time long or short
dim(res.cat.2)
head(res.cat.2)


# annotation on mass
result3 <- data.frame(query_mass=round(as.numeric(substr(colnames(res.cat.2)[2:59],2,12)),4))
result4 <- merge(result3,nucleotide.annotation)
dim(result4)
head(result4)
result4$name.adduct <- paste(result4$adduct,result4$name,sep = '.')

head(median1)
names(median1)
result4$median.intensity <- median1[match(result4$query_mass,
                                          round(as.numeric(substr(names(median1),2,12)),4))]
dim(result4)
head(result4)
result5 <- result4[order(result4$ppm,decreasing = F),]
result6 <- result5[which(duplicated(result5$name.adduct)==F),]
dim(result6)
result7 <- result6[order(result6$median.intensity,decreasing = T),]
result8 <- result7[which(duplicated(result7$name1)==F),]
result8 <- result8[order(result8$query_mass,decreasing = F),]
dim(result8)
head(result8)
unique(result8$query_mass)

result8$mz.value <- paste("X",result8$query_mass,sep='')
dim(result8)
head(result8)



# random forest
res.cat.2 <- res.cat.2[,c("STATOS",unique(result8$mz.value))]
dim(res.cat.2)
colnames(res.cat.2)

library(randomForest)
set.seed(12345)
otu_group.forest <- randomForest(STATOS ~ ., data = res.cat.2, importance = TRUE)
otu_group.forest

# 10-fold cross validation
set.seed(12345)
otu_group.cv <- replicate(5, rfcv(res.cat.2[,-1], res.cat.2$STATOS, cv.fold = 10,step = 1.5), simplify = FALSE)
otu_group.cv

otu_group.cv <- data.frame(sapply(otu_group.cv, '[[', 'error.cv'))
otu_group.cv$otus <- rownames(otu_group.cv)
otu_group.cv <- reshape2::melt(otu_group.cv, id = 'otus')
otu_group.cv$otus <- as.numeric(as.character(otu_group.cv$otus))

library(ggplot2)
library(splines)  #用于在 geom_smooth() 中添加拟合线，或者使用 geom_line() 替代 geom_smooth() 绘制普通折线

p <- ggplot(otu_group.cv, aes(otus, value)) +
  geom_smooth(se = FALSE,  method = 'glm', formula = y~ns(x, 6)) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  labs(title = '',x = 'Number of nucleotides', y = 'Cross-validation error')

p
p + geom_vline(xintercept = 26)
importance_otu <- data.frame(importance(otu_group.forest))
importance_otu <- importance_otu[order(importance_otu$MeanDecreaseGini, decreasing = TRUE), ]
head(importance_otu)
dim(importance_otu)

otu_select <- rownames(importance_otu)[1:nrow(importance_otu)]
otu_group2 <- res.cat.2[,c("STATOS",otu_select)]
otu_group2$STATOS <- as.factor(otu_group2$STATOS)
dim(otu_group2)
head(otu_group2)
table(otu_group2$STATOS)


accuracy1 <- data.frame()
for (i in 1:1000) {
  set.seed(i)
  otu_group2.forest <- randomForest(STATOS ~ ., data = otu_group2, importance = TRUE)
  table(otu_group2.forest$predicted,otu_group2$STATOS)
  accuracy <- sum(diag(table(otu_group2.forest$predicted,otu_group2$STATOS)))/sum(table(otu_group2.forest$predicted,otu_group2$STATOS))
  accuracy
  accuracy2 <- data.frame(Accuracy=accuracy,I=i)
  accuracy1 <- rbind(accuracy1,accuracy2)
}
head(accuracy1)
i=1
table(accuracy1$Accuracy)


res.cat.22$newclass <- otu_group2.forest$predicted
dim(res.cat.22)
head(res.cat.22)
# fit survival curve
fit <- survfit(Surv(OS,STATOS)~newclass,data =res.cat.22)
fit
##绘制生存曲线##
ggsurvplot(fit, data = res.cat.22, conf.int = F, # 增加置信区间
           risk.table = TRUE, pval=T) #不同时间点风险人数表，增加P值

res.cox <- coxph(Surv(OS,STATOS)~newclass,data =res.cat.22)
summary(res.cox)
cox.zph(res.cox)



importance_otu <- data.frame(importance(otu_group2.forest))
importance_otu <- importance_otu[order(importance_otu$MeanDecreaseGini, decreasing = TRUE), ]
head(importance_otu)
dim(importance_otu)

rownames(importance_otu) <- result8[match(rownames(importance_otu),result8$mz.value),7]
importance_otu$Nucleotide <- rownames(importance_otu)
head(importance_otu)
dim(importance_otu)

n41 <- result8$name1
n41 <- c(n34[1:3],"dGDP",n34[4:14],"AMP","dGMP",n34[15:18],"Guanine",n34[19:21],"Pseudouridylic acid",
         n34[22:28],"Uridine","Thymidine",n34[29:34])
n41
#write.csv(n40,"41 nucleotides in RF classifier For GC.csv")
setdiff(n41,n34)

n34 <- importance_otu$Nucleotide
n34

importance_otu$Nucleotide <- factor(importance_otu$Nucleotide,levels = importance_otu$Nucleotide)
ggplot(data=importance_otu,aes(x=Nucleotide,y=MeanDecreaseGini))+
  geom_bar(stat="identity")+
  theme(axis.text.x = element_text(size = 10,vjust=0.5,hjust=1,angle=90))




###### compare to other clinical features G staging
clin1 <- read.csv("HER2 gastric cancer n82.csv")
dim(clin1)
head(clin1)
clin1$Gstage1 <- ifelse(clin1$Gstage<3,"1/2","3/4")
table(clin1$Gstage1,clin1$outcome)
sum(diag(table(clin1$Gstage1,clin1$outcome)))/sum(table(clin1$Gstage1,clin1$outcome))
29/60
16/21


# fit survival curve
fit <- survfit(Surv(time,outcome)~Gstage1,data =clin1)
fit
##绘制生存曲线##
ggsurvplot(fit, data = clin1, conf.int = F, # 增加置信区间
           risk.table = TRUE, pval=T) #不同时间点风险人数表，增加P值

res.cox <- coxph(Surv(time,outcome)~Gstage1,data =clin1)
summary(res.cox)
cox.zph(res.cox)




################# correlation of 34 nucleotides abundance and chemotherapy responses
dim(otu_group2)
head(otu_group2)

res <- data.frame()
for (i in 2:35) {
  cor1 <- cor.test(otu_group2[,i],as.numeric(otu_group2$STATOS),method = "kendall")
  res1 <- data.frame(corr=cor1$estimate,pvalue=cor1$p.value)
  res <- rbind(res,res1)
}
dim(res)
head(res)
res$nucleotide <- colnames(otu_group2)[2:35]
res$Nucleotide <- result8[match(res$nucleotide,result8$mz.value),6]
res2 <- res[res$pvalue<0.05,]
res2
write.csv(res2,"8 Significantly correlated nucleotides with chemotherapy in advanced gastric cancer.csv")

library(reshape2)
library(plyr)
otu_group3 <- otu_group2[,c("STATOS",res2$nucleotide)]
colnames(otu_group3)[2:9] <- res2$Nucleotide
head(otu_group3)
dim(otu_group3)
otu_group4 <- table(otu_group3$STATOS,otu_group3[,9])
colnames(otu_group3)
rownames(otu_group4) <- c("Sensitive","Resistant")
colnames(otu_group4) <- c("Low","High")
otu_group4

otu_group5 <- melt(otu_group4)
colnames(otu_group5) <- c("Response","Abundance","Count")
otu_group5

ggplot(data=otu_group5,aes(x=Abundance,y=Count,fill=Response))+
  geom_bar(stat="identity")+
  geom_text(aes(label=Count))+theme(axis.title.x = element_text(size=15),
                                    axis.title.y = element_text(size=15),
                                    axis.text.x = element_text(size=10),
                                    axis.text.y = element_text(size=10))


######################## 41 nucleotides (figure 7) involved pathway/module 
setwd("C:/MyRdata8/Summary")

head(n41)
n40 <- data.frame(name=n41)
dim(n40)
head(n40)

pathway <- read.csv("nucleotide-metabolic module.csv")
pathway <- pathway[,-1]
colnames(pathway)[1] <- "module"
dim(pathway)
head(pathway)
intersect(pathway$name,n40$name)

nucleotide2 <- read.csv("Nucleotide metabolism1.csv")
dim(nucleotide2)
head(nucleotide2)

n40$pathway <- nucleotide2[match(n40$name,nucleotide2$name1),2]
head(n40)
dim(n40)

n40.1 <- data.frame(name=n40$name)
n40.1

head(pathway)
pathway1 <- pathway[,c(1,4)]
head(pathway1)
n40.2 <- merge(n40.1,pathway1)
colnames(n40.2)[2] <- "pathway"
dim(n40.2)

n40.3 <- rbind(n40,n40.2)
n40.3$value <- rep("1",nrow(n40.3))
dim(n40.3)
head(n40.3)

unique(n40.3$pathway)
unique(n40.3$name)


library(ggplot2)
n40.3$name <- factor(n40.3$name,levels = rev(unique(n40.3$name)))
n40.3$pathway <- factor(n40.3$pathway,levels = unique(n40.3$pathway)[c(1,2,4,3,6,5,7,9,12,13,10,11,15,8,14,16)])

ggplot(n40.3,aes(x=pathway,y=name))+
  geom_tile(aes(fill=value),color="black")+
  scale_fill_manual(values = c("orange"))+
  xlab("Nucleotide")+ylab("Cancer")+guides(fill=F)+
  theme(axis.text.x = element_text(size = 7,vjust=1,hjust=1,angle=45),
        axis.text.y=element_text(size = 7))






library(reshape2)
dim(n40.3)
head(n40.3)
n40.3 <- n40.3[,-3]

n40.4 <- dcast(n40.3,pathway~name)
dim(n40.4)
n40.5 <- melt(n40.4)
dim(n40.5)
head(n40.5)
n40.5$value <- ifelse(n40.5$value>0,1,0)

n40.5$variable <- factor(n40.5$variable,levels = unique(n40.3$name))
n40.5$pathway <- factor(n40.5$pathway,levels = rev(unique(n40.3$pathway)))
n40.5$value <- as.factor(n40.5$value)

ggplot(n40.5,aes(x=variable,y=pathway))+
  geom_tile(aes(fill=value),color="black")+
  scale_fill_manual(values = c("white","orange"))+
  xlab("Nucleotide")+ylab("Cancer")+guides(fill=F)+
  theme(axis.text.x = element_text(size = 8,vjust=1,hjust=1,angle=90),
        axis.text.y=element_text(size = 8))


######################## 41 nucleotides (figure 7) involved drug correlation network 
setwd("C:/MyRdata8/Summary")
drug4 <- read.csv("24 nucleotides associated chemotherapy drug pathway1.csv")
dim(drug4)
head(drug4)
drug5 <- drug4[,c(4:6)]
dim(drug5)
head(drug5)

n41.1 <- data.frame(name=n41)
n41.1

n41.2 <- merge(n41.1,drug5,all.x=T)
head(n41.2)
n41.3 <- n41.2
dim(n41.3)
head(n41.3)
unique(n41.3$name)

library(ggplot2)
n41.3$name <- factor(n41.3$name,levels = rev(n41))
head(n41.3)

ggplot(n41.3,aes(x=Drug,y=name))+
  geom_tile(aes(fill=Correlation.type),color="black")+
  scale_fill_manual(values = c("purple","light blue","gold","cyan"),
                    na.translate=T)+
  xlab("Nucleotide")+ylab("Cancer")+guides(fill=F)+
  theme(axis.text.x = element_text(size = 7,vjust=1,hjust=1,angle=45),
        axis.text.y=element_text(size = 7))



############### summary on distribution of 105 detectable nucleotides in cancers #####################

sum <- data.frame(mass=c(colnames(nucleotide.profile.lung.primary.tumor.n85),
                         colnames(nucleotide.profile.lung.NAC.tumor.n77),
                         colnames(nucleotide.profile.lunge.plattenepithel.n238),
                         colnames(nucleotide.profile.all.primary.EAC.n102),
                         colnames(nucleotide.profile.Neoadjuvant.treated.EACs.n144),
                         colnames(nucleotide.profile.Chromophobe.RCC.n100),
                         colnames(nucleotide.profile.Clear.cell.RCC.n476),
                         colnames(nucleotide.profile.Papillary.RCC.n101),
                         colnames(nucleotide.profile.ACC.n72),
                         colnames(nucleotide.profile.Pancreas.tumor.n107),
                         colnames(nucleotide.profile.Gastric.Cancer.n246),
                         colnames(nucleotide.profile.HER2.n82)),
                  intensity=c(median1.n85,median1.n77,median1.n238,median1.n102,
                              median1.n144,median1.n100,median1.n476,median1.n101,
                              median1.n72,median1.n107,median1.n246,median1.n82),
                  tumor=c(rep("Lung Primary Tumor",ncol(nucleotide.profile.lung.primary.tumor.n85)),
                          rep("Lung NAC Tumor",ncol(nucleotide.profile.lung.NAC.tumor.n77)),
                          rep("Lung squamous cell carcinoma",ncol(nucleotide.profile.lunge.plattenepithel.n238)),
                          rep("Primary EAC",ncol(nucleotide.profile.all.primary.EAC.n102)),
                          rep("Neoadjuvant treated EAC",ncol(nucleotide.profile.Neoadjuvant.treated.EACs.n144)),
                          rep("Chromophobe RCC",ncol(nucleotide.profile.Chromophobe.RCC.n100)),
                          rep("Clear cell RCC",ncol(nucleotide.profile.Clear.cell.RCC.n476)),
                          rep("Papillary RCC",ncol(nucleotide.profile.Papillary.RCC.n101)),
                          rep("Adrenocortical Carcinoma",ncol(nucleotide.profile.ACC.n72)),
                          rep("Pancreatic cancer",ncol(nucleotide.profile.Pancreas.tumor.n107)),
                          rep("Gastric cancer",ncol(nucleotide.profile.Gastric.Cancer.n246)),
                          rep("HER2-Gastric cancer",ncol(nucleotide.profile.HER2.n82))))


dim(sum)
head(sum)
sum$query_mass <- round(as.numeric(substr(sum$mass,2,12)),4)
head(sum)

dim(nucleotide.annotation)
head(nucleotide.annotation)

sum1 <- merge(sum,nucleotide.annotation,by="query_mass")
dim(sum1)
head(sum1)

library(dplyr)
library(ggplot2)
setwd("C:/MyRdata8/Summary")
#write.csv(sum1,"annotated and detectable 105 nucleotides in all cancers.csv")
sum2 <- read.csv("annotated and detectable 105 nucleotides in all cancers.csv")
sum3 <- sum2[,-1]
head(sum3)
dim(sum3)

sum4 <- sum3[,c(4,9)]
sum5 <- distinct(sum4)
dim(sum4)
dim(sum5)
head(sum5)
length(unique(sum5$name))


###################### number of measured cancer type for each nucleotide
count1 <- by(sum5,sum5$name,count)
as.numeric(count1)
names(count1)

count2 <- data.frame(name=names(count1),count=as.numeric(count1))
dim(count2)
head(count2)

count3 <- count2[order(count2$count,decreasing = T),]
count3$name <- factor(count3$name,levels = count3$name)

head(count3)
dim(count3)

ggplot(data = count3, mapping = aes(x = name, y = count,group = 1)) + geom_line()+geom_point()+
  theme(axis.text.x = element_blank(),axis.ticks.x=element_blank(),
        axis.title.y = element_text(size=10),
        axis.text.y=element_text(size=10))+
  scale_y_continuous(limits=c(0,12), breaks=seq(0,12,2))+
  xlab('')+ylab('Count')


#################### number of measured nucleotides in each cancer
count4 <- by(sum5,sum5$tumor,count)
as.numeric(count4)
names(count4)

count5 <- data.frame(name=names(count4),count=as.numeric(count4))
dim(count5)
head(count5)
count6 <- count5[order(count5$count,decreasing = F),]
count6$name <- factor(count6$name,levels = count6$name)
dim(count6)
head(count6)

ggplot(data = count6, mapping = aes(x = name, y = count)) + geom_bar(stat = 'identity')+
  theme(axis.text.x = element_text(angle=90,vjust=0.5,hjust=1))+
  scale_y_continuous(limits=c(0,100), breaks=seq(0,100,20))

ggplot(data = count6, mapping = aes(x = name, y = count,group = 1)) + geom_line()+geom_point()+
  theme(axis.text.x = element_blank(),axis.ticks.x=element_blank(),
        axis.title.y = element_text(size=15),
        axis.text.y=element_text(size=15))+
  scale_y_continuous(limits=c(0,100), breaks=seq(0,100,20))+
  xlab('')+ylab('Count')


#################### heatmap shows measurement details for each cancer and nucleotide
dim(sum5)
head(sum5)
count6 <- count5[order(count5$count,decreasing = T),]
count6$name <- factor(count6$name,levels = count6$name)
sum5$value <- rep(1,776)
sum5$value <- as.factor(sum5$value)
sum5$tumor <- factor(sum5$tumor,levels = count6$name)
sum5$name <- factor(sum5$name,levels=count3$name)
head(sum5)
levels(sum5$name)
levels(sum5$tumor)
sum5.5 <- sum5


ggplot(sum5,aes(x=name,y=tumor))+
  geom_tile(aes(fill=value),color="black")+
  scale_fill_manual(values = c("light green"))+
  xlab("Nucleotide")+ylab("Cancer")+guides(fill=F)+
  theme(axis.text.x = element_text(size = 6,vjust=0.5,hjust=1,angle=90),
        axis.text.y=element_text(size = 8))

ggplot(sum5,aes(x=name,y=tumor))+
  geom_tile(aes(fill=value),color="gray",linewidth=0.1)+
  scale_fill_manual(values = c("black"))+
  xlab("Nucleotide")+ylab("Cancer")+guides(fill=F)+
  theme(axis.text.x = element_text(colour="white",size = 14,vjust=0.5,hjust=1,angle=90),
        axis.text.y=element_text(colour="white",size = 16),
        axis.title = element_text(colour = "white",size=20),
        plot.background = element_rect(colour = "black",fill = "black"),
        panel.background = element_rect(colour = "black",fill = "black"),
        panel.border = element_rect(fill=NA,colour = "white",size=1),
        panel.grid = element_line(colour = "gray"))


############### summary KEGG pathways showed 29 common and 105 nucleotides measured in all cancers  ######################
library(pathview)
setwd("C:/MyRdata8/Summary")
sum2 <- read.csv("annotated and detectable 105 nucleotides in all cancers.csv")
sum3 <- sum2[,-1]
head(sum3)
dim(sum3)
length(unique(sum3$name)) #105 nucleotides

sum4.1 <- unique(sum3[which(sum3$tumor==unique(sum3$tumor)[1]),9])
sum4.2 <- unique(sum3[which(sum3$tumor==unique(sum3$tumor)[2]),9])
sum4.3 <- unique(sum3[which(sum3$tumor==unique(sum3$tumor)[3]),9])
sum4.4 <- unique(sum3[which(sum3$tumor==unique(sum3$tumor)[4]),9])
sum4.5 <- unique(sum3[which(sum3$tumor==unique(sum3$tumor)[5]),9])
sum4.6 <- unique(sum3[which(sum3$tumor==unique(sum3$tumor)[6]),9])
sum4.7 <- unique(sum3[which(sum3$tumor==unique(sum3$tumor)[7]),9])
sum4.8 <- unique(sum3[which(sum3$tumor==unique(sum3$tumor)[8]),9])
sum4.9 <- unique(sum3[which(sum3$tumor==unique(sum3$tumor)[9]),9])
sum4.10 <- unique(sum3[which(sum3$tumor==unique(sum3$tumor)[10]),9])
sum4.11 <- unique(sum3[which(sum3$tumor==unique(sum3$tumor)[11]),9])
sum4.12 <- unique(sum3[which(sum3$tumor==unique(sum3$tumor)[12]),9])

sum5 <- intersect(intersect(intersect(intersect(intersect(intersect(intersect(intersect(intersect(intersect(intersect(sum4.1,sum4.2),sum4.3),sum4.4),sum4.5),
                                                                              sum4.6),sum4.7),sum4.8),sum4.9),sum4.10),sum4.11),sum4.12)
sum5

length(sum5)
sum.29.com <- sum5


sum.29 <- sum3[match(sum5,sum3$name),6]
sum.29
sum.105 <- sum3[match(unique(sum3$name),sum3$name),6]
sum.105


head(nucleotide)
dim(nucleotide)

sum.54 <- setdiff(nucleotide$KEGG.ID,sum.105)
sum.54

sum.159 <- c(rep(1,105),rep(0,54))
sum.159
names(sum.159) <- c(sum.105,sum.54)
sum.159

purine.gene <- read.csv("purine gene list.csv")
purine.gene$Name

pyrimidine.gene <- read.csv("pyrimidine gene list.csv")
pyrimidine.gene


sum.29
sum.105
sum.54
sum.159.1 <- sum.159
sum.159.2 <- data.frame(n105=sum.159,n29=sum.159.1)
sum.159.2$n29[match(sum.29,rownames(sum.159.2))] <- rep(-1,length(match(sum.29,rownames(sum.159.2))))
sum.159.2

length(as.character(purine.gene$Name))
purine.gene1 <- data.frame(gene=rep(1,128))
rownames(purine.gene1) <- as.character(purine.gene$Name)

pathview(gene.data=as.character(purine.gene$Name),cpd.data = sum.159.2[,1:2],
         pathway.id = '00230',species = "hsa",out.suffix="purine.gene.cpd2",
         keys.align = "y", kegg.native = T,key.pos = "topright",na.col = "#E5F5E0",
         limit=list(gene = 1,cpd=1),low=list(cpd="red"),mid = list(cpd="gray"),
         high = list(gene="#E5F5E0",cpd="green"),
         discrete = list(gene=T,cpd=T),
         bins = list(gene=1,cpd=2))


pathview(gene.data=as.character(pyrimidine.gene$Name),cpd.data = sum.159.2[,1:2],
         pathway.id = '00240',species = "hsa",out.suffix="pyrimidine.gene.cpd2",
         keys.align = "y", kegg.native = T,key.pos = "topright",na.col = "#E5F5E0",
         limit=list(gene = 1,cpd=1),low=list(cpd="red"),mid = list(cpd="gray"),
         high = list(gene="#E5F5E0",cpd="green"),
         discrete = list(gene=T,cpd=T),
         bins = list(gene=1,cpd=2))


############### summary overall survival curves in all cancer cohorts ###################
clin <- rbind(res.cat.lung.primary.tumor.n85[,c("OS","STATOS")],res.cat.lung.NAC.tumor.n77[,c("OS","STATOS")],
              res.cat.lunge.plattenepithel.n238[,c("OS","STATOS")],res.cat.all.primary.EAC.n102[,c("OS","STATOS")],
              res.cat.Neoadjuvant.treated.EACs.n144[,c("OS","STATOS")],res.cat.Chromophobe.RCC.n100[,c("OS","STATOS")],
              res.cat.Clear.cell.RCC.n476[,c("OS","STATOS")],res.cat.Papillary.RCC.n101[,c("OS","STATOS")],
              res.cat.ACC.n72[,c("OS","STATOS")],res.cat.Pancreas.tumor.n107[,c("OS","STATOS")],
              res.cat.Gastric.Cancer.n246[,c("OS","STATOS")],res.cat.HER2.n82[,c("OS","STATOS")])
dim(clin)
head(clin)
clin$Cancer <- c(rep("Lung Primary Tumor",85),rep("Lung NAC Tumor",77),rep("Lung squamous cell carcinoma",238),
                 rep("Primary EAC",102),rep("Neoadjuvant treated EAC",144),rep("Chromophobe RCC",100),
                 rep("Clear cell RCC",476),rep("Papillary RCC",101),rep("Adrenocortical Carcinoma",72),
                 rep("Pancreatic cancer",107),rep("Gastric cancer",246),rep("HER2-Gastric cancer",82))
dim(clin)
head(clin)

# fit survival curve
fit <- survfit(Surv(OS,STATOS)~Cancer,data =clin)
fit
survdiff(Surv(OS,STATOS)~Cancer,data =clin) #采用log-rank 检验分析生存率差异

##绘制生存曲线##
ggsurvplot(fit, data = clin, conf.int = F, # 增加置信区间
           risk.table = TRUE, pval=T) #不同时间点风险人数表，增加P值

res.cox <- coxph(Surv(OS,STATOS)~Cancer,data =clin)
summary(res.cox)
cox.zph(res.cox)



############################ summary on univariate cox analysis results ########################
sum <- rbind(univariate4.lung.primary.tumor.n85,univariate4.lung.NAC.tumor.n77,univariate4.lunge.plattenepithel.n238,
             univariate4.all.primary.EAC.n102,univariate4.Neoadjuvant.treated.EACs.n144,
             univariate4.Chromophobe.RCC.n100,univariate4.Clear.cell.RCC.n476,
             univariate4.Papillary.RCC.n101,univariate4.ACC.n72,univariate4.Pancreas.tumor.n107,
             univariate4.Gastric.Cancer.n246,univariate4.HER2.n82)
sum$Tumor <- c(rep("Lung Primary Tumor",nrow(univariate4.lung.primary.tumor.n85)),
               rep("Lung NAC Tumor",nrow(univariate4.lung.NAC.tumor.n77)),
               rep("Lung squamous cell carcinoma",nrow(univariate4.lunge.plattenepithel.n238)),
               rep("Primary EAC",nrow(univariate4.all.primary.EAC.n102)),
               rep("Neoadjuvant treated EAC",nrow(univariate4.Neoadjuvant.treated.EACs.n144)),
               rep("Chromophobe RCC",nrow(univariate4.Chromophobe.RCC.n100)),
               rep("Clear cell RCC",nrow(univariate4.Clear.cell.RCC.n476)),
               rep("Papillary RCC",nrow(univariate4.Papillary.RCC.n101)),
               rep("Adrenocortical Carcinoma",nrow(univariate4.ACC.n72)),
               rep("Pancreatic cancer",nrow(univariate4.Pancreas.tumor.n107)),
               rep("Gastric cancer",nrow(univariate4.Gastric.Cancer.n246)),
               rep("HER2-Gastric cancer",nrow(univariate4.HER2.n82)))
dim(sum)
head(sum)
univariate.sum <- sum
setwd("C:/MyRdata8/Summary")
#write.csv(sum,"summary on single metabolite from univariate cox analysis5.csv")
sum <- read.csv("summary on single metabolite from univariate cox analysis5.csv")
sum <- sum[,-1]
dim(sum)
head(sum)

################################# Heatmap diaplay univariate analysis

library(ggplot2)
sum1 <- data.frame(name=c(sum$name),
                   Tumor=c(sum$Tumor),
                   Pathway=c(sum$Pathway),
                   mass=c(sum$query_mass),
                   HR=c(sum$HR),
                   pvalue=c(sum$sctest.p))
dim(sum1)
head(sum1)
sum1$Prognosis <- ifelse(sum1$HR<1,"Good","Poor")
dim(sum1)
head(sum1)
sum2 <- sum1[,c(1,2,7)]
head(sum2)
dim(sum2)
sum2$name1 <- paste(sum2$Tumor,sum2$name,sep = '.')

dim(sum5.5)
sum6 <- sum5.5
head(sum6)
colnames(sum6)[1] <- "Tumor"
head(sum6)
dim(sum6)
sum6$name1 <- paste(sum6$Tumor,sum6$name,sep = '.')
sum6$value <- sum2$Prognosis[match(sum6$name1,sum2$name1)]
sum6$value[which(is.na(sum6$value)==T)] <- "no"

ggplot(sum6,aes(x=name,y=Tumor))+
  geom_tile(aes(fill=value),color="black")+
  scale_fill_manual(values = c("green","white","red"))+
  xlab("Nucleotide")+ylab("Cancer")+guides(fill=F)+
  theme(axis.text.x = element_text(size = 6,vjust=0.5,hjust=1,angle=90),
        axis.text.y=element_text(size = 8))





################################### prognostic value of single metabolite in cancers
library(forestplot)
dim(sum)
head(sum)
sum2 <- sum

sum2$HR1 <- paste(round(sum2$HR,3),"(",round(sum2$HR.confint.lower,3),"-",round(sum2$HR.confint.upper,3),")",sep = '')
sum2$pvalue <- round(sum2$sctest.p,3)
dim(sum2)
head(sum2)

setwd("C:/MyRdata8/Summary/single metabolite on univariate cox analysis")


name1 <- unique(sum2$name)
length(name1)
i=1

sum3 <- sum2[which(sum2$name==name1[i]),]
dim(sum3)
head(sum3)

labeltext <- as.matrix(rbind(c("Tumor","HR","pvalue"),sum3[,c(24,25,26)]))
png(filename=paste("Forest plot on hazard ratio of ",name1[i]," in tumors.png",sep = ''),
    width=4000, height=2000,res=300)
forestplot(labeltext, 
           mean = c(NA,sum3$HR),
           lower = c(NA,sum3$HR.confint.lower),
           upper = c(NA,sum3$HR.confint.upper),
           zero = 1,lwd.zero = 2,
           clip = c(0,10),
           graph.pos = 2,
           boxsize = 0.2,
           align = "l",
           xlog = FALSE,
           col=fpColors(line = "black", zero = "black",
                        box="black"),
           lty.ci = 7,
           lwd.ci = 1,
           ci.vertices.height = 0.05,
           lwd.xaxis=2,
           title=paste("Forest plot on hazard ratio of ",name1[i]," in tumors",sep = ''),
           xlab="Hazard Ratio",
           txt_gp = fpTxtGp(ticks = gpar(cex = 0.9), xlab = gpar(cex = 0.8), cex = 0.8),
           lineheight = unit(8,"mm")) 
dev.off()
i=i+1
i


########################### pathway-metabolite-prognosis network/circos plot
setwd("C:/MyRdata8/Summary")
dim(sum1)
head(sum1)
sum1$pvalue2 <- (-log10(sum1$pvalue))
head(sum1)
dim(sum1)


### circos plot
library(dplyr)
library(circlize)

setwd("C:/MyRdata8/Summary")
sum <- read.csv("summary on single metabolite from univariate cox analysis5.csv")
head(sum)
dim(sum)
sum1 <- sum[,c(15,25,4)]
dim(sum1)
head(sum1)
sum1$color <- ifelse(sum1$HR<1,"red","blue")
sum1$value <- rep(1,nrow(sum1))

head(sum1)
sum2 <- sum1[,-3]
head(sum2)

tumor <- data.frame(name=names(by(sum2,sum2$Tumor,count)),count=as.numeric(by(sum2,sum2$Tumor,count)))
tumor <- tumor[order(tumor$count),]
tumor

nucleotide1 <- data.frame(name=names(by(sum2,sum2$name,count)),count=as.numeric(by(sum2,sum2$name,count)))
nucleotide1 <- nucleotide1[order(nucleotide1$count),]
dim(nucleotide1)
head(nucleotide1)

grid.col <- c(1:12,rep("grey",nrow(nucleotide1)))
names(grid.col) <- c(tumor$name,nucleotide1$name)
grid.col

par(mfrow = c(1, 1))
chordDiagram(sum2,order = c(tumor$name,nucleotide1$name),
             annotationTrack = "grid", 
             preAllocateTracks = list(track.height = min(strwidth(c(tumor$name,nucleotide1$name)))),
             directional = 1,
             grid.col=grid.col,
             col = sum1$color,
             transparency = 0.5)
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", cex = 0.5,niceFacing = TRUE, adj = c(0, 0.5))}, bg.border = NA)
circos.clear()






############################ summary on multivariate cox analysis results ############################
setwd("C:/MyRdata8/Summary")

sum <- rbind(multivariate.lung.primary.tumor.n85.sig,multivariate.lung.NAC.tumor.n77.sig,multivariate.lunge.plattenepithel.n238.sig,
             multivariate.all.primary.EAC.n102.sig,multivariate.Neoadjuvant.treated.EACs.n144.sig,
             multivariate.Chromophobe.RCC.n100.sig,multivariate.Clear.cell.RCC.n476.sig,
             multivariate.Papillary.RCC.n101.sig,multivariate.ACC.n72.sig,multivariate.Pancreas.tumor.n107.sig,
             multivariate.Gastric.Cancer.n246.sig,multivariate.HER2.n82.sig)
sum$Tumor <- c(rep("Lung Primary Tumor",nrow(multivariate.lung.primary.tumor.n85.sig)),
               rep("Lung NAC Tumor",nrow(multivariate.lung.NAC.tumor.n77.sig)),
               rep("Lung squamous cell carcinoma",nrow(multivariate.lunge.plattenepithel.n238.sig)),
               rep("Primary EAC",nrow(multivariate.all.primary.EAC.n102.sig)),
               rep("Neoadjuvant treated EAC",nrow(multivariate.Neoadjuvant.treated.EACs.n144.sig)),
               rep("Chromophobe RCC",nrow(multivariate.Chromophobe.RCC.n100.sig)),
               rep("Clear cell RCC",nrow(multivariate.Clear.cell.RCC.n476.sig)),
               rep("Papillary RCC",nrow(multivariate.Papillary.RCC.n101.sig)),
               rep("Adrenocortical Carcinoma",nrow(multivariate.ACC.n72.sig)),
               rep("Pancreatic cancer",nrow(multivariate.Pancreas.tumor.n107.sig)),
               rep("Gastric cancer",nrow(multivariate.Gastric.Cancer.n246.sig)),
               rep("HER2-Gastric cancer",nrow(multivariate.HER2.n82.sig)))

dim(sum)
head(sum)
setwd("C:/MyRdata8/Summary")
#write.csv(sum,"summary on single metabolite from multivariate cox analysis5.csv")
sum <- read.csv("summary on single metabolite from multivariate cox analysis5.csv")
sum <- sum[,-1]
dim(sum)
head(sum)


####################################### heatmap by ggplot2
library(ggplot2)
sum1 <- data.frame(name=c(sum$name),
                   Tumor=c(sum$Tumor),
                   Pathway=c(sum$Pathway),
                   mass=c(sum$query_mass),
                   HR=c(sum[,2]),
                   pvalue=c(sum[,3]))
dim(sum1)
head(sum1)
sum1$Prognosis <- ifelse(sum1$HR<1,"Good","Poor")
dim(sum1)
head(sum1)
sum2 <- sum1[,c(1,2,7)]
head(sum2)
dim(sum2)
sum2$name1 <- paste(sum2$Tumor,sum2$name,sep = '.')


sum6 <- sum5.5
head(sum6)
colnames(sum6)[1] <- "Tumor"
head(sum6)
dim(sum6)
sum6$name1 <- paste(sum6$Tumor,sum6$name,sep = '.')
sum6$value <- sum2$Prognosis[match(sum6$name1,sum2$name1)]
sum6$value[which(is.na(sum6$value)==T)] <- "no"

ggplot(sum6,aes(x=name,y=Tumor))+
  geom_tile(aes(fill=value),color="black")+
  scale_fill_manual(values = c("green","white","red"))+
  xlab("Nucleotide")+ylab("Cancer")+guides(fill=F)+
  theme(axis.text.x = element_text(size = 6,vjust=0.5,hjust=1,angle=90),
        axis.text.y=element_text(size = 8)) # width:1200,hight:400



################################### prognostic value of single metabolite in cancers
library(forestplot)
dim(sum)
head(sum)
sum2 <- sum

sum2$HR1 <- paste(round(sum2$exp.coef.,3),"(",round(sum2$lower..95,3),"-",round(sum2$upper..95,3),")",sep = '')
sum2$pvalue <- round(sum2$Pr...z..,3)
dim(sum2)
head(sum2)

setwd("C:/MyRdata8/Summary/single metabolite on multivariate cox analysis")

name1 <- unique(sum2$name)
length(name1)
i=18

sum3 <- sum2[which(sum2$name==name1[i]),]
dim(sum3)
head(sum3)

labeltext <- as.matrix(rbind(c("Tumor","HR","pvalue"),sum3[,c(28,29,30)]))
png(filename=paste("Forest plot on hazard ratio of ",name1[i]," in tumors.png",sep = ''),
    width=4000, height=2000,res=300)
forestplot(labeltext, 
           mean = c(NA,sum3$exp.coef.),
           lower = c(NA,sum3$lower..95),
           upper = c(NA,sum3$upper..95),
           zero = 1,lwd.zero = 2,
           clip = c(0,10),
           graph.pos = 2,
           boxsize = 0.2,
           align = "l",
           xlog = FALSE,
           col=fpColors(line = "black", zero = "black",
                        box="black"),
           lty.ci = 7,
           lwd.ci = 1,
           ci.vertices.height = 0.05,
           lwd.xaxis=2,
           title=paste("Forest plot on hazard ratio of ",name1[i]," in tumors",sep = ''),
           xlab="Hazard Ratio",
           txt_gp = fpTxtGp(ticks = gpar(cex = 0.9), xlab = gpar(cex = 0.8), cex = 0.8),
           lineheight = unit(8,"mm")) 
dev.off()
i=i+1
i



forestplot(labeltext, 
           mean = c(NA,sum3$exp.coef.),
           lower = c(NA,sum3$lower..95),
           upper = c(NA,sum3$upper..95),
           zero = 1,lwd.zero = 2,
           clip = c(0,10),
           graph.pos = 2,
           boxsize = 0.2,
           align = "l",
           xlog = FALSE,
           col=fpColors(line = "black", zero = "black",
                        box="black"),
           lty.ci = 7,
           lwd.ci = 1,
           ci.vertices.height = 0.05,
           lwd.xaxis=2,
           title=paste("Forest plot on hazard ratio of ",name1[i]," in tumors",sep = ''),
           xlab="Hazard Ratio",
           txt_gp = fpTxtGp(ticks = gpar(cex = 0.9), xlab = gpar(cex = 0.8), cex = 0.9),
           lineheight = unit(8,"mm")) 



########################### pathway-metabolite-prognosis network/circos plot

head(sum1)
dim(sum1)
#write.csv(sum1,"All significant nucleotides for prognosis in various tumors from multivariate cox analysis.csv")


### circos plot
library(dplyr)
library(circlize)

setwd("C:/MyRdata8/Summary")
sum1 <- read.csv("All significant nucleotides for prognosis in various tumors from multivariate cox analysis.csv")
sum1 <- sum1[,-1]
dim(sum1)
head(sum1)
sum1$color <- ifelse(sum1$Prognosis=="Good","green","red")
sum1$value <- rep(1,nrow(sum1))

sum2 <- sum1[,c(2,1,9)]
head(sum2)
dim(sum2)

tumor <- data.frame(name=names(by(sum2,sum2$Tumor,count)),count=as.numeric(by(sum2,sum2$Tumor,count)))
tumor <- tumor[order(tumor$count),]
tumor

nucleotide1 <- data.frame(name=names(by(sum2,sum2$name,count)),count=as.numeric(by(sum2,sum2$name,count)))
nucleotide1 <- nucleotide1[order(nucleotide1$count),]
dim(nucleotide1)
head(nucleotide1)


grid.col <- c("black","red","blue","yellow","green","orange",
              "purple","brown","cyan","grey", "yellowgreen",
              rep("grey",nrow(nucleotide1)))
names(grid.col) <- c(tumor$name,nucleotide1$name)
grid.col

par(mfrow = c(1, 1))
chordDiagram(sum2,order = c(tumor$name,nucleotide1$name),
             annotationTrack = "grid", 
             preAllocateTracks = list(track.height = min(strwidth(c(tumor$name,nucleotide1$name)))),
             directional = -1,
             grid.col=grid.col,
             col = sum1$color,
             transparency = 0.5) 
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", cex = 0.7,niceFacing = TRUE, adj = c(0, 0.5))}, bg.border = NA)
circos.clear()

rand_color(3)
head(colors(),657)



########################### overview on diagnostic/predictive performance of cox model
setwd("C:/MyRdata8/Summary")
library(ggplot2)
library(reshape)

sum <- rbind(prognosis.n82,prognosis.n246,prognosis.n107,prognosis.n72,prognosis.n101,
             prognosis.n476,prognosis.n100,prognosis.n144,prognosis.n102,prognosis.n238,
             prognosis.n77,prognosis.n85)

dim(sum)
head(sum)
#write.csv(sum,"overview on diagnostic predictive performance of cox model.csv")

sum$Tumor <- factor(sum$Tumor,levels = c("Lung Primary Tumor","Lung NAC Tumor","Lung squamous cell carcinoma",
                                         "Primary EAC","Neoadjuvant treated EAC","Chromophobe RCC",
                                         "Clear cell RCC","Papillary RCC","Adrenocortical Carcinoma",
                                         "Pancreatic cancer","Gastric cancer","HER2-Gastric cancer"))

ggplot(data=sum,aes(x=Tumor,y=AUC,group=Prognosis,colour = Prognosis))+geom_line()+geom_point()+
  ylab("AUC value")+
  theme(axis.text.x = element_text(hjust=1,vjust=0.5,angle=90))




####################### summary independent nucleotide-involved drug pathway based on SMPDB database ########################

### SMPDB pathway
setwd("C:/MyRdata8/Summary/SMPDB/smpdb_pathways.csv")

smpdb_pathways <- read.csv("smpdb_pathways.csv")
dim(smpdb_pathways)
head(smpdb_pathways)
factor(smpdb_pathways$Subject)

drug <- smpdb_pathways[smpdb_pathways$Subject=="Drug Action"|smpdb_pathways$Subject=="Drug Metabolism",]
dim(drug)
head(drug)
#write.csv(drug,"all drug-related pathway.csv")

### cancer related drug pathway
index <- c()
for (i in 1:nrow(drug)) {
  if (("cancer" %in% strsplit(drug$Description[i],' ')[[1]])|("cancers" %in% strsplit(drug$Description[i],' ')[[1]])|("anticancer" %in% strsplit(drug$Description[i],' ')[[1]])) {
    index <- c(index,i)
  }
}
#drug.cancer <- drug[-index,] 
drug.cancer <- drug[index,] 
dim(drug.cancer)
head(drug.cancer)
#drug.cancer <- drug

### all anticancer drug related pathways and metabolites
setwd("C:/MyRdata8/Summary/SMPDB/smpdb_metabolites.csv")

drug1 <- read.csv(paste(drug.cancer$SMPDB.ID[1],"_metabolites.csv",sep = ''))
for (i in 2:length(drug.cancer$SMPDB.ID)) {
  drug2 <- read.csv(paste(drug.cancer$SMPDB.ID[i],"_metabolites.csv",sep = ''))
  drug1 <- rbind(drug1,drug2)
}

dim(drug1)
head(drug1)


### drug pathway related Nucleotide
setwd("C:/MyRdata8/Summary/SMPDB")

nucleotide2 <- read.csv("Nucleotide metabolism1.csv")
dim(nucleotide2)
head(nucleotide2)

nucleotide.drug1 <- data.frame()
for (i in 1:131) {
  nucleotide.drug <- drug1[which(drug1$HMDB.ID==nucleotide2$HMDB.ID[i]),]
  nucleotide.drug1 <- rbind(nucleotide.drug1,nucleotide.drug)
}
dim(nucleotide.drug1)
head(nucleotide.drug1)

nucleotide.drug1$KEGG.pathway <- nucleotide2[match(nucleotide.drug1$HMDB.ID,nucleotide2$HMDB.ID),2]
nucleotide.drug1$name <- nucleotide2[match(nucleotide.drug1$HMDB.ID,nucleotide2$HMDB.ID),4]
dim(nucleotide.drug1)

nucleotide.drug2 <- nucleotide.drug1[,c(1:5,17,16,6:11)]
dim(nucleotide.drug2)
head(nucleotide.drug2)
length(unique(c(nucleotide.drug2$HMDB.ID))) #5
length(unique(c(nucleotide.drug2$Pathway.Name))) #11
unique(c(nucleotide.drug2$name))


### nucleotide-drug pathway network
#write.csv(nucleotide.drug2,"60 Nucleotides involved drug pathways.csv")
#write.csv(nucleotide.drug2,"5 Nucleotides involved anticancer drug pathways.csv")


### drug pathways related to significant nucleotides from univariate analysis
setwd("C:/MyRdata8/Summary")
sum <- read.csv("summary on single metabolite from univariate cox analysis5.csv")
sum <- sum[,-1]
dim(sum)
head(sum)
colnames(sum)[10] <- "HMDB.ID"

com1 <- intersect(unique(c(sum$HMDB.ID)),unique(c(nucleotide.drug2$HMDB.ID)))
length(com1)
com1

sum1 <- merge(x=sum,y=nucleotide.drug2,by="HMDB.ID")
dim(sum1)
head(sum1)
length(unique(c(sum1$HMDB.ID)))
length(unique(c(sum1$Pathway.Name)))
length(unique(c(sum1$Tumor)))
#write.csv(sum1,"anticancer drug pathways related to significant nucleotides from univariate analysis.csv")


### drug pathways related to independent nucleotides from multivariate analysis
setwd("C:/MyRdata8/Summary")
sum5 <- read.csv("summary on single metabolite from multivariate cox analysis5.csv")
sum5 <- sum5[,-1]
colnames(sum5)[14] <- "HMDB.ID"
head(sum5)
dim(sum5)

com1 <- intersect(unique(c(sum5$HMDB.ID)),unique(c(nucleotide.drug2$HMDB.ID)))
length(com1)

sum6 <- merge(x=sum5,y=nucleotide.drug2,by="HMDB.ID")
dim(sum6)
head(sum6)
length(unique(c(sum6$HMDB.ID)))
length(unique(c(sum6$Pathway.Name)))
length(unique(c(sum6$Tumor)))
#write.csv(sum6,"anticancer drug pathways related to independent nucleotides from multivariate analysis.csv")



########## circular plot shows drug pathway analysis on independent prognostic factor
setwd("C:/MyRdata8/Summary")

pathway <- read.csv("anticancer drug pathways related to independent nucleotides from multivariate analysis1.csv")
dim(pathway)
head(pathway)
pathway$name.x[which(pathway$name.x=="3,5-ADP")] <- "ADP"
pathway <- distinct(pathway)

pathway1 <- pathway[,c(1,2,6)]
head(pathway1)
dim(pathway1)


library(circlize)

#rand_color(3)
grid.col <- structure(c(rep(1, 3), rep(2, 11), rep(3, 6),rep(4,2)),
                      names = c(unique(pathway1$name.x),unique(pathway1$name.y)[c(1:6,18:19)])[c(1:14,19,20,15,16,17:18,21,22)])
grid.col

group <- grid.col
group

chordDiagram(pathway1,order =c(unique(pathway1$name.x),unique(pathway1$name.y)[c(1:6,18:19)])[c(1:14,19,20,15,16,17:18,21,22)],
             annotationTrack = "grid", 
             preAllocateTracks = list(track.height = min(strwidth(c(pathway1$name.x,pathway1$name.y)))),
             directional = 0,
             group = group,
             grid.col=grid.col,
             col = pathway$color,
             transparency = 0.5)
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise",  cex = 1,niceFacing = TRUE, adj = c(0, 0.5))}, bg.border = NA)
circos.clear()




############################################## co-expressed network of nucleotide related to drug pathway

setwd("C:/Database1 HMDB annotation")
library(dplyr)
annotation <- read.csv("355372 Annotation on all mass values.csv")
dim(annotation)
head(annotation)
annotation <- annotation[,-1]
dim(annotation)
head(annotation)



######### ACC
dim(table3.n72)
table3.n72[1:5,1:5]

# median intensity on each mass
median1 <- data.frame(intensity=apply(table3.n72,MARGIN=2,FUN=median))
median1$query_mass <- round(as.numeric(substr(rownames(median1),2,12)),4)
head(median1)
dim(median1)

cor1 <- cor(table3.n72)
dim(cor1)
cor1[1:5,1:5]

dim(annotation)
head(annotation)

mass1 <- c("X408.0127") #adp
mass1 <- c("X424.9794") #udp

nucleotide2 <-data.frame(corr=cor1[,mass1])
dim(nucleotide2)
head(nucleotide2)
nucleotide2$query_mass <- round(as.numeric(substr(rownames(nucleotide2),2,12)),4)
nucleotide2 <- nucleotide2[abs(nucleotide2$corr)>0.5&abs(nucleotide2$corr)<1,]

nucleotide3 <- merge(nucleotide2,annotation)
nucleotide3 <- nucleotide3[which(is.na(nucleotide3$compound_id)==F),]
dim(nucleotide3)
nucleotide4 <- merge(median1,nucleotide3)
dim(nucleotide4)
head(nucleotide4)
nucleotide4$name.adduct <- paste(nucleotide4$adduct,nucleotide4$compound_id,sep = '.')
nucleotide4 <- nucleotide4[order(nucleotide4$ppm,decreasing = F),]
nucleotide4 <- nucleotide4[which(duplicated(nucleotide4$name.adduct)==F),]
nucleotide4 <- nucleotide4[order(nucleotide4$intensity,decreasing = T),]
nucleotide4 <- nucleotide4[which(duplicated(nucleotide4$compound_id)==F),]
head(nucleotide4)
dim(nucleotide4)

table4 <- table3.n72[,match(nucleotide4$query_mass,round(as.numeric(substr(colnames(table3.n72),2,12)),4))]
colnames(table4) <- nucleotide4$compound_id
dim(table4)
table4[1:5,1:5]

adp.n72 <- table4
udp.n72 <- table4

#write.csv(colnames(adp.n72),"pathway analysis on adp.n72.csv")
#write.csv(colnames(udp.n72),"pathway analysis on udp.n72.csv")


######### Clear cell RCC
dim(table3.n476)
table3.n476[1:5,1:5]

# median intensity on each mass
median1 <- data.frame(intensity=apply(table3.n476,MARGIN=2,FUN=median))
median1$query_mass <- round(as.numeric(substr(rownames(median1),2,12)),4)
head(median1)
dim(median1)

cor1 <- cor(table3.n476)
dim(cor1)
cor1[1:5,1:5]

dim(annotation)
head(annotation)

mass1 <- c("X426.0220333") #adp

nucleotide2 <-data.frame(corr=cor1[,mass1])
dim(nucleotide2)
head(nucleotide2)
nucleotide2$query_mass <- round(as.numeric(substr(rownames(nucleotide2),2,12)),4)
nucleotide2 <- nucleotide2[abs(nucleotide2$corr)>0.5&abs(nucleotide2$corr)<1,]

nucleotide3 <- merge(nucleotide2,annotation)
nucleotide3 <- nucleotide3[which(is.na(nucleotide3$compound_id)==F),]
dim(nucleotide3)
nucleotide4 <- merge(median1,nucleotide3)
dim(nucleotide4)
head(nucleotide4)
nucleotide4$name.adduct <- paste(nucleotide4$adduct,nucleotide4$compound_id,sep = '.')
nucleotide4 <- nucleotide4[order(nucleotide4$ppm,decreasing = F),]
nucleotide4 <- nucleotide4[which(duplicated(nucleotide4$name.adduct)==F),]
nucleotide4 <- nucleotide4[order(nucleotide4$intensity,decreasing = T),]
nucleotide4 <- nucleotide4[which(duplicated(nucleotide4$compound_id)==F),]
head(nucleotide4)
dim(nucleotide4)

table4 <- table3.n476[,match(nucleotide4$query_mass,round(as.numeric(substr(colnames(table3.n476),2,12)),4))]
colnames(table4) <- nucleotide4$compound_id
dim(table4)
table4[1:5,1:5]

adp.n476 <- table4
#write.csv(colnames(adp.n476),"pathway analysis on adp.n476.csv")



######## HER2 Test Deviations 
dim(table3.n82)
table3.n82[1:5,1:5]

# median intensity on each mass
median1 <- data.frame(intensity=apply(table3.n82,MARGIN=2,FUN=median))
median1$query_mass <- round(as.numeric(substr(rownames(median1),2,12)),4)
head(median1)
dim(median1)

cor1 <- cor(table3.n82)
dim(cor1)
cor1[1:5,1:5]

mass1 <- c("X408.0117") #adp

nucleotide2 <-data.frame(corr=cor1[,mass1])
dim(nucleotide2)
head(nucleotide2)
nucleotide2$query_mass <- round(as.numeric(substr(rownames(nucleotide2),2,12)),4)
nucleotide2 <- nucleotide2[abs(nucleotide2$corr)>0.5&abs(nucleotide2$corr)<1,]

nucleotide3 <- merge(nucleotide2,annotation)
nucleotide3 <- nucleotide3[which(is.na(nucleotide3$compound_id)==F),]
dim(nucleotide3)
nucleotide4 <- merge(median1,nucleotide3)
dim(nucleotide4)
head(nucleotide4)
nucleotide4$name.adduct <- paste(nucleotide4$adduct,nucleotide4$compound_id,sep = '.')
nucleotide4 <- nucleotide4[order(nucleotide4$ppm,decreasing = F),]
nucleotide4 <- nucleotide4[which(duplicated(nucleotide4$name.adduct)==F),]
nucleotide4 <- nucleotide4[order(nucleotide4$intensity,decreasing = T),]
nucleotide4 <- nucleotide4[which(duplicated(nucleotide4$compound_id)==F),]
head(nucleotide4)
dim(nucleotide4)

table4 <- table3.n82[,match(nucleotide4$query_mass,round(as.numeric(substr(colnames(table3.n82),2,12)),4))]
colnames(table4) <- nucleotide4$compound_id
dim(table4)
table4[1:5,1:5]

adp.n82 <- table4
#write.csv(colnames(adp.n82),"pathway analysis on adp.n82.csv")



########## Lung primary tumor (n=85) 
dim(table3.n85)
table3.n85[1:5,1:5]

# median intensity on each mass
median1 <- data.frame(intensity=apply(table3.n85,MARGIN=2,FUN=median))
median1$query_mass <- round(as.numeric(substr(rownames(median1),2,12)),4)
head(median1)
dim(median1)

cor1 <- cor(table3.n85)
dim(cor1)
cor1[1:5,1:5]

mass1 <- c("X424.977") #udp
mass1 <- c("X505.9817") #paps

nucleotide2 <-data.frame(corr=cor1[,mass1])
dim(nucleotide2)
head(nucleotide2)
nucleotide2$query_mass <- round(as.numeric(substr(rownames(nucleotide2),2,12)),4)
nucleotide2 <- nucleotide2[abs(nucleotide2$corr)>0.4&abs(nucleotide2$corr)<1,]

nucleotide3 <- merge(nucleotide2,annotation)
nucleotide3 <- nucleotide3[which(is.na(nucleotide3$compound_id)==F),]
dim(nucleotide3)
nucleotide4 <- merge(median1,nucleotide3)
dim(nucleotide4)
head(nucleotide4)
nucleotide4$name.adduct <- paste(nucleotide4$adduct,nucleotide4$compound_id,sep = '.')
nucleotide4 <- nucleotide4[order(nucleotide4$ppm,decreasing = F),]
nucleotide4 <- nucleotide4[which(duplicated(nucleotide4$name.adduct)==F),]
nucleotide4 <- nucleotide4[order(nucleotide4$intensity,decreasing = T),]
nucleotide4 <- nucleotide4[which(duplicated(nucleotide4$compound_id)==F),]
head(nucleotide4)
dim(nucleotide4)

table4 <- table3.n85[,match(nucleotide4$query_mass,round(as.numeric(substr(colnames(table3.n85),2,12)),4))]
colnames(table4) <- nucleotide4$compound_id
dim(table4)
table4[1:5,1:5]

udp.n85 <- table4
paps.n85 <- table4
#write.csv(colnames(udp.n85),"pathway analysis on udp.n85.csv")
#write.csv(colnames(paps.n85),"pathway analysis on paps.n85.csv")


############ lunge Plattenepithel (n=238)
dim(table3.n238)
table3.n238[1:5,1:5]

# median intensity on each mass
median1 <- data.frame(intensity=apply(table3.n238,MARGIN=2,FUN=median))
median1$query_mass <- round(as.numeric(substr(rownames(median1),2,12)),4)
head(median1)
dim(median1)

cor1 <- cor(table3.n238)
dim(cor1)
cor1[1:5,1:5]

mass1 <- c("X384.9846068") #udp

nucleotide2 <-data.frame(corr=cor1[,mass1])
dim(nucleotide2)
head(nucleotide2)
nucleotide2$query_mass <- round(as.numeric(substr(rownames(nucleotide2),2,12)),4)
nucleotide2 <- nucleotide2[abs(nucleotide2$corr)>0.5&abs(nucleotide2$corr)<1,]

nucleotide3 <- merge(nucleotide2,annotation)
nucleotide3 <- nucleotide3[which(is.na(nucleotide3$compound_id)==F),]
dim(nucleotide3)
nucleotide4 <- merge(median1,nucleotide3)
dim(nucleotide4)
head(nucleotide4)
nucleotide4$name.adduct <- paste(nucleotide4$adduct,nucleotide4$compound_id,sep = '.')
nucleotide4 <- nucleotide4[order(nucleotide4$ppm,decreasing = F),]
nucleotide4 <- nucleotide4[which(duplicated(nucleotide4$name.adduct)==F),]
nucleotide4 <- nucleotide4[order(nucleotide4$intensity,decreasing = T),]
nucleotide4 <- nucleotide4[which(duplicated(nucleotide4$compound_id)==F),]
head(nucleotide4)
dim(nucleotide4)

table4 <- table3.n238[,match(nucleotide4$query_mass,round(as.numeric(substr(colnames(table3.n238),2,12)),4))]
colnames(table4) <- nucleotide4$compound_id
dim(table4)
table4[1:5,1:5]

udp.n238 <- table4
#write.csv(colnames(udp.n238),"pathway analysis on udp.n238.csv")

########## Papillary RCC (n=101)
dim(table3.n101)
table3.n101[1:5,1:5]

# median intensity on each mass
median1 <- data.frame(intensity=apply(table3.n101,MARGIN=2,FUN=median))
median1$query_mass <- round(as.numeric(substr(rownames(median1),2,12)),4)
head(median1)
dim(median1)

cor1 <- cor(table3.n101)
dim(cor1)
cor1[1:5,1:5]

mass1 <- c("X384.9842") #udp

nucleotide2 <-data.frame(corr=cor1[,mass1])
dim(nucleotide2)
head(nucleotide2)
nucleotide2$query_mass <- round(as.numeric(substr(rownames(nucleotide2),2,12)),4)
nucleotide2 <- nucleotide2[abs(nucleotide2$corr)>0.5&abs(nucleotide2$corr)<1,]

nucleotide3 <- merge(nucleotide2,annotation)
nucleotide3 <- nucleotide3[which(is.na(nucleotide3$compound_id)==F),]
dim(nucleotide3)
nucleotide4 <- merge(median1,nucleotide3)
dim(nucleotide4)
head(nucleotide4)
nucleotide4$name.adduct <- paste(nucleotide4$adduct,nucleotide4$compound_id,sep = '.')
nucleotide4 <- nucleotide4[order(nucleotide4$ppm,decreasing = F),]
nucleotide4 <- nucleotide4[which(duplicated(nucleotide4$name.adduct)==F),]
nucleotide4 <- nucleotide4[order(nucleotide4$intensity,decreasing = T),]
nucleotide4 <- nucleotide4[which(duplicated(nucleotide4$compound_id)==F),]
head(nucleotide4)
dim(nucleotide4)

table4 <- table3.n101[,match(nucleotide4$query_mass,round(as.numeric(substr(colnames(table3.n101),2,12)),4))]
colnames(table4) <- nucleotide4$compound_id
dim(table4)
table4[1:5,1:5]

udp.n101 <- table4
#write.csv(colnames(udp.n101),"pathway analysis on udp.n101.csv")


########## circular plot shows pathway analysis on expression correlated metabolites
setwd("C:/MyRdata8/Summary")

pathway1 <- read.csv("circular plot shows pathway analysis on expression correlated metabolites1.csv")
dim(pathway1)
head(pathway1)

library(circlize)

#rand_color(3)
grid.col = structure(c(1:8,"#FDC0C3FF",rep("#129488FF", 12)),
                     names = c(unique(pathway1$Nucleotide),unique(pathway1$Pathway)[1:12]))
grid.col

group = structure(c(rep(2, 3), rep(3, 6), rep(4, 12)),
                  names = c(unique(pathway1$Nucleotide),unique(pathway1$Pathway)[1:12]))

chordDiagram(pathway1,
             annotationTrack = "grid", 
             preAllocateTracks = list(track.height = min(strwidth(c(pathway1$Nucleotide,pathway1$Pathway)))),
             directional = -1,
             group = group,
             grid.col=grid.col,
             transparency = 0.5)
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", cex=1,niceFacing = TRUE, adj = c(0, 0.5))}, bg.border = NA)
circos.clear()




####################### summary analysis on energy charge and prognosis using mixed adducts ##########################

setwd("C:/MyRdata8/Summary/Annotation")
annotation2 <- read.csv("355372 Annotation on all mass values.csv")
annotation2 <- annotation2[,-1]
dim(annotation2)
head(annotation2)

atp <- annotation2[which(annotation2$compound_id=="HMDB0000538"),]
dim(atp)
head(atp)

adp <- annotation2[which(annotation2$compound_id=="HMDB0001341"),]
dim(adp)
head(adp)

amp <- annotation2[which(annotation2$compound_id=="HMDB0000045"),]
dim(amp)
head(amp)


############################## Lunge Plattenepithel (n=238)
library(survival)
library(survminer)
setwd("C:/MyRdata8/Lunge Plattenepithel (n=238)")
table <- read.csv("Lunge Plattenepithel (n=238).csv")
dim(table)
table[1:5,1:10]

table1 <- table[,c(1,2,3,5,6,21,53:11414)]
dim(table1)
table1[1:5,1:10]
rownames(table1) <- table1[,1]
table2 <- table1[,-1]
dim(table2)
table2[1:5,1:7]

table3 <- table2[,-c(1:5)]
dim(table3)
table3[1:5,1:7]


#ATP
com1 <- intersect(round(as.numeric(substr(colnames(table3),2,10)),4),round(as.numeric(atp$query_mass),4))
com1

n238.atp <- table3[1:5,match(com1,round(as.numeric(substr(colnames(table3),2,10)),4))]
n238.atp #X505.9840315
atp.n238 <- atp[match(com1,round(as.numeric(atp$query_mass),4)),]
atp.n238

#ADP
com1 <- intersect(round(as.numeric(substr(colnames(table3),2,10)),4),round(as.numeric(adp$query_mass),4))
com1
com2 <- intersect(round(as.numeric(substr(colnames(table3),2,10)),3),round(as.numeric(adp$query_mass),3))
com2


n238.adp <- table3[1:5,match(com2,round(as.numeric(substr(colnames(table3),2,10)),3))]
n238.adp #X426.0225535
adp.n238 <- adp[match(com2,round(as.numeric(adp$query_mass),3)),]
adp.n238


#AMP
com1 <- intersect(round(as.numeric(substr(colnames(table3),2,10)),4),round(as.numeric(amp$query_mass),4))
com1
com2 <- intersect(round(as.numeric(substr(colnames(table3),2,10)),3),round(as.numeric(amp$query_mass),3))
com2

n238.amp <- table3[1:5,match(com2,round(as.numeric(substr(colnames(table3),2,10)),3))]
n238.amp #X368.0379518
amp.n238 <- amp[match(com2,round(as.numeric(amp$query_mass),3)),]
amp.n238


#cox analysis/survival curve
dim(table2)
table2[1:5,1:7]
table4 <- table2[,c("STATOS","OS","X505.9840315","X426.0225535","X368.0379518")]
head(table4)
table4$energy <- (table4$X505.9840315+(table4$X426.0225535)*0.5)/(table4$X505.9840315+table4$X426.0225535+table4$X368.0379518)
table4$energy[which(is.na(table4$energy)==T)] <- 0.5
table5 <- table4
dim(table5)
head(table5)
energy.n238 <- table5[,c(1,2,6)]
head(energy.n238)

res.cut <- surv_cutpoint(table5, time = "OS", event = "STATOS",
                         minprop = 0.1,
                         variables = colnames(table5)[3:ncol(table5)])
res.cat <- surv_categorize(res.cut,labels = c(0, 1))
dim(res.cat)
res.cat[1:5,1:6]
res.cat.n238 <- res.cat

result <- data.frame()
for (i in c(3:ncol(res.cat))) {
  res.cox <- coxph(Surv(OS,STATOS)~res.cat[,i],data =res.cat)
  summ <- summary(res.cox)
  hypo <- cox.zph(res.cox)
  hypo1 <- hypo$table[-nrow(hypo$table),]
  result1 <- data.frame(mz.value=colnames(res.cat)[i],
                        HR = signif(summ$coef[2], digits=3),
                        HR.confint.lower = signif(summ$conf.int[,"lower .95"], 3),
                        HR.confint.upper = signif(summ$conf.int[,"upper .95"],3),
                        logtest.p=summ$logtest[3],
                        waldtest.p=summ$waldtest[3],
                        sctest.p=summ$sctest[3],
                        hypo.p=hypo1[3])
  result <- rbind(result,result1)
}
result

fit <- survfit(Surv(OS,STATOS)~energy,data =res.cat)
ggsurvplot(fit, data = res.cat, conf.int = F, # 增加置信区间
           risk.table = TRUE, pval=T) #不同时间点风险人数表，增加P值


res.cox <- coxph(Surv(OS,STATOS)~energy,data =res.cat)
summary(res.cox)
cox.zph(res.cox)


############################### all primary EAC (n=102)
library(survival)
library(survminer)

setwd("C:/MyRdata8/All primary EAC (n=102)")
table <- read.csv("All primary EAC (n=102).csv")
dim(table)
table[1:5,1:15]

table1 <- table[,c(1,2,3,4,5,13:6959)]
dim(table1)
table1[1:5,1:7]

rownames(table1) <- table1[,1]
table2 <- table1[,-1]
dim(table2)
table2[1:5,1:8]

table3 <- table2[,-c(1:6)]
dim(table3)
table3[1:5,1:7]


#ATP
com1 <- intersect(round(as.numeric(substr(colnames(table3),2,10)),4),round(as.numeric(atp$query_mass),4))
com1
com2 <- intersect(round(as.numeric(substr(colnames(table3),2,10)),3),round(as.numeric(atp$query_mass),3))
com2


n102.atp <- table3[1:5,match(com1,round(as.numeric(substr(colnames(table3),2,10)),4))]
n102.atp #X543.9423328
atp.n102 <- atp[match(com1,round(as.numeric(atp$query_mass),4)),]
atp.n102

#ADP
com1 <- intersect(round(as.numeric(substr(colnames(table3),2,10)),4),round(as.numeric(adp$query_mass),4))
com1
com2 <- intersect(round(as.numeric(substr(colnames(table3),2,10)),3),round(as.numeric(adp$query_mass),3))
com2

n102.adp <- table3[1:5,match(com2,round(as.numeric(substr(colnames(table3),2,10)),3))]
n102.adp #X408.0119107
adp.n102 <- adp[match(com2,round(as.numeric(adp$query_mass),3)),]
adp.n102


#AMP
com1 <- intersect(round(as.numeric(substr(colnames(table3),2,10)),4),round(as.numeric(amp$query_mass),4))
com1
com2 <- intersect(round(as.numeric(substr(colnames(table3),2,10)),3),round(as.numeric(amp$query_mass),3))
com2

n102.amp <- table3[1:5,match(com2,round(as.numeric(substr(colnames(table3),2,10)),3))]
n102.amp #X346.0555279
amp.n102 <- amp[match(com2,round(as.numeric(amp$query_mass),3)),]
amp.n102


#cox analysis/survival curve
dim(table2)
table2[1:5,1:7]
table4 <- table2[,c("STATOS","OS","X543.9423328","X408.0119107","X346.0555279")]
head(table4)
table4$energy <- (table4$X543.9423328+(table4$X408.0119107)*0.5)/(table4$X543.9423328+table4$X408.0119107+table4$X346.0555279)
table4$energy
table5 <- table4
dim(table5)
head(table5)
energy.n102 <- table5[,c(1,2,6)]
head(energy.n102)

res.cut <- surv_cutpoint(table5, time = "OS", event = "STATOS",
                         minprop = 0.1,
                         variables = colnames(table5)[3:ncol(table5)])
res.cat <- surv_categorize(res.cut,labels = c(0, 1))
dim(res.cat)
res.cat[1:5,1:6]
res.cat.n102 <- res.cat

result <- data.frame()
for (i in c(3:ncol(res.cat))) {
  res.cox <- coxph(Surv(OS,STATOS)~res.cat[,i],data =res.cat)
  summ <- summary(res.cox)
  hypo <- cox.zph(res.cox)
  hypo1 <- hypo$table[-nrow(hypo$table),]
  result1 <- data.frame(mz.value=colnames(res.cat)[i],
                        HR = signif(summ$coef[2], digits=3),
                        HR.confint.lower = signif(summ$conf.int[,"lower .95"], 3),
                        HR.confint.upper = signif(summ$conf.int[,"upper .95"],3),
                        logtest.p=summ$logtest[3],
                        waldtest.p=summ$waldtest[3],
                        sctest.p=summ$sctest[3],
                        hypo.p=hypo1[3])
  result <- rbind(result,result1)
}
result


fit <- survfit(Surv(OS,STATOS)~energy,data =res.cat)
ggsurvplot(fit, data = res.cat, conf.int = F, # 增加置信区间
           risk.table = TRUE, pval=T) #不同时间点风险人数表，增加P值

res.cox <- coxph(Surv(OS,STATOS)~energy,data =res.cat)
summary(res.cox)
cox.zph(res.cox)

############################### Neoadjuvant treated EACs (n=144) 
library(survival)
library(survminer)

setwd("C:/MyRdata8/Neoadjuvant treated EACs (n=144)")

table <- read.csv("Neoadjuvant treated EACs (n=144).csv")
dim(table)
table[1:5,1:7]

table1 <- table
dim(table1)
table1[1:5,1:7]

rownames(table1) <- table1[,1]
table2 <- table1[,-1]
dim(table2)
table2[1:5,1:12]

table3 <- table2[,-c(1:11)]
dim(table3)
table3[1:5,1:7]



#ATP
com1 <- intersect(round(as.numeric(substr(colnames(table3),2,10)),4),round(as.numeric(atp$query_mass),4))
com1
com2 <- intersect(round(as.numeric(substr(colnames(table3),2,10)),3),round(as.numeric(atp$query_mass),3))
com2
com3 <- intersect(round(as.numeric(substr(colnames(table3),2,10)),2),round(as.numeric(atp$query_mass),2))
com3

n144.atp <- table3[1:5,match(com2,round(as.numeric(substr(colnames(table3),2,10)),3))]
n144.atp #X505.9845464
atp.n144 <- atp[match(com2,round(as.numeric(atp$query_mass),3)),]
atp.n144

#ADP
com1 <- intersect(round(as.numeric(substr(colnames(table3),2,10)),4),round(as.numeric(adp$query_mass),4))
com1
com2 <- intersect(round(as.numeric(substr(colnames(table3),2,10)),3),round(as.numeric(adp$query_mass),3))
com2

n144.adp <- table3[1:5,match(com2,round(as.numeric(substr(colnames(table3),2,10)),3))]
n144.adp #X448.0049327
adp.n144 <- adp[match(com2,round(as.numeric(adp$query_mass),3)),]
adp.n144


#AMP
com1 <- intersect(round(as.numeric(substr(colnames(table3),2,10)),4),round(as.numeric(amp$query_mass),4))
com1
com2 <- intersect(round(as.numeric(substr(colnames(table3),2,10)),3),round(as.numeric(amp$query_mass),3))
com2

n144.amp <- table3[1:5,match(com2,round(as.numeric(substr(colnames(table3),2,10)),3))]
n144.amp #X328.0461861
amp.n144 <- amp[match(com2,round(as.numeric(amp$query_mass),3)),]
amp.n144


#cox analysis/survival curve
dim(table2)
table2[1:5,1:7]
table4 <- table2[,c("STATOS","OS","X505.9845464","X448.0049327","X328.0461861")]
head(table4)
table4$energy <- (table4$X505.9845464+(table4$X448.0049327)*0.5)/(table4$X505.9845464+table4$X448.0049327+table4$X328.0461861)
table4$energy
table5 <- table4
dim(table5)
head(table5)
energy.n144 <- table5[,c(1,2,6)] 
head(energy.n144)

res.cut <- surv_cutpoint(table5, time = "OS", event = "STATOS",
                         minprop = 0.1,
                         variables = colnames(table5)[3:ncol(table5)])
res.cat <- surv_categorize(res.cut,labels = c(0, 1))
dim(res.cat)
res.cat[1:5,1:6]
res.cat.n144 <- res.cat

result <- data.frame()
for (i in c(3:ncol(res.cat))) {
  res.cox <- coxph(Surv(OS,STATOS)~res.cat[,i],data =res.cat)
  summ <- summary(res.cox)
  hypo <- cox.zph(res.cox)
  hypo1 <- hypo$table[-nrow(hypo$table),]
  result1 <- data.frame(mz.value=colnames(res.cat)[i],
                        HR = signif(summ$coef[2], digits=3),
                        HR.confint.lower = signif(summ$conf.int[,"lower .95"], 3),
                        HR.confint.upper = signif(summ$conf.int[,"upper .95"],3),
                        logtest.p=summ$logtest[3],
                        waldtest.p=summ$waldtest[3],
                        sctest.p=summ$sctest[3],
                        hypo.p=hypo1[3])
  result <- rbind(result,result1)
}
result

fit <- survfit(Surv(OS,STATOS)~energy,data =res.cat)
ggsurvplot(fit, data = res.cat, conf.int = F, # 增加置信区间
           risk.table = TRUE, pval=T) #不同时间点风险人数表，增加P值

res.cox <- coxph(Surv(OS,STATOS)~energy,data =res.cat)
summary(res.cox)
cox.zph(res.cox)


############################## Adrenocortical Carcinoma (ACC) (n=72)
library(survival)
library(survminer)

setwd("C:/MyRdata8/Adrenocortical Carcinoma (ACC) (n=72)")
table1 <- read.csv("Adrenocortical Carcinoma (ACC) (n=72).csv")
dim(table1)
table1[1:5,1:7]
rownames(table1) <- table1[,1]
table2 <- table1[,-1]
dim(table2)
table2[1:5,1:7]

table3 <- table2[,-c(1:2)]
dim(table3)
table3[1:5,1:7]


#ATP
com1 <- intersect(round(as.numeric(substr(colnames(table3),2,10)),4),round(as.numeric(atp$query_mass),4))
com1
com2 <- intersect(round(as.numeric(substr(colnames(table3),2,10)),3),round(as.numeric(atp$query_mass),3))
com2
com3 <- intersect(round(as.numeric(substr(colnames(table3),2,10)),2),round(as.numeric(atp$query_mass),2))
com3

n72.atp <- table3[1:5,match(com2,round(as.numeric(substr(colnames(table3),2,10)),3))]
n72.atp #X543.9427
atp.n72 <- atp[match(com2,round(as.numeric(atp$query_mass),3)),]
atp.n72

#ADP
com1 <- intersect(round(as.numeric(substr(colnames(table3),2,10)),4),round(as.numeric(adp$query_mass),4))
com1
com2 <- intersect(round(as.numeric(substr(colnames(table3),2,10)),3),round(as.numeric(adp$query_mass),3))
com2

n72.adp <- table3[1:5,match(com2,round(as.numeric(substr(colnames(table3),2,10)),3))]
n72.adp #X408.0127
adp.n72 <- adp[match(com2,round(as.numeric(adp$query_mass),3)),]
adp.n72


#AMP
com1 <- intersect(round(as.numeric(substr(colnames(table3),2,10)),4),round(as.numeric(amp$query_mass),4))
com1
com2 <- intersect(round(as.numeric(substr(colnames(table3),2,10)),3),round(as.numeric(amp$query_mass),3))
com2

n72.amp <- table3[1:5,match(com2,round(as.numeric(substr(colnames(table3),2,10)),3))]
n72.amp #X328.0462
amp.n72 <- amp[match(com2,round(as.numeric(amp$query_mass),3)),]
amp.n72


#cox analysis/survival curve
dim(table2)
table2[1:5,1:7]
table4 <- table2[,c("STATOS","OS","X543.9427","X408.0127","X328.046")]
head(table4)
table4$energy <- (table4$X543.9427+(table4$X408.0127)*0.5)/(table4$X543.9427+table4$X408.0127+table4$X328.046)
table4$energy
table5 <- table4
dim(table5)
head(table5)
energy.n72 <- table5[,c(1,2,6)]
head(energy.n72)

res.cut <- surv_cutpoint(table5, time = "OS", event = "STATOS",
                         minprop = 0.1,
                         variables = colnames(table5)[3:ncol(table5)])
res.cat <- surv_categorize(res.cut,labels = c(0, 1))
dim(res.cat)
res.cat[1:5,1:6]
res.cat.n72 <- res.cat

result <- data.frame()
for (i in c(3:ncol(res.cat))) {
  res.cox <- coxph(Surv(OS,STATOS)~res.cat[,i],data =res.cat)
  summ <- summary(res.cox)
  hypo <- cox.zph(res.cox)
  hypo1 <- hypo$table[-nrow(hypo$table),]
  result1 <- data.frame(mz.value=colnames(res.cat)[i],
                        HR = signif(summ$coef[2], digits=3),
                        HR.confint.lower = signif(summ$conf.int[,"lower .95"], 3),
                        HR.confint.upper = signif(summ$conf.int[,"upper .95"],3),
                        logtest.p=summ$logtest[3],
                        waldtest.p=summ$waldtest[3],
                        sctest.p=summ$sctest[3],
                        hypo.p=hypo1[3])
  result <- rbind(result,result1)
}
result

fit <- survfit(Surv(OS,STATOS)~energy,data =res.cat)
ggsurvplot(fit, data = res.cat, conf.int = F, # 增加置信区间
           risk.table = TRUE, pval=T) #不同时间点风险人数表，增加P值

res.cox <- coxph(Surv(OS,STATOS)~energy,data =res.cat)
summary(res.cox)
cox.zph(res.cox)

############################## Pancreas tumor (n=107)
library(survival)
library(survminer)

setwd("C:/MyRdata8/Pancreas tumor (n=107)")
table <- read.csv("Pancreas tumor (n=107).csv")
dim(table)
table[1:5,1:10]

table1 <- table[,c(1,2,3,8,9,14,15,17:5253)]
dim(table1)
table1[1:5,1:10]
rownames(table1) <- table1[,1]
table2 <- table1[,-1]
dim(table2)
table2[1:5,1:10]

table3 <- table2[,-c(1:6)]
dim(table3)
table3[1:5,1:7]

#ATP
com1 <- intersect(round(as.numeric(substr(colnames(table3),2,10)),4),round(as.numeric(atp$query_mass),4))
com1
com2 <- intersect(round(as.numeric(substr(colnames(table3),2,10)),3),round(as.numeric(atp$query_mass),3))
com2
com3 <- intersect(round(as.numeric(substr(colnames(table3),2,10)),2),round(as.numeric(atp$query_mass),2))
com3

n107.atp <- table3[1:5,match(com1,round(as.numeric(substr(colnames(table3),2,10)),4))]
n107.atp #X543.9405
atp.n107 <- atp[match(com1,round(as.numeric(atp$query_mass),4)),]
atp.n107

#ADP
com1 <- intersect(round(as.numeric(substr(colnames(table3),2,10)),4),round(as.numeric(adp$query_mass),4))
com1
com2 <- intersect(round(as.numeric(substr(colnames(table3),2,10)),3),round(as.numeric(adp$query_mass),3))
com2

n107.adp <- table3[1:5,match(com1,round(as.numeric(substr(colnames(table3),2,10)),4))]
n107.adp #X408.0125
adp.n107 <- adp[match(com1,round(as.numeric(adp$query_mass),4)),]
adp.n107


#AMP
com1 <- intersect(round(as.numeric(substr(colnames(table3),2,10)),4),round(as.numeric(amp$query_mass),4))
com1
com2 <- intersect(round(as.numeric(substr(colnames(table3),2,10)),3),round(as.numeric(amp$query_mass),3))
com2

n107.amp <- table3[1:5,match(com1,round(as.numeric(substr(colnames(table3),2,10)),4))]
n107.amp #X346.0555
amp.n107 <- amp[match(com1,round(as.numeric(amp$query_mass),4)),]
amp.n107


#cox analysis/survival curve
dim(table2)
table2[1:5,1:7]
table4 <- table2[,c("STATOS","OS","X543.9405","X408.0125","X346.0555")]
head(table4)
table4$energy <- (table4$X543.9405+(table4$X408.0125)*0.5)/(table4$X543.9405+table4$X408.0125+table4$X346.0555)
table4$energy
table5 <- table4[,-3]
dim(table5)
head(table5)
energy.n107 <- table5[,c(1,2,5)]
head(energy.n107)

res.cut <- surv_cutpoint(table5, time = "OS", event = "STATOS",
                         minprop = 0.2,
                         variables = colnames(table5)[3:ncol(table5)])
res.cat <- surv_categorize(res.cut,labels = c(0, 1))
dim(res.cat)
res.cat[1:5,1:5]
res.cat.n107 <- res.cat

result <- data.frame()
for (i in c(3:ncol(res.cat))) {
  res.cox <- coxph(Surv(OS,STATOS)~energy,data =res.cat)
  summ <- summary(res.cox)
  hypo <- cox.zph(res.cox)
  hypo1 <- hypo$table[-nrow(hypo$table),]
  result1 <- data.frame(mz.value=colnames(res.cat)[i],
                        HR = signif(summ$coef[2], digits=3),
                        HR.confint.lower = signif(summ$conf.int[,"lower .95"], 3),
                        HR.confint.upper = signif(summ$conf.int[,"upper .95"],3),
                        logtest.p=summ$logtest[3],
                        waldtest.p=summ$waldtest[3],
                        sctest.p=summ$sctest[3],
                        hypo.p=hypo1[3])
  result <- rbind(result,result1)
}
result

fit <- survfit(Surv(OS,STATOS)~energy,data =res.cat)
fit
ggsurvplot(fit, data = res.cat, conf.int = F, # 增加置信区间
           risk.table = TRUE, pval=T) #不同时间点风险人数表，增加P值

res.cox <- coxph(Surv(OS,STATOS)~energy,data =res.cat)
summary(res.cox)
cox.zph(res.cox)

############################ summary prognosis value of energy charge on cancers

### with Pancreas tumor
sum4 <- rbind(energy.n238,energy.n102,energy.n144,energy.n72,energy.n107)
dim(sum4)
head(sum4)

res.cut <- surv_cutpoint(sum4, time = "OS", event = "STATOS",
                         minprop = 0.1,
                         variables = "energy")
res.cat <- surv_categorize(res.cut,labels = c(0, 1))
dim(res.cat)
res.cat[1:5,1:3]

res.cox <- coxph(Surv(OS,STATOS)~energy,data =res.cat)
summary(res.cox)
cox.zph(res.cox)

fit <- survfit(Surv(OS,STATOS)~energy,data =res.cat)
fit
ggsurvplot(fit, data = res.cat, conf.int = F, # 增加置信区间
           risk.table = TRUE, pval=T) #不同时间点风险人数表，增加P值


### with Pancreas tumor
sum4 <- rbind(res.cat.n238[,c("OS","STATOS","energy")],res.cat.n102[,c("OS","STATOS","energy")],res.cat.n144[,c("OS","STATOS","energy")],res.cat.n72[,c("OS","STATOS","energy")],res.cat.n107[,c("OS","STATOS","energy")])
dim(sum4)
head(sum4)

res.cat <- sum4

res.cox <- coxph(Surv(OS,STATOS)~energy,data =res.cat)
summary(res.cox)
cox.zph(res.cox)

fit <- survfit(Surv(OS,STATOS)~energy,data =res.cat)
fit
ggsurvplot(fit, data = res.cat, conf.int = F, # 增加置信区间
           risk.table = TRUE, pval=T) #不同时间点风险人数表，增加P值


### without Pancreas tumor
sum5 <- rbind(energy.n238,energy.n102,energy.n144,energy.n72)
dim(sum5)
head(sum5)

res.cut <- surv_cutpoint(sum5, time = "OS", event = "STATOS",
                         minprop = 0.2,
                         variables = "energy")
res.cat <- surv_categorize(res.cut,labels = c(0, 1))
dim(res.cat)
res.cat[1:5,1:3]

res.cox <- coxph(Surv(OS,STATOS)~energy,data =res.cat)
summary(res.cox)
cox.zph(res.cox)

fit <- survfit(Surv(OS,STATOS)~energy,data =res.cat)
ggsurvplot(fit, data = res.cat, conf.int = F, # 增加置信区间
           risk.table = TRUE, pval=T) #不同时间点风险人数表，增加P值


### without Pancreas tumor
sum4 <- rbind(res.cat.n238[,c("OS","STATOS","energy")],res.cat.n102[,c("OS","STATOS","energy")],res.cat.n144[,c("OS","STATOS","energy")],res.cat.n72[,c("OS","STATOS","energy")])
dim(sum4)
head(sum4)

res.cat <- sum4

res.cox <- coxph(Surv(OS,STATOS)~energy,data =res.cat)
summary(res.cox)
cox.zph(res.cox)

fit <- survfit(Surv(OS,STATOS)~energy,data =res.cat)
fit
ggsurvplot(fit, data = res.cat, conf.int = F, # 增加置信区间
           risk.table = TRUE, pval=T) #不同时间点风险人数表，增加P值



####### Pan-cancer cluster analysis on presence or absence of 105 nucleotides in all patients and survival analysis #########

setwd("C:/MyRdata8/Summary")
sum2 <- read.csv("annotated and detectable 105 nucleotides in all cancers.csv")
sum3 <- sum2[,-1]
head(sum3)
dim(sum3)
length(unique(sum3$name))
unique(sum3$tumor)

sum7 <- data.frame()
for (i in 1:12) {
  sum4 <- sum3[which(sum3$tumor==unique(sum3$tumor)[i]),]
  sum4$name.adduct <- paste(sum4$adduct,sum4$name,sep = '.')
  sum4 <- sum4[order(sum4$ppm,decreasing = F),]
  head(sum4)
  dim(sum4)
  sum5 <- sum4[which(duplicated(sum4$name.adduct)==F),]
  dim(sum5)
  head(sum5)
  sum5 <- sum5[order(sum5$intensity,decreasing = T),]
  sum6 <- sum5[which(duplicated(sum5$name)==F),]
  dim(sum6)
  head(sum6)
  sum7 <- rbind(sum7,sum6)
}
dim(sum7)
head(sum7)


tumor <- c("Lung Primary Tumor","Lung NAC Tumor","Lung squamous cell carcinoma",
           "Primary EAC","Neoadjuvant treated EAC","Chromophobe RCC",
           "Clear cell RCC","Papillary RCC","Adrenocortical Carcinoma",
           "Pancreatic cancer","Gastric cancer","HER2-Gastric cancer")
tumor


table4.1 <- nucleotide.profile.lung.primary.tumor.n85[,sum7[which(sum7$tumor==tumor[1]),]$mass]
colnames(table4.1) <- sum7[which(sum7$tumor==tumor[1]),]$name
table4.1$ID <- rownames(table4.1)
dim(table4.1)

table4.2 <- nucleotide.profile.lung.NAC.tumor.n77[,sum7[which(sum7$tumor==tumor[2]),]$mass]
colnames(table4.2) <- sum7[which(sum7$tumor==tumor[2]),]$name
table4.2$ID <- rownames(table4.2)
dim(table4.2)

table4.3 <- nucleotide.profile.lunge.plattenepithel.n238[,sum7[which(sum7$tumor==tumor[3]),]$mass]
colnames(table4.3) <- sum7[which(sum7$tumor==tumor[3]),]$name
table4.3$ID <- rownames(table4.3)
dim(table4.3)

table4.4 <- nucleotide.profile.all.primary.EAC.n102[,sum7[which(sum7$tumor==tumor[4]),]$mass]
colnames(table4.4) <- sum7[which(sum7$tumor==tumor[4]),]$name
table4.4$ID <- rownames(table4.4)
dim(table4.4)

table4.5 <- nucleotide.profile.Neoadjuvant.treated.EACs.n144[,sum7[which(sum7$tumor==tumor[5]),]$mass]
colnames(table4.5) <- sum7[which(sum7$tumor==tumor[5]),]$name
table4.5$ID <- rownames(table4.5)
dim(table4.5)

table4.6 <- nucleotide.profile.Chromophobe.RCC.n100[,sum7[which(sum7$tumor==tumor[6]),]$mass]
colnames(table4.6) <- sum7[which(sum7$tumor==tumor[6]),]$name
table4.6$ID <- rownames(table4.6)
dim(table4.6)

table4.7 <- nucleotide.profile.Clear.cell.RCC.n476[,sum7[which(sum7$tumor==tumor[7]),]$mass]
colnames(table4.7) <- sum7[which(sum7$tumor==tumor[7]),]$name
table4.7$ID <- rownames(table4.7)
dim(table4.7)

table4.8 <- nucleotide.profile.Papillary.RCC.n101[,sum7[which(sum7$tumor==tumor[8]),]$mass]
colnames(table4.8) <- sum7[which(sum7$tumor==tumor[8]),]$name
table4.8$ID <- rownames(table4.8)
dim(table4.8)

table4.9 <- nucleotide.profile.ACC.n72[,sum7[which(sum7$tumor==tumor[9]),]$mass]
colnames(table4.9) <- sum7[which(sum7$tumor==tumor[9]),]$name
table4.9$ID <- rownames(table4.9)
dim(table4.9)

table4.10 <- nucleotide.profile.Pancreas.tumor.n107[,sum7[which(sum7$tumor==tumor[10]),]$mass]
colnames(table4.10) <- sum7[which(sum7$tumor==tumor[10]),]$name
table4.10$ID <- rownames(table4.10)
dim(table4.10)

table4.11 <- nucleotide.profile.Gastric.Cancer.n246[,sum7[which(sum7$tumor==tumor[11]),]$mass]
colnames(table4.11) <- sum7[which(sum7$tumor==tumor[11]),]$name
table4.11$ID <- rownames(table4.11)
dim(table4.11)

table4.12 <- nucleotide.profile.HER2.n82[,sum7[which(sum7$tumor==tumor[12]),]$mass]
colnames(table4.12) <- sum7[which(sum7$tumor==tumor[12]),]$name
table4.12$ID <- rownames(table4.12)
dim(table4.12)

table4.1 <- merge(table4.1,table4.2,all = T)
table4.1 <- merge(table4.1,table4.3,all = T)
table4.1 <- merge(table4.1,table4.4,all = T)
table4.1 <- merge(table4.1,table4.5,all = T)
table4.1 <- merge(table4.1,table4.6,all = T)
table4.1 <- merge(table4.1,table4.7,all = T)
table4.1 <- merge(table4.1,table4.8,all = T)
table4.1 <- merge(table4.1,table4.9,all = T)
table4.1 <- merge(table4.1,table4.10,all = T)
table4.1 <- merge(table4.1,table4.11,all = T)
table4.1 <- merge(table4.1,table4.12,all = T)
dim(table4.1)
table4.1[1:5,1:5]
table5 <- table4.1
rownames(table5) <- table5$ID
table6 <- table5[,-match("ID",colnames(table5))]
dim(table6)
table6[1:5,1:5]


clin <- rbind(res.cat.lung.primary.tumor.n85[,c("OS","STATOS")],res.cat.lung.NAC.tumor.n77[,c("OS","STATOS")],
              res.cat.lunge.plattenepithel.n238[,c("OS","STATOS")],res.cat.all.primary.EAC.n102[,c("OS","STATOS")],
              res.cat.Neoadjuvant.treated.EACs.n144[,c("OS","STATOS")],res.cat.Chromophobe.RCC.n100[,c("OS","STATOS")],
              res.cat.Clear.cell.RCC.n476[,c("OS","STATOS")],res.cat.Papillary.RCC.n101[,c("OS","STATOS")],
              res.cat.ACC.n72[,c("OS","STATOS")],res.cat.Pancreas.tumor.n107[,c("OS","STATOS")],
              res.cat.Gastric.Cancer.n246[,c("OS","STATOS")],res.cat.HER2.n82[,c("OS","STATOS")])
dim(clin)
head(clin)
clin$Cancer <- c(rep("Lung Primary Tumor",85),rep("Lung NAC Tumor",77),rep("Lung squamous cell carcinoma",238),
                 rep("Primary EAC",102),rep("Neoadjuvant treated EAC",144),rep("Chromophobe RCC",100),
                 rep("Clear cell RCC",476),rep("Papillary RCC",101),rep("Adrenocortical Carcinoma",72),
                 rep("Pancreatic cancer",107),rep("Gastric cancer",246),rep("HER2-Gastric cancer",82))
dim(clin)
head(clin)


table7 <- table6[match(rownames(clin),rownames(table6)),]
table7[1:5,1:5]
dim(table7)
#write.csv(table7,"Metabolite profiles of 105 nucleotides in 1830 patients1.csv")
#write.csv(t(table7),"Metabolite profiles of 105 nucleotides in 1830 patients2.csv")
table7[is.na(table7)] <- 0
table8 <- ifelse(table7>0,1,0)
table8[1:5,1:5]

table9 <- cbind(clin,table8)
table9[1:5,1:5]
dim(table9)

table10 <- table9
#write.csv(table9,"Metabolite presence or absence of 105 nucleotides in 1830 patients.csv")
#write.csv(t(table9),"Metabolite presence or absence of 105 nucleotides in 1830 patients1.csv")

#table10 <- read.csv("Metabolite presence or absence of 105 nucleotides in 1830 patients.csv")
#rownames(table10) <- table10[,1]
#table10 <- table10[,-1]
table10[1:5,1:5]
dim(table10)

table10$STATOS1 <- ifelse(table10$OS < median(table10$OS),1,0) #1=short/poor survival or high risk,0=long/good survival or low risk  
table10$STATOS1
dim(table10)

table11 <- as.data.frame(t(table10[,-c(1:3,109)]))
dim(table11)

annotation_col <- data.frame(Type = table10$Cancer,STATOS=table10$STATOS,Risk=table10$STATOS1,
                             Cluster=sum55$Cluster)
rownames(annotation_col) <- rownames(table10)
head(annotation_col)
dim(annotation_col)
annotation_col$STATOS <- as.factor(annotation_col$STATOS)
annotation_col$Risk <- as.factor(annotation_col$Risk)
annotation_col$Cluster <- as.factor(annotation_col$Cluster)
str(annotation_col)

anncol = list(Risk=c("0"="gray","1"="purple"),
              STATOS=c("0"="gray","1"="red"),
              Cluster=c("1"="red","2"="green","3"="blue"))
anncol

library(pheatmap)
clustering.method <- c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty","median","centroid")
clustering.distance.rows <- c("euclidean", "maximum", "manhattan", "canberra",  "minkowski", "correlation")
clustering.distance.cols <- c("euclidean", "maximum", "manhattan", "canberra",  "minkowski", "correlation")

pheatmap1 <- pheatmap(table11,cluster_cols =T,cluster_rows =T,annotation_col = annotation_col,annotation_colors=anncol,
                      clustering_method = "ward.D",
                      clustering_distance_rows = "euclidean",
                      clustering_distance_cols = "euclidean",
                      cutree_cols=3,
                      fontsize = 6,
                      show_colnames = F,
                      show_rownames = T)

cluster.col <- pheatmap1$tree_col

#对聚类树进行分簇
cut <- cutree(cluster.col,3)

cut1 <- names(cut)[cut==1]
cut2 <- names(cut)[cut==2]
cut3 <- names(cut)[cut==3]


### pan-cancer subgroup constitution
sub1 <- table9[cut1,]
sub2 <- table9[cut2,]
sub3 <- table9[cut3,]
dim(sub1)
dim(sub2)
dim(sub3)
head(sub3)

sub11 <- by(sub1,sub1$Cancer,count)
sub11 <- data.frame(name=names(sub11),count=as.numeric(sub11))
dim(sub11)
head(sub11)
sub22 <- by(sub2,sub2$Cancer,count)
sub22 <- data.frame(name=names(sub22),count=as.numeric(sub22))
dim(sub22)
head(sub22)
sub33 <- by(sub3,sub3$Cancer,count)
sub33 <- data.frame(name=names(sub33),count=as.numeric(sub33))
dim(sub33)
head(sub33)

sub <- by(table9,table9$Cancer,count)
sub <- data.frame(name=names(sub),count=as.numeric(sub))
colnames(sub)[2] <- "count1"
sub

sub111 <- merge(sub11,sub)
sub111$count2 <- sub111$count/sub111$count1
sub111$count3 <- sub111$count2/sum(sub111$count2)
sub111

sub222 <- merge(sub22,sub)
sub222$count2 <- sub222$count/sub222$count1
sub222$count3 <- sub222$count2/sum(sub222$count2)
sub222

sub333 <- merge(sub33,sub)
sub333$count2 <- sub333$count/sub333$count1
sub333$count3 <- sub333$count2/sum(sub333$count2)
sub333

sub.2 <- sub111
sub.2 <- sub222
sub.2 <- sub333
label_value <- paste(round(sub.2$count3/sum(sub.2$count3) * 100, 1), '%', sep = '')
label_value
ggplot(data=sub.2,aes(x='Content',y=count3,fill=name))+geom_bar(stat = 'identity',position = 'stack')+
  coord_polar(theta = 'y')+
  labs(x='',y='',title = '')+
  theme(axis.text = element_blank(),axis.ticks = element_blank(),panel.grid = element_blank(),panel.border = element_blank())+
  geom_text(aes(y = count3/2 + c(0, cumsum(count3)[-length(count3)]), x =1.2, label = label_value)) 



############# survival analysis on patients group from different clusters
dim(table10)
table10[1:5,1:5]
sum55 <- table10
sum55$Cluster <- rep(1,1830)
sum55$Cluster[match(cut1,rownames(sum55))] <- rep(3,length(match(cut1,rownames(sum55))))
sum55$Cluster[match(cut2,rownames(sum55))] <- rep(2,length(match(cut2,rownames(sum55))))
sum55$Cluster[match(cut3,rownames(sum55))] <- rep(1,length(match(cut3,rownames(sum55))))
head(sum55)


# fit survival curve
fit <- survfit(Surv(OS,STATOS)~Cluster,data =sum55)
fit
survdiff(Surv(OS,STATOS)~Cluster,data =sum55) #采用log-rank 检验分析生存率差异

##绘制生存曲线##
ggsurvplot(fit, data = sum55, palette = c("red","green","blue"),conf.int = F, # 增加置信区间
           risk.table = TRUE, pval=T) #不同时间点风险人数表，增加P值

res.cox <- coxph(Surv(OS,STATOS)~Cluster,data =sum55)
summary(res.cox)
cox.zph(res.cox)


### survival analysis on any two subgroups
sum66 <- sum55[c(which(sum55$Cluster==1),which(sum55$Cluster==2)),] #*
sum66 <- sum55[c(which(sum55$Cluster==1),which(sum55$Cluster==3)),] #***
sum66 <- sum55[c(which(sum55$Cluster==2),which(sum55$Cluster==3)),] #*


# fit survival curve
fit <- survfit(Surv(OS,STATOS)~Cluster,data =sum66)
fit
survdiff(Surv(OS,STATOS)~Cluster,data =sum66) #采用log-rank 检验分析生存率差异

##绘制生存曲线##
ggsurvplot(fit, data = sum66,palette = c("green","blue"), conf.int = F, # 增加置信区间
           risk.table = TRUE, pval=T) #不同时间点风险人数表，增加P值

res.cox <- coxph(Surv(OS,STATOS)~Cluster,data =sum66)
summary(res.cox)
cox.zph(res.cox)




####### Pan-cancer cluster analysis on all patients based on 29 common nucleotides and survival analysis #######################

sum.29.com1 <- data.frame(name=sum.29.com)
sum.29.com1

dim(univariate.sum)
head(univariate.sum)

# significant 29 nucleotides in all cancers
univariate.sum.29 <- merge(sum.29.com1,univariate.sum)
dim(univariate.sum.29)
head(univariate.sum.29)


setwd("C:/MyRdata8/Summary")
sum2 <- read.csv("annotated and detectable 105 nucleotides in all cancers.csv")
sum3 <- sum2[,-1]
head(sum3)
dim(sum3)
length(unique(sum3$name))
unique(sum3$tumor)


tumor <- c("Lung Primary Tumor","Lung NAC Tumor","Lung squamous cell carcinoma",
           "Primary EAC","Neoadjuvant treated EAC","Chromophobe RCC",
           "Clear cell RCC","Papillary RCC","Adrenocortical Carcinoma",
           "Pancreatic cancer","Gastric cancer","HER2-Gastric cancer")
tumor


# consider univariate significance and intensity
sum10 <- data.frame()
for (i in 1:12) {
  sum4 <- sum3[which(sum3$tumor==tumor[i]),]
  univariate.sum.29.1 <- univariate.sum.29[which(univariate.sum.29$Tumor==tumor[i]),c(1,3,24)]
  sum5 <- data.frame(name=setdiff(sum.29.com,univariate.sum.29.1$name))
  sum6 <- merge(sum5,sum4)
  sum6$name.adduct <- paste(sum6$adduct,sum6$name,sep = '.')
  sum6 <- sum6[order(sum6$ppm,decreasing = F),]
  sum7 <- sum6[which(duplicated(sum6$name.adduct)==F),]
  sum7 <- sum7[order(sum7$intensity,decreasing = T),]
  sum8 <- sum7[which(duplicated(sum7$name)==F),]
  sum9 <- data.frame(name=c(sum8$name,univariate.sum.29.1$name),
                     mass=c(sum8$mass,univariate.sum.29.1$mz.value),
                     tumor=c(sum8$tumor,univariate.sum.29.1$Tumor))
  sum10 <- rbind(sum10,sum9)
}
dim(sum10)
head(sum10)
#write.csv(sum10,"annotation on 29 nucleotides in all tumor considering univariate significance and intensity.csv")

# No run: only consider intensity for annotation
sum10 <- data.frame()
for (i in 1:12) {
  sum6 <- sum3[which(sum3$tumor==tumor[i]),]
  sum6$name.adduct <- paste(sum6$adduct,sum6$name,sep = '.')
  sum6 <- sum6[order(sum6$ppm,decreasing = F),]
  sum7 <- sum6[which(duplicated(sum6$name.adduct)==F),]
  sum7 <- sum7[order(sum7$intensity,decreasing = T),]
  sum8 <- sum7[which(duplicated(sum7$name)==F),]
  sum9 <- sum8[match(sum.28.com,sum8$name),c(9,2,4)]
  sum10 <- rbind(sum10,sum9)
}
dim(sum10)
head(sum10)
#write.csv(sum10,"annotation on 28 nucleotides in all tumor only considering intensity.csv")


############################## metabolic profiles of 29 common nucleotides by surv_cutpoint for each tumor :No run
#surv_cutpoint minprop = 0.01
sum10.n85 <- res.cat.lung.primary.tumor.n85[,sum10[sum10$tumor==tumor[1],2]]
colnames(sum10.n85) <- sum10[sum10$tumor==tumor[1],1]
dim(sum10.n85)

sum10.n77 <- res.cat.lung.NAC.tumor.n77[,sum10[sum10$tumor==tumor[2],2]]
colnames(sum10.n77) <- sum10[sum10$tumor==tumor[2],1]

sum10.n238 <- res.cat.lunge.plattenepithel.n238[,sum10[sum10$tumor==tumor[3],2]]
colnames(sum10.n238) <- sum10[sum10$tumor==tumor[3],1]

sum10.n102 <- res.cat.all.primary.EAC.n102[,sum10[sum10$tumor==tumor[4],2]]
colnames(sum10.n102) <- sum10[sum10$tumor==tumor[4],1]

sum10.n144 <- res.cat.Neoadjuvant.treated.EACs.n144[,sum10[sum10$tumor==tumor[5],2]]
colnames(sum10.n144) <- sum10[sum10$tumor==tumor[5],1]


dim(res.cat.Chromophobe.RCC.n100)
rcc <- sum10[sum10$tumor==tumor[6],]
com <- data.frame(mass=intersect(unique(rcc$mass),colnames(res.cat.Chromophobe.RCC.n100)))
com1 <- merge(com,rcc)
sum10.n100 <- res.cat.Chromophobe.RCC.n100[,com1$mass]
colnames(sum10.n100) <- com1$name
dim(sum10.n100)


dim(res.cat.Clear.cell.RCC.n476)
rcc <- sum10[sum10$tumor==tumor[7],]
com <- data.frame(mass=intersect(unique(rcc$mass),colnames(res.cat.Clear.cell.RCC.n476)))
com1 <- merge(com,rcc)
sum10.n476 <- res.cat.Clear.cell.RCC.n476[,com1$mass]
colnames(sum10.n476) <- com1$name
dim(sum10.n476)


dim(res.cat.Papillary.RCC.n101)
rcc <- sum10[sum10$tumor==tumor[8],]
com <- data.frame(mass=intersect(unique(rcc$mass),colnames(res.cat.Papillary.RCC.n101)))
com1 <- merge(com,rcc)
sum10.n101 <- res.cat.Papillary.RCC.n101[,com1$mass]
colnames(sum10.n101) <- com1$name
dim(sum10.n101)


sum10.n72 <- res.cat.ACC.n72[,sum10[sum10$tumor==tumor[9],2]]
colnames(sum10.n72) <- sum10[sum10$tumor==tumor[9],1]


dim(res.cat.Pancreas.tumor.n107)
pan <- sum10[sum10$tumor==tumor[10],]
com <- data.frame(mass=intersect(unique(pan$mass),colnames(res.cat.Pancreas.tumor.n107)))
com1 <- merge(com,pan)
sum10.n107 <- res.cat.Pancreas.tumor.n107[,com1$mass]
colnames(sum10.n107) <- com1$name
dim(sum10.n107)


sum10.n246 <- res.cat.Gastric.Cancer.n246[,sum10[sum10$tumor==tumor[11],2]]
colnames(sum10.n246) <- sum10[sum10$tumor==tumor[11],1]


dim(res.cat.HER2.n82)
her <- sum10[sum10$tumor==tumor[12],]
com <- data.frame(mass=intersect(unique(her$mass),colnames(res.cat.HER2.n82)))
com1 <- merge(com,her)
sum10.n82 <- res.cat.HER2.n82[,com1$mass]
colnames(sum10.n82) <- com1$name
dim(sum10.n82)


com2 <- intersect(intersect(intersect(intersect(colnames(sum10.n82),colnames(sum10.n107)),colnames(sum10.n101)),colnames(sum10.n476)),colnames(sum10.n100))
length(com2)

sum10.n85 <- sum10.n85[,com2]
sum10.n77 <- sum10.n77[,com2]
sum10.n238 <- sum10.n238[,com2]
sum10.n102 <- sum10.n102[,com2]
sum10.n144 <- sum10.n144[,com2]
sum10.n100 <- sum10.n100[,com2]
sum10.n476 <- sum10.n476[,com2]
sum10.n101 <- sum10.n101[,com2]
sum10.n72 <- sum10.n72[,com2]
sum10.n107 <- sum10.n107[,com2]
sum10.n246 <- sum10.n246[,com2]
sum10.n82 <- sum10.n82[,com2]

sum11 <- rbind(sum10.n85,sum10.n77,sum10.n238,sum10.n102,
               sum10.n144,sum10.n100,sum10.n476,sum10.n101,
               sum10.n72,sum10.n107,sum10.n246,sum10.n82)
dim(sum11)
head(sum11)
head(nucleotide)
colnames(sum11) <- nucleotide[match(colnames(sum11),nucleotide$name),5]
profile.19.com <- sum11

#write.csv(sum11,"Profile data of 19 nucleotides in all patients.csv")
sum11 <- read.csv("Profile data of 19 nucleotides in all patients.csv")
rownames(sum11) <- sum11[,1]
sum11 <- sum11[,-1]
dim(sum11)
head(sum11)
head(nucleotide)
colnames(sum11) <- nucleotide[match(colnames(sum11),nucleotide$compound_id),4]
dim(sum11)
head(sum11)
sum13 <- t(sum11)
dim(sum13)



############################## metabolic profiles of 28 common nucleotides by scale for each tumor

sum11.1 <- nucleotide.profile.lung.primary.tumor.n85[,sum10[sum10$tumor==tumor[1],2]]
colnames(sum11.1) <- sum10[sum10$tumor==tumor[1],1]
sum11.1 <- sum11.1[,match(sum.29.com,colnames(sum11.1))]
sum11.1 <- scale(sum11.1,center = T,scale = T)
head(sum11.1)

sum11.2 <- nucleotide.profile.lung.NAC.tumor.n77[,sum10[sum10$tumor==tumor[2],2]]
colnames(sum11.2) <- sum10[sum10$tumor==tumor[2],1]
sum11.2 <- sum11.2[,match(sum.29.com,colnames(sum11.2))]
sum11.2 <- scale(sum11.2,center = T,scale = T)
head(sum11.2)

sum11.3 <- nucleotide.profile.lunge.plattenepithel.n238[,sum10[sum10$tumor==tumor[3],2]]
colnames(sum11.3) <- sum10[sum10$tumor==tumor[3],1]
sum11.3 <- sum11.3[,match(sum.29.com,colnames(sum11.3))]
sum11.3 <- scale(sum11.3,center = T,scale = T)
head(sum11.3)

sum11.4 <- nucleotide.profile.all.primary.EAC.n102[,sum10[sum10$tumor==tumor[4],2]]
colnames(sum11.4) <- sum10[sum10$tumor==tumor[4],1]
sum11.4 <- sum11.4[,match(sum.29.com,colnames(sum11.4))]
sum11.4 <- scale(sum11.4,center = T,scale = T)
head(sum11.4)

sum11.5 <- nucleotide.profile.Neoadjuvant.treated.EACs.n144[,sum10[sum10$tumor==tumor[5],2]]
colnames(sum11.5) <- sum10[sum10$tumor==tumor[5],1]
sum11.5 <- sum11.5[,match(sum.29.com,colnames(sum11.5))]
sum11.5 <- scale(sum11.5,center = T,scale = T)
head(sum11.5)

sum11.6 <- nucleotide.profile.Chromophobe.RCC.n100[,sum10[sum10$tumor==tumor[6],2]]
colnames(sum11.6) <- sum10[sum10$tumor==tumor[6],1]
sum11.6 <- sum11.6[,match(sum.29.com,colnames(sum11.6))]
sum11.6 <- scale(sum11.6,center = T,scale = T)
head(sum11.6)

sum11.7 <- nucleotide.profile.Clear.cell.RCC.n476[,sum10[sum10$tumor==tumor[7],2]]
colnames(sum11.7) <- sum10[sum10$tumor==tumor[7],1]
sum11.7 <- sum11.7[,match(sum.29.com,colnames(sum11.7))]
sum11.7 <- scale(sum11.7,center = T,scale = T)
head(sum11.7)

sum11.8 <- nucleotide.profile.Papillary.RCC.n101[,sum10[sum10$tumor==tumor[8],2]]
colnames(sum11.8) <- sum10[sum10$tumor==tumor[8],1]
sum11.8 <- sum11.8[,match(sum.29.com,colnames(sum11.8))]
sum11.8 <- scale(sum11.8,center = T,scale = T)
head(sum11.8)

sum11.9 <- nucleotide.profile.ACC.n72[,sum10[sum10$tumor==tumor[9],2]]
colnames(sum11.9) <- sum10[sum10$tumor==tumor[9],1]
sum11.9 <- sum11.9[,match(sum.29.com,colnames(sum11.9))]
sum11.9 <- scale(sum11.9,center = T,scale = T)
head(sum11.9)

sum11.10 <- nucleotide.profile.Pancreas.tumor.n107[,sum10[sum10$tumor==tumor[10],2]]
colnames(sum11.10) <- sum10[sum10$tumor==tumor[10],1]
sum11.10 <- sum11.10[,match(sum.29.com,colnames(sum11.10))]
sum11.10 <- scale(sum11.10,center = T,scale = T)
head(sum11.10)

sum11.11 <- nucleotide.profile.Gastric.Cancer.n246[,sum10[sum10$tumor==tumor[11],2]]
colnames(sum11.11) <- sum10[sum10$tumor==tumor[11],1]
sum11.11 <- sum11.11[,match(sum.29.com,colnames(sum11.11))]
sum11.11 <- scale(sum11.11,center = T,scale = T)
head(sum11.11)

sum11.12 <- nucleotide.profile.HER2.n82[,sum10[sum10$tumor==tumor[12],2]]
colnames(sum11.12) <- sum10[sum10$tumor==tumor[12],1]
sum11.12 <- sum11.12[,match(sum.29.com,colnames(sum11.12))]
sum11.12 <- scale(sum11.12,center = T,scale = T)
head(sum11.12)


sum11 <- rbind(sum11.1,sum11.2,sum11.3,sum11.4,sum11.5,sum11.6,sum11.7,sum11.8,
               sum11.9,sum11.10,sum11.11,sum11.12)
dim(sum11)
head(sum11)
profile.29.com <- sum11

#No run: fill NaN values
sum11[is.nan(sum11)==T] <- (-2)
which(is.nan(sum11[,17])==T)

# remove nucleotide with NaN values
sum11 <- sum11[,-c(2,3,7,8,13,14,15,16,17)]
dim(sum11)
head(sum11)


#### cutoff defined in all patients by surv_cutpoint after scale in each cancer type for each nucleotide
head(nucleotide)
colnames(sum11) <- nucleotide[match(colnames(sum11),nucleotide$name),6]
head(sum11)
dim(sum11)

#No run: remove duplicated nucleotides with same profile data
index1 <- which(duplicated(sum11[1,])==T)
for (i in 2:1830) {
  index <- which(duplicated(sum11[i,])==T)
  index1 <- intersect(index1,index)
}
length(index1)
sum11 <- sum11[,-index1]
dim(sum11)


# cutoff defined in all patients by surv_cutpoint
sum12 <- cbind(clin,sum11)
sum12 <- sum12[,-3]
dim(sum12)
head(sum12)

res.cut <- surv_cutpoint(sum12, time = "OS", event = "STATOS",
                         minprop = 0.1,
                         variables = colnames(sum12)[3:ncol(sum12)])
res.cat <- surv_categorize(res.cut,labels = c(0, 1))
head(res.cat)
sum12 <- res.cat

sum13 <- as.data.frame(t(sum12[,-c(1,2)]))
rownames(sum13) <- nucleotide[match(rownames(sum13),nucleotide$compound_id),4]
dim(sum13)
sum13[1:5,1:5]



#################################################### heatmap
annotation_col <- data.frame(Type = clin$Cancer)
rownames(annotation_col) <- colnames(sum13)
head(annotation_col)
dim(annotation_col)

sum13 <- distinct(as.data.frame(sum13))
dim(sum13)
pheatmap1 <- pheatmap(sum13,cluster_cols =T,cluster_rows =T,annotation_col = annotation_col,
                      clustering_method = "ward.D2",
                      clustering_distance_rows = "euclidean",
                      clustering_distance_cols = "euclidean",
                      cutree_cols=4,
                      show_colnames = F,
                      show_rownames = T)

cluster.col <- pheatmap1$tree_col

#对聚类树进行分簇
cut <- cutree(cluster.col,4)

cut1 <- names(cut)[cut==1]
cut2 <- names(cut)[cut==2]
cut3 <- names(cut)[cut==3]
cut4 <- names(cut)[cut==4]

############# survival analysis on patients group from different clusters
dim(clin)
head(clin)
sum55 <- clin
sum55$Cluster <- rep(1,1830)
sum55$Cluster[match(cut1,rownames(sum55))] <- rep(1,length(match(cut1,rownames(sum55))))
sum55$Cluster[match(cut2,rownames(sum55))] <- rep(2,length(match(cut2,rownames(sum55))))
sum55$Cluster[match(cut3,rownames(sum55))] <- rep(3,length(match(cut3,rownames(sum55))))
sum55$Cluster[match(cut4,rownames(sum55))] <- rep(4,length(match(cut4,rownames(sum55))))
head(sum55)


# fit survival curve
fit <- survfit(Surv(OS,STATOS)~Cluster,data =sum55)
fit
survdiff(Surv(OS,STATOS)~Cluster,data =sum55) #采用log-rank 检验分析生存率差异

##绘制生存曲线##
ggsurvplot(fit, data = sum55, conf.int = F, # 增加置信区间
           risk.table = TRUE, pval=T) #不同时间点风险人数表，增加P值

res.cox <- coxph(Surv(OS,STATOS)~Cluster,data =sum55)
summary(res.cox)
cox.zph(res.cox)


### survival analysis on any two subgroups
sum66 <- sum55[c(which(sum55$Cluster==1),which(sum55$Cluster==2)),] #***
sum66 <- sum55[c(which(sum55$Cluster==1),which(sum55$Cluster==3)),] #***
sum66 <- sum55[c(which(sum55$Cluster==2),which(sum55$Cluster==3)),] #***
sum66 <- sum55[c(which(sum55$Cluster==1),which(sum55$Cluster==4)),] #***
sum66 <- sum55[c(which(sum55$Cluster==2),which(sum55$Cluster==4)),] #***
sum66 <- sum55[c(which(sum55$Cluster==3),which(sum55$Cluster==4)),] #***

# fit survival curve
fit <- survfit(Surv(OS,STATOS)~Cluster,data =sum66)
fit
survdiff(Surv(OS,STATOS)~Cluster,data =sum66) #采用log-rank 检验分析生存率差异

##绘制生存曲线##
ggsurvplot(fit, data = sum66, conf.int = F, # 增加置信区间
           risk.table = TRUE, pval=T) #不同时间点风险人数表，增加P值

res.cox <- coxph(Surv(OS,STATOS)~Cluster,data =sum66)
summary(res.cox)
cox.zph(res.cox)



########## pan-cancer subgroup constitution
sub1 <- clin[cut1,]
sub2 <- clin[cut2,]
sub3 <- clin[cut3,]
sub4 <- clin[cut4,]
dim(sub1)
dim(sub2)
dim(sub3)
head(sub3)

sub11 <- by(sub1,sub1$Cancer,count)
sub11 <- data.frame(name=names(sub11),count=as.numeric(sub11))
dim(sub11)
head(sub11)
sub22 <- by(sub2,sub2$Cancer,count)
sub22 <- data.frame(name=names(sub22),count=as.numeric(sub22))
dim(sub22)
head(sub22)
sub33 <- by(sub3,sub3$Cancer,count)
sub33 <- data.frame(name=names(sub33),count=as.numeric(sub33))
dim(sub33)
head(sub33)
sub44 <- by(sub4,sub4$Cancer,count)
sub44 <- data.frame(name=names(sub44),count=as.numeric(sub44))
dim(sub44)
head(sub44)

sub <- by(clin,clin$Cancer,count)
sub <- data.frame(name=names(sub),count=as.numeric(sub))
colnames(sub)[2] <- "count1"
sub

sub111 <- merge(sub11,sub)
sub111$count2 <- sub111$count/sub111$count1
sub111$count3 <- sub111$count2/sum(sub111$count2)
sub111

sub222 <- merge(sub22,sub)
sub222$count2 <- sub222$count/sub222$count1
sub222$count3 <- sub222$count2/sum(sub222$count2)
sub222

sub333 <- merge(sub33,sub)
sub333$count2 <- sub333$count/sub333$count1
sub333$count3 <- sub333$count2/sum(sub333$count2)
sub333

sub444 <- merge(sub44,sub)
sub444$count2 <- sub444$count/sub444$count1
sub444$count3 <- sub444$count2/sum(sub444$count2)
sub444

sub.2 <- sub111
sub.2 <- sub222
sub.2 <- sub333
sub.2 <- sub444
label_value <- paste(round(sub.2$count3/sum(sub.2$count3) * 100, 1), '%', sep = '')
label_value
ggplot(data=sub.2,aes(x='Content',y=count3,fill=name))+geom_bar(stat = 'identity',position = 'stack')+
  coord_polar(theta = 'y')+
  labs(x='',y='',title = '')+
  theme(axis.text = element_blank(),axis.ticks = element_blank(),panel.grid = element_blank(),panel.border = element_blank())+
  geom_text(aes(y = count3/2 + c(0, cumsum(count3)[-length(count3)]), x =1.2, label = label_value)) 






############################################### consensus clustering in all patients for each nucleotide

setwd("C:/MyRdata8/Summary")
library(survival)
library(survminer)
library(ConsensusClusterPlus)
library(dplyr)

dim(profile.28.com)
head(profile.28.com)
sum1 <- profile.28.com
sum1[is.nan(sum1)==T] <- (-2)
dim(sum1)
head(sum1)
head(nucleotide)
colnames(sum1) <- nucleotide[match(colnames(sum1),nucleotide$name),5]

sum6 <- cbind(clin[,1:2],sum1)
head(sum6)
res.cut <- surv_cutpoint(sum6, time = "OS", event = "STATOS",
                         minprop = 0.1,
                         variables = colnames(sum6)[3:ncol(sum6)])
res.cat <- surv_categorize(res.cut,labels = c(0, 1))
head(res.cat)
sum7 <- res.cat
head(sum7)
dim(sum7)

sum8 <- as.matrix(distinct(as.data.frame(t(sum7[,-c(1,2)])))) #17
sum8 <- as.matrix(t(sum7[,-c(1,2)])) #28
dim(sum8)

annotation_col <- data.frame(Type = clin$Cancer)
rownames(annotation_col) <- colnames(sum8)
head(annotation_col)
dim(annotation_col)


title = 'C:/MyRdata8/Summary/ConsensusClusterPlus'
results <- ConsensusClusterPlus(sum8,maxK = 6,reps = 50,pItem = 0.8,pFeature = 1,title = title,
                                clusterAlg = 'hc',distance = 'euclidean',plot = 'png')


results[[2]][['consensusClass']][1:50]

icl <- calcICL(results,title=title,plot = 'png')


############# survival analysis on patients group from different clusters
dim(sum7)
head(sum7)
sum55 <- sum7
sum55$Cluster <- results[[2]][['consensusClass']]
head(sum55)
dim(sum55)

# fit survival curve
fit <- survfit(Surv(OS,STATOS)~Cluster,data =sum55)
fit
survdiff(Surv(OS,STATOS)~Cluster,data =sum55) #采用log-rank 检验分析生存率差异

##绘制生存曲线##
ggsurvplot(fit, data = sum55, conf.int = F, # 增加置信区间
           risk.table = TRUE, pval=T) #不同时间点风险人数表，增加P值


### survival analysis on any two subgroups
sum66 <- sum55[c(which(sum55$Cluster==1),which(sum55$Cluster==2)),] #***
sum66 <- sum55[c(which(sum55$Cluster==1),which(sum55$Cluster==3)),] #***
sum66 <- sum55[c(which(sum55$Cluster==2),which(sum55$Cluster==3)),] #***



# fit survival curve
fit <- survfit(Surv(OS,STATOS)~Cluster,data =sum66)
fit
survdiff(Surv(OS,STATOS)~Cluster,data =sum66) #采用log-rank 检验分析生存率差异

##绘制生存曲线##
ggsurvplot(fit, data = sum66, conf.int = F, # 增加置信区间
           risk.table = TRUE, pval=T) #不同时间点风险人数表，增加P值




############### pan-cancer subgroup constitution
cut1 <- which(results[[2]][['consensusClass']]==1)
cut2 <- which(results[[2]][['consensusClass']]==2)
sub1 <- clin[cut1,]
sub2 <- clin[cut2,]
sub3 <- clin[cut3,]
dim(sub1)
dim(sub2)
dim(sub3)
head(sub1)

sub11 <- by(sub1,sub1$Cancer,count)
sub11 <- data.frame(name=names(sub11),count=as.numeric(sub11))
dim(sub11)
head(sub11)
sub22 <- by(sub2,sub2$Cancer,count)
sub22 <- data.frame(name=names(sub22),count=as.numeric(sub22))
dim(sub22)
head(sub22)
sub33 <- by(sub3,sub3$Cancer,count)
sub33 <- data.frame(name=names(sub33),count=as.numeric(sub33))
dim(sub33)
head(sub33)

sub <- by(clin,clin$Cancer,count)
sub <- data.frame(name=names(sub),count=as.numeric(sub))
colnames(sub)[2] <- "count1"
sub

sub111 <- merge(sub11,sub)
sub111$count2 <- sub111$count/sub111$count1
sub111$count3 <- sub111$count2/sum(sub111$count2)
sub111

sub222 <- merge(sub22,sub)
sub222$count2 <- sub222$count/sub222$count1
sub222$count3 <- sub222$count2/sum(sub222$count2)
sub222

sub333 <- merge(sub33,sub)
sub333$count2 <- sub333$count/sub333$count1
sub333$count3 <- sub333$count2/sum(sub333$count2)
sub333

sub.2 <- sub111
sub.2 <- sub222
label_value <- paste(round(sub.2$count3/sum(sub.2$count3) * 100, 1), '%', sep = '')
label_value
ggplot(data=sub.2,aes(x='Content',y=count3,fill=name))+geom_bar(stat = 'identity',position = 'stack')+
  coord_polar(theta = 'y')+
  labs(x='',y='',title = '')+
  theme(axis.text = element_blank(),axis.ticks = element_blank(),panel.grid = element_blank(),panel.border = element_blank())+
  geom_text(aes(y = count3/2 + c(0, cumsum(count3)[-length(count3)]), x =1.2, label = label_value)) 




############################ Pan-cancer prognosis analysis on all nucleotides and all patients ################

####################################### method1:combine all data after scale each tumor
dim(profile.29.com)
head(profile.29.com)
profile.29.com[is.nan(profile.29.com)==T] <- (-2)

dim(clin)
head(clin)
head(nucleotide)
profile.29.com1 <- cbind(clin[,c(1,2)],profile.29.com)
colnames(profile.29.com1)[3:31] <- nucleotide[match(colnames(profile.29.com1)[3:31],nucleotide$name),6]
dim(profile.29.com1)
head(profile.29.com1)

setwd("C:/MyRdata8/Summary")
#write.csv(profile.29.com1,"profile.29.com1.csv")
profile.29.com1 <- read.csv("profile.29.com1.csv")
rownames(profile.29.com1) <- profile.29.com1[,1]
profile.29.com1 <- profile.29.com1[,-c(1,4)]
dim(profile.29.com1)
head(profile.29.com1)

res.cut <- surv_cutpoint(profile.29.com1, time = "OS", event = "STATOS",
                         minprop = 0.1,
                         variables = colnames(profile.29.com1)[3:ncol(profile.29.com1)])
res.cat <- surv_categorize(res.cut,labels = c(0, 1))
head(res.cat)
dim(res.cat)

### survival curve
# fit survival curve
head(nucleotide)
fit <- survfit(Surv(OS,STATOS)~HMDB0000133,data =res.cat)
fit
survdiff(Surv(OS,STATOS)~HMDB0000133,data =res.cat) #采用log-rank 检验分析生存率差异

##绘制生存曲线##
ggsurvplot(fit, data = res.cat, conf.int = F, # 增加置信区间
           risk.table = TRUE, pval=T,
           palette = c("black","red")) #不同时间点风险人数表，增加P值

res.cox <- coxph(Surv(OS,STATOS)~HMDB0000133,data =res.cat)
summary(res.cox)
cox.zph(res.cox)



### univariate cox analysis to calculate p value, HR
result <- data.frame()
for (i in 3:ncol(res.cat)) {
  res.cox <- coxph(Surv(OS,STATOS)~res.cat[,i],data =res.cat)
  summ <- summary(res.cox)
  hypo <- cox.zph(res.cox)
  hypo1 <- hypo$table[-nrow(hypo$table),]
  result1 <- data.frame(HMDB=colnames(res.cat)[i],
                        HR = signif(summ$coef[2], digits=3),
                        HR.confint.lower = signif(summ$conf.int[,"lower .95"], 3),
                        HR.confint.upper = signif(summ$conf.int[,"upper .95"],3),
                        logtest.p=summ$logtest[3],
                        waldtest.p=summ$waldtest[3],
                        sctest.p=summ$sctest[3],
                        hypo.p=hypo1[3])
  result <- rbind(result,result1)
}
head(result)
dim(result)

# screen by sctest p value
result2 <- result[result$sctest.p<0.05,]
result2 <- result2[result2$hypo.p>0.05,]
dim(result2)
head(result2)

dim(nucleotide)
head(nucleotide)

result2$Nucleotide <- nucleotide[match(result2$HMDB,nucleotide$compound_id),4]
dim(result2)
head(result2)

result3 <- result
result3$Nucleotide <- nucleotide[match(result3$HMDB,nucleotide$compound_id),4]
dim(result3)
head(result3)

# forest plot
library(forestplot)
result6 <- result2
result6 <- result3
dim(result6)
head(result6)
which(duplicated(result6$HMDB)==T)

result6$HR1 <- paste(result6$HR,"(",result6$HR.confint.lower,"-",result6$HR.confint.upper,")",sep = '')
result6$pvalue <- ifelse(result6$sctest.p<0.001,"<0.001",ifelse(result6$sctest.p<0.01,"<0.01",ifelse(result6$sctest.p<0.01,"<0.05","NS")))
dim(result6)
head(result6)

labeltext <- as.matrix(rbind(c("Nucleotide","HR","pvalue"),result6[,c(9,10,11)]))
labeltext

forestplot(labeltext, 
           mean = c(NA,result6$HR),
           lower = c(NA,result6$HR.confint.lower),
           upper = c(NA,result6$HR.confint.upper),
           zero = 1,lwd.zero = 2,
           clip = c(0,10),
           graph.pos = 2,
           boxsize = 0.5,
           align = "l",
           xlog = FALSE,
           col=fpColors(line = "black", zero = "black",
                        box="black"),
           lty.ci = 7,
           lwd.ci = 1,
           ci.vertices.height = 0.15,
           lwd.xaxis=3,
           xlab="Hazard Ratio",
           txt_gp = fpTxtGp(ticks = gpar(cex = 0.9), xlab = gpar(cex = 0.8), cex = 0.8),
           lineheight = unit(8,"mm")) 



### multivariate analysis
head(res.cat)
dim(res.cat)
head(nucleotide)

res.cat1 <- res.cat[,c("OS","STATOS",result2$HMDB)]
dim(res.cat1)
#colnames(res.cat1)[3:ncol(res.cat1)] <- nucleotide[match(colnames(res.cat1)[3:ncol(res.cat1)],nucleotide$compound_id),4]
head(res.cat1)
dim(res.cat1)

res.cox1 <- coxph(Surv(OS, STATOS) ~ .,data =res.cat1)
res.cox2 <- step(res.cox1)
ggforest(model=res.cox2,data=res.cat1)

#等比例风险假定
hypo <- cox.zph(res.cox2) 
hypo1 <- hypo$table[-nrow(hypo$table),]
dim(hypo1)
hypo2 <- hypo1[hypo1[,3]>0.05,]
dim(hypo2)
hypo2


### forest plot
summ <- as.data.frame(cbind(summary(res.cox2)$coefficients[,c(2,5)],summary(res.cox2)$conf.int[,c(3,4)]))
summ$HMDB <- rownames(summ)
summ$name <- nucleotide[match(summ$HMDB,nucleotide$compound_id),5]
head(summ)
dim(summ)


summ1 <- summ[summ$`Pr(>|z|)`<0.05,]
dim(summ1)
head(summ1)
#write.csv(summ1,"11 independent pan-cancer nucleotides.csv")

nucleotide13 <- c(summ1$HMDB,"HMDB0000542","HMDB0000288","HMDB0001271","HMDB0000058","HMDB0001314","HMDB0000061","HMDB0000960")
#n11 <- c(summ1$name,"8-Hydroxyadenine","Pseudouridylic acid","dGDP")
#write.csv(n11,"11 pan-cancer nucleotides.csv")

library(forestplot)
result6 <- read.csv("11 independent pan-cancer nucleotides1.csv")
result6 <- result6[,-1]
dim(result6)
head(result6)
n11 <- result6
n11
result7 <- data.frame(Feature=c(result6$name),
                      HR=c(result6[,1]),
                      HR.lower=c(result6[,3]),
                      HR.upper=c(result6[,4]),
                      pvalue=round(c(result6[,2]),3))
dim(result7)
head(result7)

result7$HR1 <- paste(round(result7$HR,2),"(",round(result7$HR.lower,2),"-",round(result7$HR.upper,2),")",sep = '')
result7$pvalue1 <- ifelse(result7$pvalue<0.001,"***",ifelse(result7$pvalue<0.01,"**",ifelse(result7$pvalue<0.05,"*"," ")))
result7$pvalue2 <- paste(result7$pvalue,result7$pvalue1,sep = '')
dim(result7)
head(result7)

labeltext <- as.matrix(rbind(c("Feature","HR","pvalue"),result7[,c(1,6,8)]))
labeltext

forestplot(labeltext, 
           mean = c(NA,result7$HR),
           lower = c(NA,result7$HR.lower),
           upper = c(NA,result7$HR.upper),
           zero = 1,lwd.zero = 2,
           clip = c(0,5),
           graph.pos = 2,
           boxsize = 0.4,
           align = "l",
           xlog = F,
           col=fpColors(line = "black", zero = "black",
                        box="black"),
           lty.ci = 7,
           lwd.ci = 1,
           ci.vertices.height = 0.15,
           lwd.xaxis=2,
           xlab="Hazard Ratio",
           txt_gp = fpTxtGp(ticks = gpar(cex = 0.9), xlab = gpar(cex = 0.8), cex = 0.8),
           lineheight = unit(8,"mm")) 


######################### KEGG map distribution of 8 nucleotides
length(nucleotide13)
nucleotide13

head(nucleotide)
Nucleotide111 <- nucleotide[match(nucleotide13,nucleotide$compound_id),1]
Nucleotide111
Nucleotide112 <- nucleotide[match(nucleotide13,nucleotide$compound_id),5]
Nucleotide112
nucleotide$KEGG.ID


sum.144 <- setdiff(nucleotide$KEGG.ID,Nucleotide111)
length(sum.144)

sum.159 <- c(rep(0,144),rep(1,15))
sum.159
names(sum.159) <- c(sum.144,Nucleotide111)
sum.159

library(pathview)

purine.gene <- read.csv("purine gene list.csv")
purine.gene$Name

pyrimidine.gene <- read.csv("pyrimidine gene list.csv")
pyrimidine.gene

pathview(gene.data=as.character(purine.gene$Name),cpd.data = sum.159,
         pathway.id = '00230',species = "hsa",out.suffix="purine.gene.cpd.8.1",
         keys.align = "y", kegg.native = T,key.pos = "topright",na.col = "#E5F5E0",
         limit=list(gene = 1,cpd=1),mid = list(cpd="white"),
         high = list(gene="#E5F5E0",cpd="green"),
         discrete = list(gene=T,cpd=T),
         bins = list(gene=1,cpd=2))


pathview(gene.data=as.character(pyrimidine.gene$Name),cpd.data = sum.159,
         pathway.id = '00240',species = "hsa",out.suffix="pyrimidine.gene.cpd.8",
         keys.align = "y", kegg.native = T,key.pos = "topright",na.col = "#E5F5E0",
         limit=list(gene = 1,cpd=1),mid = list(cpd="white"),
         high = list(gene="#E5F5E0",cpd="green"),
         discrete = list(gene=T,cpd=T),
         bins = list(gene=1,cpd=2))







### predict survival outcome by ROC curve
library(pROC)
gfit <- roc(STATOS ~predict(res.cox2), data = res.cat1)

plot(gfit,
     print.auc=TRUE, #输出AUC值
     print.thres=TRUE, #输出cut-off值
     #main = "ROC Curve in Lung NAC tumor (n=77)", #设置图形的标题
     col= "red", #曲线颜色
     print.thres.col="black", #cut-off值字体的颜色
     identity.col="blue", #对角线颜色
     identity.lty=1,identity.lwd=1)



### Survival curve on risk score of cox model
head(res.cat1)
dim(res.cat1)
res.cat3 <- res.cat1[,c(1,2,4,6,9,10,12,15,17,18,20,24,27)]
dim(res.cat3)
res.cat3[1:5,]

# use Cox model to predict good/poor prognosis or low/high risk groups by ROC curve
res.cat3$STATOS1 <- ifelse(res.cat3$OS < median(res.cat3$OS),0,1) #0=short/poor survival or high risk,1=long/good survival or low risk  
res.cat3$STATOS1

gfit1 <- roc(STATOS1 ~predict(res.cox2), data = res.cat3)
cutoff <- gfit1$thresholds[which.max(gfit1$sensitivities+gfit1$specificities)]
cutoff 
plot(gfit1,
     print.auc=TRUE, #输出AUC值
     print.thres=TRUE, #输出cut-off值
     #main = "ROC Curve in Lung NAC tumor (n=77)", #设置图形的标题
     col= "red", #曲线颜色
     print.thres.col="black", #cut-off值字体的颜色
     identity.col="blue", #对角线颜色
     identity.lty=1,identity.lwd=1)




### cutoff-based cox model to stratify patients into different prognostic risk groups by survival curve
res.cat3[1:5,1:ncol(res.cat3)]
res.cat3$risk <- ifelse(predict(res.cox2) < cutoff,'low','high')
# fit survival curve
fit <- survfit(Surv(OS,STATOS)~risk,data =res.cat3)
fit
##绘制生存曲线##
ggsurvplot(fit, data = res.cat3, conf.int = F, # 增加置信区间
           risk.table = TRUE, pval=T) #不同时间点风险人数表，增加P值



### ROC curve on the prediction of 1/3/5-year survival rate
library(car)
library(rms)
library(pROC)
library(timeROC)
library(ggDCA)

head(res.cat1)
res.cat4 <- res.cat1[,c(1,2,4,6,9,10,12,15,17,18,20,24,27)]
res.cat4[1:5,1:ncol(res.cat4)]

res.cox4 <- cph(Surv(OS, STATOS) ~ `Methylmalonic acid` + `2',3'-CUMP` + AMP + UMP + `2',3'-Cyclic GMP` + 
                  `8-Hydroxyadenine`+Guanosine+`Cyclic AMP`+dAMP+dGDP+GDP,
                x=T,y=T,surv=T,data =res.cat4)
print(res.cox4)

pred_f_training <- predict(res.cox4,res.cat4,type="lp")#!!!type="lp",是他没错
data_table <- data.frame(time=res.cat4[,"OS"],status=res.cat4[,"STATOS"],score=pred_f_training)


time_roc_res <- timeROC(
  T = data_table$time,
  delta = data_table$status,
  marker = data_table$score,
  cause = 1,
  weighting="marginal",
  times = c(12, 36, 60),
  ROC = TRUE,
  iid = TRUE
)


time_ROC_df <- data.frame(
  TP_1year = time_roc_res$TP[, 1],
  FP_1year = time_roc_res$FP[, 1],
  TP_3year = time_roc_res$TP[, 2],
  FP_3year = time_roc_res$FP[, 2],
  TP_5year = time_roc_res$TP[, 3],
  FP_5year = time_roc_res$FP[, 3]
)


ggplot(data = time_ROC_df) +
  geom_line(aes(x = FP_1year, y = TP_1year), linewidth = 1, color = "#BC3C29FF") +
  geom_line(aes(x = FP_3year, y = TP_3year), linewidth = 1, color = "#0072B5FF") +
  geom_line(aes(x = FP_5year, y = TP_5year), linewidth = 1, color = "#E18727FF") +
  geom_abline(slope = 1, intercept = 0, color = "grey", linewidth = 1, linetype = 2) +
  theme_bw() +
  annotate("text",
           x = 0.75, y = 0.25, size = 4.5,
           label = paste0("AUC at 1 year = ", sprintf("%.3f", time_roc_res$AUC[[1]])), color = "#BC3C29FF"
  ) +
  annotate("text",
           x = 0.75, y = 0.15, size = 4.5,
           label = paste0("AUC at 3 years = ", sprintf("%.3f", time_roc_res$AUC[[2]])), color = "#0072B5FF"
  ) +
  annotate("text",
           x = 0.75, y = 0.05, size = 4.5,
           label = paste0("AUC at 5 years = ", sprintf("%.3f", time_roc_res$AUC[[3]])), color = "#E18727FF"
  ) +
  labs(x = "False positive rate", y = "True positive rate") +
  theme(
    axis.text = element_text(face = "bold", size = 11, color = "black"),
    axis.title.x = element_text(face = "bold", size = 14, color = "black"),
    axis.title.y = element_text(face = "bold", size = 14, color = "black"))


####################################### clustering analysis on independent nucleotides
head(res.cat)
dim(res.cat)

res.cat1 <- res.cat[,rownames(summ1)]
dim(res.cat1)

res.cat2 <- cbind(clin,res.cat1)
dim(res.cat2)
head(res.cat2)

res.cat2$STATOS1 <- ifelse(res.cat2$OS < median(res.cat2$OS),1,0) #1=short/poor survival or high risk,0=long/good survival or low risk  
res.cat2$STATOS1
dim(res.cat2)
head(res.cat2)

annotation_col <- data.frame(Type = res.cat2$Cancer,STATOS=res.cat2$STATOS,Risk=res.cat2$STATOS1,
                             Cluster=sum55$Cluster)
rownames(annotation_col) <- rownames(res.cat2)
head(annotation_col)
dim(annotation_col)
annotation_col$STATOS <- as.factor(annotation_col$STATOS)
annotation_col$Risk <- as.factor(annotation_col$Risk)
annotation_col$Cluster <- as.factor(annotation_col$Cluster)
str(annotation_col)

anncol = list(Risk=c("0"="gray","1"="purple"),
              STATOS=c("0"="gray","1"="red"),
              Cluster=c("1"="red","2"="green","3"="blue"))
anncol


res.cat3 <- res.cat2[,-c(1:3,12)]
res.cat4 <- distinct(as.data.frame(t(res.cat3)))
head(nucleotide)
rownames(res.cat4) <- nucleotide[match(rownames(res.cat4),nucleotide$compound_id),4]
dim(res.cat4)

library(pheatmap)
pheatmap1 <- pheatmap(res.cat4,cluster_cols =T,cluster_rows =T,
                      annotation_col = annotation_col,annotation_colors=anncol,
                      clustering_method = "ward.D2",
                      clustering_distance_rows = "euclidean",
                      clustering_distance_cols = "euclidean",
                      cutree_cols=3,
                      show_colnames = F,
                      show_rownames = T)

cluster.col <- pheatmap1$tree_col

#对聚类树进行分簇
cut <- cutree(cluster.col,3)

cut1 <- names(cut)[cut==1]
cut2 <- names(cut)[cut==2]
cut3 <- names(cut)[cut==3]
length(cut1)
length(cut2)
length(cut3)

############# survival analysis on patients group from different clusters
dim(clin)
head(clin)
sum55 <- clin
sum55$Cluster <- rep(1,1830)
sum55$Cluster[match(cut3,rownames(sum55))] <- rep(1,length(match(cut3,rownames(sum55))))
sum55$Cluster[match(cut1,rownames(sum55))] <- rep(2,length(match(cut1,rownames(sum55))))
sum55$Cluster[match(cut2,rownames(sum55))] <- rep(3,length(match(cut2,rownames(sum55))))
head(sum55)


# fit survival curve
fit <- survfit(Surv(OS,STATOS)~Cluster,data =sum55)
fit
survdiff(Surv(OS,STATOS)~Cluster,data =sum55) #采用log-rank 检验分析生存率差异

##绘制生存曲线##
ggsurvplot(fit, data = sum55, palette = c("red","green","blue"),conf.int = F, # 增加置信区间
           risk.table = TRUE, pval=T) #不同时间点风险人数表，增加P值

res.cox <- coxph(Surv(OS,STATOS)~Cluster,data =sum55)
summary(res.cox)
cox.zph(res.cox)


### survival analysis on any two subgroups
sum66 <- sum55[c(which(sum55$Cluster==1),which(sum55$Cluster==2)),] #***
sum66 <- sum55[c(which(sum55$Cluster==1),which(sum55$Cluster==3)),] #***
sum66 <- sum55[c(which(sum55$Cluster==2),which(sum55$Cluster==3)),] #***

# fit survival curve
fit <- survfit(Surv(OS,STATOS)~Cluster,data =sum66)
fit
survdiff(Surv(OS,STATOS)~Cluster,data =sum66) #采用log-rank 检验分析生存率差异

##绘制生存曲线##
ggsurvplot(fit, data = sum66, palette = c("green","blue"),conf.int = F, # 增加置信区间
           risk.table = TRUE, pval=T) #不同时间点风险人数表，增加P值

res.cox <- coxph(Surv(OS,STATOS)~Cluster,data =sum66)
summary(res.cox)
cox.zph(res.cox)


########## pan-cancer subgroup constitution
sub1 <- clin[cut1,]
sub2 <- clin[cut2,]
sub3 <- clin[cut3,]
dim(sub1)
dim(sub2)
dim(sub3)
head(sub3)

sub11 <- by(sub1,sub1$Cancer,count)
sub11 <- data.frame(name=names(sub11),count=as.numeric(sub11))
dim(sub11)
head(sub11)
sub22 <- by(sub2,sub2$Cancer,count)
sub22 <- data.frame(name=names(sub22),count=as.numeric(sub22))
dim(sub22)
head(sub22)
sub33 <- by(sub3,sub3$Cancer,count)
sub33 <- data.frame(name=names(sub33),count=as.numeric(sub33))
dim(sub33)
head(sub33)


sub <- by(clin,clin$Cancer,count)
sub <- data.frame(name=names(sub),count=as.numeric(sub))
colnames(sub)[2] <- "count1"
sub

sub111 <- merge(sub11,sub)
sub111$count2 <- sub111$count/sub111$count1
sub111$count3 <- sub111$count2/sum(sub111$count2)
sub111

sub222 <- merge(sub22,sub)
sub222$count2 <- sub222$count/sub222$count1
sub222$count3 <- sub222$count2/sum(sub222$count2)
sub222

sub333 <- merge(sub33,sub)
sub333$count2 <- sub333$count/sub333$count1
sub333$count3 <- sub333$count2/sum(sub333$count2)
sub333


sub.2 <- sub111
sub.2 <- sub222
sub.2 <- sub333
label_value <- paste(round(sub.2$count3/sum(sub.2$count3) * 100, 1), '%', sep = '')
label_value
ggplot(data=sub.2,aes(x='Content',y=count3,fill=name))+geom_bar(stat = 'identity',position = 'stack')+
  coord_polar(theta = 'y')+
  labs(x='',y='',title = '')+
  theme(axis.text = element_blank(),axis.ticks = element_blank(),panel.grid = element_blank(),panel.border = element_blank())+
  geom_text(aes(y = count3/2 + c(0, cumsum(count3)[-length(count3)]), x =1.2, label = label_value)) 


############ boxplot show differential analysis on significant nucleotide
length(cut1)
length(cut2)
length(cut3)
head(nucleotide)
dim(profile.29.com1)
head(profile.29.com1)

profile.29.com2 <- profile.29.com1[,summ1$HMDB]
dim(profile.29.com2)
colnames(profile.29.com2) <- nucleotide[match(colnames(profile.29.com2),nucleotide$compound_id),4]
head(profile.29.com2)
sum55 <- profile.29.com2
sum55$Cluster <- rep(1,1830)
sum55$Cluster[match(cut3,rownames(sum55))] <- rep("Cluster 1",length(match(cut3,rownames(sum55))))
sum55$Cluster[match(cut1,rownames(sum55))] <- rep("Cluster 2",length(match(cut1,rownames(sum55))))
sum55$Cluster[match(cut2,rownames(sum55))] <- rep("Cluster 3",length(match(cut2,rownames(sum55))))
head(sum55)

library(reshape2)
sum66 <- melt(sum55,id.vars = c("Cluster"))
colnames(sum66) <- c("Cluster","Nucleotide","Abundance")
dim(sum66)
head(sum66)

sum66 <- sum66[sum66$Abundance<4,]
sum66 <- sum66[sum66$Abundance> (-4),]
dim(sum66)

sum66$Cluster <- factor(sum66$Cluster,levels = c("Cluster 1","Cluster 2","Cluster 3"))

ggplot(sum66,aes(x=Nucleotide,y=Abundance,fill=Cluster))+stat_boxplot(geom="errorbar")+
  geom_boxplot(outlier.shape = NA)+
  ylab('Relative abundance')


# multiple comparison
sum77 <- sum66[sum66$Nucleotide==unique(sum66$Nucleotide)[8],]
dim(sum77)
head(sum77)

aov1 <- aov(Abundance~Cluster, sum77)

tuk <- TukeyHSD(aov1)
tuk

sum <- summary(aov1)
p <- sum[[1]][["Pr(>F)"]][1]
p




################# RF algorithm on independent nucleotides to predict survival outcome/rate and risk stratification:No run
dim(res.cat)
head(res.cat)

res.cat1 <- res.cat[,-1]
res.cat1$STATOS <- as.factor(res.cat1$STATOS)
dim(res.cat1)
head(res.cat1)

library(randomForest)

otu_group.forest <- randomForest(STATOS ~ ., data = res.cat1, importance = TRUE)
otu_group.forest


####交叉验证帮助选择特定数量
#5 次重复十折交叉验证
set.seed(12345)
otu_group.cv <- replicate(5, rfcv(res.cat1[,-1], res.cat1$STATOS, cv.fold = 10,step = 1.5), simplify = FALSE)
otu_group.cv

#提取验证结果绘图
otu_group.cv <- data.frame(sapply(otu_group.cv, '[[', 'error.cv'))
otu_group.cv$otus <- rownames(otu_group.cv)
otu_group.cv <- reshape2::melt(otu_group.cv, id = 'otus')
otu_group.cv$otus <- as.numeric(as.character(otu_group.cv$otus))

#拟合线图
library(ggplot2)
library(splines)  #用于在 geom_smooth() 中添加拟合线，或者使用 geom_line() 替代 geom_smooth() 绘制普通折线

p <- ggplot(otu_group.cv, aes(otus, value)) +
  geom_smooth(se = FALSE,  method = 'glm', formula = y~ns(x, 6)) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  labs(title = '',x = 'Number of OTUs', y = 'Cross-validation error')

p

#大约提取前 30 个重要的 OTUs
p + geom_vline(xintercept = 5)


#根据 OTUs 重要性排序后选择，例如根据“Mean Decrease Accuracy”指标
importance_otu <- data.frame(importance(otu_group.forest))
importance_otu <- importance_otu[order(importance_otu$MeanDecreaseAccuracy, decreasing = TRUE), ]
head(importance_otu)


res.cat1[1:5,1:10]
otu_select <- rownames(importance_otu)[1:5]
otu_group2 <- res.cat1[,c("STATOS",otu_select)]
dim(otu_group2)
head(otu_group2)

otu_group2.forest <- randomForest(STATOS ~ ., data = otu_group2, importance = TRUE)
otu_group2.forest
otu_group2.forest$confusion
otu_group2.forest$predicted



####################################### method2:combine all data after surv_cutpoint on each tumor:No run
dim(profile.19.com)
head(profile.19.com)

dim(clin)
head(clin)

profile.19.com1 <- cbind(clin[,1:2],profile.19.com)
dim(profile.19.com1)
head(profile.19.com1)
res.cat <- profile.19.com1

### univariate cox analysis to calculate p value, HR
result <- data.frame()
for (i in 3:ncol(res.cat)) {
  res.cox <- coxph(Surv(OS,STATOS)~res.cat[,i],data =res.cat)
  summ <- summary(res.cox)
  hypo <- cox.zph(res.cox)
  hypo1 <- hypo$table[-nrow(hypo$table),]
  result1 <- data.frame(HMDB=colnames(res.cat)[i],
                        HR = signif(summ$coef[2], digits=3),
                        HR.confint.lower = signif(summ$conf.int[,"lower .95"], 3),
                        HR.confint.upper = signif(summ$conf.int[,"upper .95"],3),
                        logtest.p=summ$logtest[3],
                        waldtest.p=summ$waldtest[3],
                        sctest.p=summ$sctest[3],
                        hypo.p=hypo1[3])
  result <- rbind(result,result1)
}
head(result)
dim(result)

# screen by sctest p value
result2 <- result[result$sctest.p<0.05,]
result2 <- result2[result2$hypo.p>0.05,]
dim(result2)
head(result2)

dim(nucleotide)
head(nucleotide)

result2$Nucleotide <- nucleotide[match(result2$HMDB,nucleotide$compound_id),4]
dim(result2)
head(result2)


# forest plot
library(forestplot)
result6 <- result2
dim(result6)
head(result6)
which(duplicated(result6$HMDB)==T)

result6$HR1 <- paste(result6$HR,"(",result6$HR.confint.lower,"-",result6$HR.confint.upper,")",sep = '')
result6$pvalue <- ifelse(result6$sctest.p<0.001,"<0.001",ifelse(result6$sctest.p<0.01,"<0.01","<0.05"))
dim(result6)
head(result6)

labeltext <- as.matrix(rbind(c("HMDB","Nucleotide","HR","pvalue"),result6[,c(1,9,10,11)]))
labeltext

forestplot(labeltext, 
           mean = c(NA,result6$HR),
           lower = c(NA,result6$HR.confint.lower),
           upper = c(NA,result6$HR.confint.upper),
           zero = 1,lwd.zero = 2,
           clip = c(0,10),
           graph.pos = 3,
           boxsize = 0.5,
           align = "l",
           xlog = FALSE,
           col=fpColors(line = "black", zero = "black",
                        box="black"),
           lty.ci = 7,
           lwd.ci = 1,
           ci.vertices.height = 0.15,
           lwd.xaxis=3,
           xlab="Hazard Ratio",
           txt_gp = fpTxtGp(ticks = gpar(cex = 0.9), xlab = gpar(cex = 0.8), cex = 0.8),
           lineheight = unit(8,"mm")) 


### multivariate analysis
head(res.cat)
dim(res.cat)
head(nucleotide)

res.cat1 <- res.cat[,c("OS","STATOS",result2$HMDB)]
dim(res.cat1)
colnames(res.cat1)[3:ncol(res.cat1)] <- nucleotide[match(colnames(res.cat1)[3:ncol(res.cat1)],nucleotide$compound_id),4]
head(res.cat1)
dim(res.cat1)

res.cox1 <- coxph(Surv(OS, STATOS) ~ .,data =res.cat1)
res.cox2 <- step(res.cox1)
ggforest(model=res.cox2,data=res.cat1)

#等比例风险假定
hypo <- cox.zph(res.cox2) 
hypo1 <- hypo$table[-nrow(hypo$table),]
dim(hypo1)
hypo2 <- hypo1[hypo1[,3]>0.05,]
dim(hypo2)
hypo2


### forest plot
summ <- as.data.frame(cbind(summary(res.cox2)$coefficients[,c(2,5)],summary(res.cox2)$conf.int[,c(3,4)]))
summ$mz.value <- rownames(summ)
head(summ)
dim(summ)


library(forestplot)

result6 <- summ
dim(result6)
head(result6)

result7 <- data.frame(Feature=c(result6$mz.value),
                      HR=c(result6[,1]),
                      HR.lower=c(result6[,3]),
                      HR.upper=c(result6[,4]),
                      pvalue=round(c(result6[,2]),3))
dim(result7)
head(result7)

result7$HR1 <- paste(round(result7$HR,2),"(",round(result7$HR.lower,2),"-",round(result7$HR.upper,2),")",sep = '')
result7$pvalue1 <- ifelse(result7$pvalue<0.001,"***",ifelse(result7$pvalue<0.01,"**",ifelse(result7$pvalue<0.05,"*"," ")))
result7$pvalue2 <- paste(result7$pvalue,result7$pvalue1,sep = '')
dim(result7)
head(result7)

labeltext <- as.matrix(rbind(c("Feature","HR","pvalue"),result7[,c(1,6,8)]))
labeltext

forestplot(labeltext, 
           mean = c(NA,result7$HR),
           lower = c(NA,result7$HR.lower),
           upper = c(NA,result7$HR.upper),
           zero = 1,lwd.zero = 2,
           clip = c(0,5),
           graph.pos = 2,
           boxsize = 0.4,
           align = "l",
           xlog = F,
           col=fpColors(line = "black", zero = "black",
                        box="black"),
           lty.ci = 7,
           lwd.ci = 1,
           ci.vertices.height = 0.15,
           lwd.xaxis=2,
           xlab="Hazard Ratio",
           txt_gp = fpTxtGp(ticks = gpar(cex = 0.9), xlab = gpar(cex = 0.8), cex = 0.8),
           lineheight = unit(8,"mm")) 






############################ Pan-cancer analysis on chemotherapy #############################

setwd("C:/MyRdata8/Summary")
#write.csv(profile.29.com1,"profile.29.com1.csv")
profile.29.com1 <- read.csv("profile.29.com1.csv")
rownames(profile.29.com1) <- profile.29.com1[,1]
profile.29.com1 <- profile.29.com1[,-1]
dim(profile.29.com1)
head(profile.29.com1)

profile.29.com1 <- profile.29.com1[-which(profile.29.com1$Cancer=="HER2-Gastric cancer"),]
profile.29.com1$Cancer <- ifelse(profile.29.com1$Cancer=="Lung NAC Tumor"|profile.29.com1$Cancer=="Neoadjuvant treated EAC","Chemotherapy","Nontreatment")
head(profile.29.com1)
dim(profile.29.com1)

################################## survival curve on Chemotherapy vs Nontreatment
# fit survival curve
fit <- survfit(Surv(OS,STATOS)~Cancer,data =profile.29.com1)
fit
survdiff(Surv(OS,STATOS)~Cancer,data =profile.29.com1) #采用log-rank 检验分析生存率差异

##绘制生存曲线##
ggsurvplot(fit, data = profile.29.com1, conf.int = F, # 增加置信区间
           risk.table = TRUE, pval=T) #不同时间点风险人数表，增加P值

res.cox <- coxph(Surv(OS,STATOS)~Cancer,data =profile.29.com1)
summary(res.cox)
cox.zph(res.cox)



################################### prognosis analysis on 29 nucleotides between chemotherapy and non-chemotherapy
dim(profile.29.com1)
head(profile.29.com1)


profile.29.com1.no <- profile.29.com1[profile.29.com1$Cancer=="Nontreatment",]
profile.29.com1.no <- profile.29.com1.no[,-3]
dim(profile.29.com1.no)
head(profile.29.com1.no)
res.cut <- surv_cutpoint(profile.29.com1.no, time = "OS", event = "STATOS",
                         minprop = 0.1,
                         variables = colnames(profile.29.com1.no)[3:ncol(profile.29.com1.no)])
res.cat.no <- surv_categorize(res.cut,labels = c(0, 1))
head(res.cat.no)

# univariate cox analysis to calculate p value, HR
result <- data.frame()
for (i in 3:ncol(res.cat.no)) {
  res.cox <- coxph(Surv(OS,STATOS)~res.cat.no[,i],data =res.cat.no)
  summ <- summary(res.cox)
  hypo <- cox.zph(res.cox)
  hypo1 <- hypo$table[-nrow(hypo$table),]
  result1 <- data.frame(HMDB=colnames(res.cat.no)[i],
                        HR = signif(summ$coef[2], digits=3),
                        HR.confint.lower = signif(summ$conf.int[,"lower .95"], 3),
                        HR.confint.upper = signif(summ$conf.int[,"upper .95"],3),
                        logtest.p=summ$logtest[3],
                        waldtest.p=summ$waldtest[3],
                        sctest.p=summ$sctest[3],
                        hypo.p=hypo1[3])
  result <- rbind(result,result1)
}
head(result)
dim(result)
result.no <- result

# screen by sctest p value
result2 <- result.no[result.no$sctest.p<0.05,]
result2.no <- result2[result2$hypo.p>0.05,]
dim(result2.no)
head(result2.no)




profile.29.com1.yes <- profile.29.com1[profile.29.com1$Cancer=="Chemotherapy",]
profile.29.com1.yes <- profile.29.com1.yes[,-3]
dim(profile.29.com1.yes)
head(profile.29.com1.yes)
res.cut <- surv_cutpoint(profile.29.com1.yes, time = "OS", event = "STATOS",
                         minprop = 0.1,
                         variables = colnames(profile.29.com1.yes)[3:ncol(profile.29.com1.yes)])
res.cat.yes <- surv_categorize(res.cut,labels = c(0, 1))
head(res.cat.yes)
# univariate cox analysis to calculate p value, HR
result <- data.frame()
for (i in 3:ncol(res.cat.yes)) {
  res.cox <- coxph(Surv(OS,STATOS)~res.cat.yes[,i],data =res.cat.yes)
  summ <- summary(res.cox)
  hypo <- cox.zph(res.cox)
  hypo1 <- hypo$table[-nrow(hypo$table),]
  result1 <- data.frame(HMDB=colnames(res.cat.yes)[i],
                        HR = signif(summ$coef[2], digits=3),
                        HR.confint.lower = signif(summ$conf.int[,"lower .95"], 3),
                        HR.confint.upper = signif(summ$conf.int[,"upper .95"],3),
                        logtest.p=summ$logtest[3],
                        waldtest.p=summ$waldtest[3],
                        sctest.p=summ$sctest[3],
                        hypo.p=hypo1[3])
  result <- rbind(result,result1)
}
head(result)
dim(result)
result.yes <- result

# screen by sctest p value
result2 <- result.yes[result.yes$sctest.p<0.05,]
result2.yes <- result2[result2$hypo.p>0.05,]
dim(result2.yes)
head(result2.yes)


length(unique(c(result2.yes$HMDB,result2.no$HMDB)))
length(intersect(result2.yes$HMDB,result2.no$HMDB))
com <- unique(c(result2.yes$HMDB,result2.no$HMDB))
com

# forest plot only with mass value
library(forestplot)

result6 <- result.no[match(unique(c(result18$name,com)),result.no$HMDB),]
result6 <- result.yes[match(unique(c(result18$name,com)),result.yes$HMDB),]
dim(result6)
head(result6)
head(nucleotide)
result6$nucleotide <- nucleotide[match(result6$HMDB,nucleotide$compound_id),5]
head(result6)
Nucleotide11 <- result6$HMDB
length(Nucleotide11)

result7 <- data.frame(Feature=c(result6$nucleotide),
                      HR=c(result6$HR),
                      HR.lower=c(result6$HR.confint.lower),
                      HR.upper=c(result6$HR.confint.upper),
                      pvalue=round(c(result6$sctest.p),3))
dim(result7)
head(result7)
result7 <- distinct(result7)
result7 <- result7[-15,]
n20 <- result7

result7$HR1 <- paste(result7$HR,"(",result7$HR.lower,"-",result7$HR.upper,")",sep = '')
dim(result7)
head(result7)

labeltext <- as.matrix(rbind(c("Feature","HR","pvalue"),result7[,c(1,2,5)]))
labeltext

forestplot(labeltext[,1], 
           mean = c(NA,result7$HR),
           lower = c(NA,result7$HR.lower),
           upper = c(NA,result7$HR.upper),
           zero = 1,lwd.zero = 2,
           clip = c(0,5),
           graph.pos = 2,
           boxsize = 0.5,
           align = "l",
           xlog = FALSE,
           col=fpColors(line = "black", zero = "black",
                        box="black"),
           lty.ci = 7,
           lwd.ci = 1,
           ci.vertices.height = 0.15,
           lwd.xaxis=2,
           xlab="Hazard Ratio",
           txt_gp = fpTxtGp(ticks = gpar(cex = 0.9), xlab = gpar(cex = 0.8), cex = 0.8),
           lineheight = unit(8,"mm")) 


######################### KEGG map distribution of 25 nucleotides
length(Nucleotide11)
Nucleotide11

head(nucleotide)
Nucleotide1111 <- nucleotide[match(Nucleotide11,nucleotide$compound_id),1]
Nucleotide1111
nucleotide$KEGG.ID


sum.136 <- setdiff(nucleotide$KEGG.ID,Nucleotide1111)
sum.136

sum.159 <- c(rep(0,134),rep(1,25))
sum.159
names(sum.159) <- c(sum.136,Nucleotide1111)
sum.159

library(pathview)

pathview(cpd.data = sum.159,pathway.id = '00230',species = "hsa",
         limit=list(cpd=1),low=list(cpd="red"),mid = list(cpd="red"),
         high = list(cpd="green"),discrete = list(cpd=T),bins = list(cpd=2))
pathview(cpd.data = sum.159,pathway.id = '00240',species = "hsa",
         limit=list(cpd=1),low=list(cpd="red"),mid = list(cpd="red"),
         high = list(cpd="green"),discrete = list(cpd=T),bins = list(cpd=2))




######################### abundance difference of significant nucleotides in two groups 
dim(profile.29.com1)
head(profile.29.com1)

profile.29.com2 <- profile.29.com1
dim(profile.29.com2)
head(profile.29.com2)

res.cut <- surv_cutpoint(profile.29.com2, time = "OS", event = "STATOS",
                         minprop = 0.1,
                         variables = colnames(profile.29.com2)[4:ncol(profile.29.com2)])
res.cat <- surv_categorize(res.cut,labels = c(0, 1))
head(res.cat)
dim(res.cat)

res.cat1 <- cbind(profile.29.com2[,c(1:3)],res.cat[,3:ncol(res.cat)])
dim(res.cat1)
head(res.cat1)

############ chisq test
dim(res.cat1)
head(res.cat1)
res.cat11 <- res.cat1
head(res.cat11)

result17 <- data.frame()
for (i in 4:ncol(res.cat11)) {
  chis <- chisq.test(table(res.cat11$Cancer,res.cat11[,i]))
  table1 <- table(res.cat11$Cancer,res.cat11[,i])
  result18 <- data.frame(name=colnames(res.cat11)[i],
                         pvalue=chis$p.value,
                         Nontreatment=table1[2,2]/table1[2,1],
                         Chemotherapy=table1[1,2]/table1[1,1],
                         diff=(table1[1,2]/table1[1,1]-table1[2,2]/table1[2,1]))
  result17 <- rbind(result17,result18)
}
dim(result17)
result17
result17$abundance <- ifelse(result17$Nontreatment>1&result17$diff>0,"high-abundant in Chemotherapy",
                             ifelse(result17$Nontreatment>1&result17$diff<0,"high-abundant in Nontreatment",
                                    ifelse(result17$Nontreatment<1&result17$diff<0,"high-abundant in Chemotherapy","high-abundant in Nontreatment")))
result17$significance <- ifelse(result17$pvalue>0.05,"NS",ifelse(result17$pvalue<0.001,"***",ifelse(result17$pvalue<0.01,"**","*")))
head(nucleotide)
result17$Name <- nucleotide[match(result17$name,nucleotide$compound_id),5]
head(result17)
#write.csv(result17,"abundance analysis of 29 nucleotides in 2 groups.csv")

result18 <- result17[result17$pvalue < 0.05,]
dim(result18)
head(result18)

unique(c(result18$name,com))
unique(result18$Name)





############################ Chemotherapy response related nucleotides in 2 cohorts by forest plot and RF ######################


setwd("C:/MyRdata8/Summary")
#write.csv(sum,"summary on single metabolite from univariate cox analysis5.csv")
sum <- read.csv("summary on single metabolite from univariate cox analysis5.csv")
sum <- sum[,-1]
univariate.sum <- sum
dim(sum)
head(sum)


### 67 common nucleotides
setwd("C:/MyRdata8/Summary")
sum2 <- read.csv("annotated and detectable 105 nucleotides in all cancers.csv")
sum3 <- sum2[,-1]
head(sum3)
dim(sum3)
length(unique(sum3$name)) #105 nucleotides
unique(sum3$tumor)
sum4.2 <- unique(sum3[which(sum3$tumor==unique(sum3$tumor)[7]),9])
sum4.4 <- unique(sum3[which(sum3$tumor==unique(sum3$tumor)[11]),9])

sum5 <- intersect(sum4.2, sum4.4)
#sum5 <- setdiff(sum5,c("FGAM"))
length(sum5)
sum5

sum.46.com1 <- data.frame(name=sum5)
sum.46.com1

dim(univariate.sum)
head(univariate.sum)


# annotation of significant 46 nucleotides in all cancers
univariate.sum.46 <- merge(sum.46.com1,univariate.sum)
dim(univariate.sum.46)
head(univariate.sum.46)


tumor <- unique(sum3$tumor)[c(7,11)]
tumor


# consider univariate significance and intensity
sum10 <- data.frame()
for (i in 1:2) {
  sum4 <- sum3[which(sum3$tumor==tumor[i]),]
  univariate.sum.46.1 <- univariate.sum.46[which(univariate.sum.46$Tumor==tumor[i]),c(1,3,24)]
  sum5 <- data.frame(name=setdiff(sum.46.com1$name,univariate.sum.46.1$name))
  sum6 <- merge(sum5,sum4)
  sum6$name.adduct <- paste(sum6$adduct,sum6$name,sep = '.')
  sum6 <- sum6[order(sum6$ppm,decreasing = F),]
  sum7 <- sum6[which(duplicated(sum6$name.adduct)==F),]
  sum7 <- sum7[order(sum7$intensity,decreasing = T),]
  sum8 <- sum7[which(duplicated(sum7$name)==F),]
  sum9 <- data.frame(name=c(sum8$name,univariate.sum.46.1$name),
                     mass=c(sum8$mass,univariate.sum.46.1$mz.value),
                     tumor=c(sum8$tumor,univariate.sum.46.1$Tumor))
  sum10 <- rbind(sum10,sum9)
}
dim(sum10)
head(sum10)


############################### Lung NAC tumor (n=77) 
sum11 <- sum10[sum10$tumor==tumor[1],]
head(sum11)
dim(sum11)

dim(table3.n77)
table3.n77[1:5,1:5]

table3.n77.1 <- table3.n77[,match(sum11$mass,colnames(table3.n77))]
colnames(table3.n77.1) <- sum11$name
dim(table3.n77.1)
table3.n77.1[1:5,1:5]
table3.n77.2 <- scale(table3.n77.1,center = T,scale = T)
head(table3.n77.2)
dim(table3.n77.2)
table3.n77.3 <- cbind(clin.n77[,c("OS","STATOS")],table3.n77.2)
table3.n77.3$response <- ifelse(clin.n77$Tumor.Regression=="Non-responder",1,0)
table3.n77.3$tumor <- rep("Lung NAC tumor",77)
dim(table3.n77.3)


############################### Neoadjuvant treated EACs (n=144)
sum11 <- sum10[sum10$tumor==tumor[2],]
head(sum11)
dim(sum11)

dim(table3.n144)
table3.n144[1:5,1:5]

table3.n144.1 <- table3.n144[,match(sum11$mass,colnames(table3.n144))]
colnames(table3.n144.1) <- sum11$name
dim(table3.n144.1)
table3.n144.1[1:5,1:5]
table3.n144.2 <- scale(table3.n144.1,center = T,scale = T)
head(table3.n144.2)
table3.n144.3 <- cbind(res.cat.Neoadjuvant.treated.EACs.n144[,c("OS","STATOS")],table3.n144.2)
table3.n144.3$response <- surv.n144$Becker.Resp.vs.Non.resp
table3.n144.3$tumor <- rep("Neoadjuvant treated EACs",144)
dim(table3.n144.3)
head(table3.n144.3)

########################## summary on 3 cancers define cutoff after scale
sum123 <- rbind(table3.n77.3,table3.n144.3)
dim(sum123)
head(sum123)

### survival curve on treatment response
# fit survival curve
fit <- survfit(Surv(OS,STATOS)~response,data =sum123)
fit # median os 67 vs 26.6 months
survdiff(Surv(OS,STATOS)~response,data =sum123) #采用log-rank 检验分析生存率差异

##绘制生存曲线##
ggsurvplot(fit, data = sum123, conf.int = F, # 增加置信区间
           risk.table = TRUE, pval=T) #不同时间点风险人数表，增加P值

res.cox <- coxph(Surv(OS,STATOS)~response,data =sum123)
summary(res.cox)
cox.zph(res.cox)

# prediction of survival by response
library(pROC)
gfit <- roc(STATOS~response, data = sum123)

plot(gfit,
     print.auc=TRUE, #输出AUC值
     print.thres=TRUE, #输出cut-off值
     #main = "ROC Curve in Lung NAC tumor (n=77)", #设置图形的标题
     col= "red", #曲线颜色
     print.thres.col="black", #cut-off值字体的颜色
     identity.col="blue", #对角线颜色
     identity.lty=1,identity.lwd=1)

table(sum123$STATOS,sum123$response)
accuracy <- sum(diag(table(sum123$STATOS,sum123$response)))/sum(table(sum123$STATOS,sum123$response))
accuracy


#################### treatment response related prognostic factor
### cutoff by surv_cutpoint
dim(sum123)
head(sum123)
sum123.good <- sum123[sum123$response==0,1:69]
dim(sum123.good)
head(sum123.good)
sum123.poor <- sum123[sum123$response==1,1:69]
dim(sum123.poor)
head(sum123.poor)

sum8 <- sum123.good
head(nucleotide)
colnames(sum8)[3:69] <- nucleotide[match(colnames(sum8)[3:69],nucleotide$name),6]
head(sum8)
dim(sum8)

res.cut <- surv_cutpoint(sum8, time = "OS", event = "STATOS",
                         minprop = 0.1,
                         variables = colnames(sum8)[3:ncol(sum8)])
res.cat <- surv_categorize(res.cut,labels = c(0, 1))
head(res.cat)

# univariate cox analysis to calculate p value, HR
result <- data.frame()
for (i in 3:ncol(res.cat)) {
  res.cox <- coxph(Surv(OS,STATOS)~res.cat[,i],data =res.cat)
  summ <- summary(res.cox)
  hypo <- cox.zph(res.cox)
  hypo1 <- hypo$table[-nrow(hypo$table),]
  result1 <- data.frame(HMDB=colnames(res.cat)[i],
                        HR = signif(summ$coef[2], digits=3),
                        HR.confint.lower = signif(summ$conf.int[,"lower .95"], 3),
                        HR.confint.upper = signif(summ$conf.int[,"upper .95"],3),
                        logtest.p=summ$logtest[3],
                        waldtest.p=summ$waldtest[3],
                        sctest.p=summ$sctest[3],
                        hypo.p=hypo1[3])
  result <- rbind(result,result1)
}
head(result)
dim(result)
result.good <- result
result.ann <- result
head(nucleotide)
result.ann$name <- nucleotide[match(result.ann$HMDB,nucleotide$compound_id),5]
head(result.ann)
dim(result.ann)

# screen by sctest p value
result2 <- result[result$sctest.p<0.05,]
result3 <- result2[result2$hypo.p>0.05,]
dim(result3)
head(result3)
result3.good <- result3


sum8 <- sum123.poor
head(nucleotide)
colnames(sum8)[3:69] <- nucleotide[match(colnames(sum8)[3:69],nucleotide$name),6]
head(sum8)
dim(sum8)

res.cut <- surv_cutpoint(sum8, time = "OS", event = "STATOS",
                         minprop = 0.1,
                         variables = colnames(sum8)[3:ncol(sum8)])
res.cat <- surv_categorize(res.cut,labels = c(0, 1))
head(res.cat)

# univariate cox analysis to calculate p value, HR
result <- data.frame()
for (i in 3:ncol(res.cat)) {
  res.cox <- coxph(Surv(OS,STATOS)~res.cat[,i],data =res.cat)
  summ <- summary(res.cox)
  hypo <- cox.zph(res.cox)
  hypo1 <- hypo$table[-nrow(hypo$table),]
  result1 <- data.frame(HMDB=colnames(res.cat)[i],
                        HR = signif(summ$coef[2], digits=3),
                        HR.confint.lower = signif(summ$conf.int[,"lower .95"], 3),
                        HR.confint.upper = signif(summ$conf.int[,"upper .95"],3),
                        logtest.p=summ$logtest[3],
                        waldtest.p=summ$waldtest[3],
                        sctest.p=summ$sctest[3],
                        hypo.p=hypo1[3])
  result <- rbind(result,result1)
}
head(result)
dim(result)
result.poor <- result

# screen by sctest p value
result2 <- result[result$sctest.p<0.05,]
result3 <- result2[result2$hypo.p>0.05,]
dim(result3)
head(result3)
result3.poor <- result3


com <- unique(c(result3.good$HMDB,result3.poor$HMDB))
com


# forest plot only with mass value
library(forestplot)

result6 <- result.good[match(unique(c(result18$name,com)),result.good$HMDB),]
result6 <- result.poor[match(unique(c(result18$name,com)),result.poor$HMDB),]
dim(result6)
head(result6)
head(nucleotide)
result6$nucleotide <- nucleotide[match(result6$HMDB,nucleotide$compound_id),5]
head(result6)
nucleotide12 <- result6$nucleotide

result7 <- data.frame(Feature=c(result6$nucleotide),
                      HR=c(result6$HR),
                      HR.lower=c(result6$HR.confint.lower),
                      HR.upper=c(result6$HR.confint.upper),
                      pvalue=round(c(result6$sctest.p),3))
dim(result7)
head(result7)
result7 <- distinct(result7)
n49 <- result7

result7$HR1 <- paste(result7$HR,"(",result7$HR.lower,"-",result7$HR.upper,")",sep = '')
dim(result7)
head(result7)

labeltext <- as.matrix(rbind(c("Feature","HR","pvalue"),result7[,c(1,2,5)]))
labeltext

forestplot(labeltext[,1], 
           mean = c(NA,result7$HR),
           lower = c(NA,result7$HR.lower),
           upper = c(NA,result7$HR.upper),
           zero = 1,lwd.zero = 2,
           clip = c(0,5),
           graph.pos = 2,
           boxsize = 0.5,
           align = "l",
           xlog = FALSE,
           col=fpColors(line = "black", zero = "black",
                        box="black"),
           lty.ci = 7,
           lwd.ci = 1,
           ci.vertices.height = 0.15,
           lwd.xaxis=2,
           xlab="Hazard Ratio",
           txt_gp = fpTxtGp(ticks = gpar(cex = 0.9), xlab = gpar(cex = 0.8), cex = 0.5),
           lineheight = unit(8,"mm")) 


######################### KEGG map distribution of 24 nucleotides
length(nucleotide12)
nucleotide12

head(nucleotide)
Nucleotide11111 <- nucleotide[match(nucleotide12,nucleotide$name),1]
Nucleotide11111
nucleotide$KEGG.ID


sum.135 <- setdiff(nucleotide$KEGG.ID,Nucleotide11111)
sum.135

sum.159 <- c(rep(0,135),rep(1,24))
sum.159
names(sum.159) <- c(sum.135,Nucleotide11111)
sum.159

library(pathview)

pathview(cpd.data = sum.159,pathway.id = '00230',species = "hsa",
         limit=list(cpd=1),low=list(cpd="red"),mid = list(cpd="red"),
         high = list(cpd="green"),discrete = list(cpd=T),bins = list(cpd=2))
pathview(cpd.data = sum.159,pathway.id = '00240',species = "hsa",
         limit=list(cpd=1),low=list(cpd="red"),mid = list(cpd="red"),
         high = list(cpd="green"),discrete = list(cpd=T),bins = list(cpd=2))



######################### abundance difference of significant nucleotides in two cohorts 
dim(sum123)
head(sum123)
head(nucleotide)
sum1234 <- sum123[,c(1,2,70,3:69)]
colnames(sum1234)[4:ncol(sum1234)] <- nucleotide[match(colnames(sum1234)[4:ncol(sum1234)],nucleotide$name),6]
dim(sum1234)
head(sum1234)

res.cut <- surv_cutpoint(sum1234, time = "OS", event = "STATOS",
                         minprop = 0.1,
                         variables = colnames(sum1234)[4:ncol(sum1234)])
res.cat <- surv_categorize(res.cut,labels = c(0, 1))
head(res.cat)
dim(res.cat)

res.cat1 <- cbind(sum1234[,c(1:3)],res.cat[,3:69])
dim(res.cat1)
head(res.cat1)



### chisq test
dim(res.cat1)
head(res.cat1)
i=4

result17 <- data.frame()
for (i in 4:70) {
  chis <- chisq.test(table(res.cat1$response,res.cat1[,i]))
  table1 <- table(res.cat1$response,res.cat1[,i])
  result18 <- data.frame(name=colnames(res.cat1)[i],
                         pvalue=chis$p.value,
                         Nonresponder=table1[2,2]/table1[2,1],
                         Responder=table1[1,2]/table1[1,1],
                         diff=(table1[2,2]/table1[2,1]-table1[1,2]/table1[1,1]))
  result17 <- rbind(result17,result18)
}
dim(result17)
result17
result17$abundance <- ifelse(result17$Nonresponder>1&result17$diff>0,"high-abundant in Nonresponder",
                             ifelse(result17$Nonresponder>1&result17$diff<0,"high-abundant in Responder",
                                    ifelse(result17$Nonresponder<1&result17$diff<0,"high-abundant in Nonresponder","high-abundant in Responder")))
result17$significance <- ifelse(result17$pvalue>0.05,"NS",ifelse(result17$pvalue<0.001,"***",ifelse(result17$pvalue<0.01,"**","*")))
result17


result18 <- result17[result17$pvalue < 0.05,]
dim(result18)
head(result18)
head(nucleotide)
result18$Name <- nucleotide[match(result18$name,nucleotide$compound_id),5]
#write.csv(result18,"abundance analysis of 67 nucleotides in 2 cohorts.csv")
unique(c(result18$name,com))



################################ RF to predict survival outcome comparing with treatment response
sum123 <- rbind(table3.n77.3,table3.n144.3)
dim(sum123)
head(sum123)

head(nucleotide)
colnames(sum123)[3:69] <- nucleotide[match(colnames(sum123)[3:69],nucleotide$name),6]


sum123.1 <- sum123[,1:69]
res.cut <- surv_cutpoint(sum123.1, time = "OS", event = "STATOS",
                         minprop = 0.1,
                         variables = colnames(sum123.1)[3:ncol(sum123.1)])
res.cat <- surv_categorize(res.cut,labels = c(0, 1))
head(res.cat)
dim(res.cat)
res.cat1 <- res.cat[,-c(1:2)]
res.cat1$STATOS <- as.factor(sum123$STATOS) # classifier predict survival outcome
#res.cat1$STATOS <- as.factor(sum123$Survival) # classifier predict survival time long or short
dim(res.cat1)
head(res.cat1)

res.cat2  <- as.data.frame(t(res.cat1))
dim(res.cat2)
res.cat3 <- distinct(res.cat2)
dim(res.cat3)
res.cat4 <- as.data.frame(t(res.cat3))
dim(res.cat4)
head(res.cat4)
res.cat4$STATOS <- as.factor(res.cat4$STATOS)

library(randomForest)
set.seed(1234567)
otu_group.forest <- randomForest(STATOS ~ ., data = res.cat4, importance = TRUE)
otu_group.forest

set.seed(1234567)
otu_group.cv <- replicate(5, rfcv(res.cat4[,-54], res.cat4$STATOS, cv.fold = 10,step = 1.5), simplify = FALSE)
otu_group.cv

otu_group.cv <- data.frame(sapply(otu_group.cv, '[[', 'error.cv'))
otu_group.cv$otus <- rownames(otu_group.cv)
otu_group.cv <- reshape2::melt(otu_group.cv, id = 'otus')
otu_group.cv$otus <- as.numeric(as.character(otu_group.cv$otus))

library(ggplot2)
library(splines)  #用于在 geom_smooth() 中添加拟合线，或者使用 geom_line() 替代 geom_smooth() 绘制普通折线

p <- ggplot(otu_group.cv, aes(otus, value)) +
  geom_smooth(se = FALSE,  method = 'glm', formula = y~ns(x, 6)) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  labs(title = '',x = 'Number of OTUs', y = 'Cross-validation error')

p
p + geom_vline(xintercept = 40)
importance_otu <- data.frame(importance(otu_group.forest))
importance_otu <- importance_otu[order(importance_otu$MeanDecreaseAccuracy, decreasing = TRUE), ]
head(importance_otu)
dim(importance_otu)

otu_select <- rownames(importance_otu)[1:53]
otu_group2 <- res.cat4[,c("STATOS",otu_select)]
dim(otu_group2)
head(otu_group2)

accuracy1 <- data.frame()
for (i in 1:1000) {
  set.seed(i)
  otu_group2.forest <- randomForest(STATOS ~ ., data = otu_group2, importance = TRUE)
  table(otu_group2.forest$predicted,otu_group2$STATOS)
  accuracy <- sum(diag(table(otu_group2.forest$predicted,otu_group2$STATOS)))/sum(table(otu_group2.forest$predicted,otu_group2$STATOS))
  accuracy
  accuracy2 <- data.frame(Accuracy=accuracy,I=i)
  accuracy1 <- rbind(accuracy1,accuracy2)
}
accuracy1
median(accuracy1$Accuracy)
table(accuracy1$Accuracy)
i=139

sum123.1$newclass <- otu_group2.forest$predicted
dim(sum123.1)
head(sum123.1)
# fit survival curve
fit <- survfit(Surv(OS,STATOS)~newclass,data =sum123.1)
fit # median os 67 vs 21 months
##绘制生存曲线##
ggsurvplot(fit, data = sum123.1, conf.int = F, # 增加置信区间
           risk.table = TRUE, pval=T) #不同时间点风险人数表，增加P值

res.cox <- coxph(Surv(OS,STATOS)~newclass,data =sum123.1)
summary(res.cox)
cox.zph(res.cox)


importance_otu <- data.frame(importance(otu_group2.forest))
importance_otu <- importance_otu[order(importance_otu$MeanDecreaseGini, decreasing = TRUE), ]
head(importance_otu)
dim(importance_otu)
importance_otu$MeanDecreaseGini

head(nucleotide)
rownames(importance_otu) <- nucleotide[match(rownames(importance_otu),nucleotide$compound_id),5]
importance_otu$Nucleotide <- rownames(importance_otu)
importance_otu$Nucleotide <- factor(importance_otu$Nucleotide,levels = importance_otu$Nucleotide)
head(importance_otu)
#write.csv(importance_otu,"67 nucleotides for classifier in cancer chemotherapy.csv")
ggplot(data=importance_otu,aes(x=Nucleotide,y=MeanDecreaseGini))+
  geom_bar(stat="identity")+
  theme(axis.text.x = element_text(size = 10,vjust=0.5,hjust=1,angle=90))+
  scale_y_continuous(breaks=seq(-4,12,2))

n53 <- importance_otu$Nucleotide


## correlation of 67 nucleotides abundance and chemotherapy responses
dim(otu_group2)
head(otu_group2)

res <- data.frame()
for (i in 2:68) {
  cor1 <- cor.test(otu_group2[,i],as.numeric(otu_group2$STATOS),method = "kendall")
  res1 <- data.frame(corr=cor1$estimate,pvalue=cor1$p.value,hmdb=colnames(otu_group2)[i])
  res <- rbind(res,res1)
}
dim(res)
head(res)
res$nucleotide <- nucleotide[match(res$hmdb,nucleotide$compound_id),4]
res2 <- res[res$pvalue<0.05,]
dim(res2)


library(reshape2)
library(plyr)
otu_group3 <- otu_group2[,c("STATOS",res2$hmdb)]
colnames(otu_group3)[2:24] <- res2$nucleotide
head(otu_group3)
dim(otu_group3)

i=i+1
otu_group4 <- table(otu_group3$STATOS,otu_group3[,i])
rownames(otu_group4) <- c("Sensitive","Resistant")
colnames(otu_group4) <- c("Low","High")
otu_group4

otu_group5 <- melt(otu_group4)
colnames(otu_group5) <- c("Response","Abundance","Count")
otu_group5

ggplot(data=otu_group5,aes(x=Abundance,y=Count,fill=Response))+
  geom_bar(stat="identity")+
  geom_text(aes(label=Count))+
  annotate("text",x=2,y=150,label=round(res2$pvalue[i-1],5))+
  annotate("text",x=2,y=170,label=res2$nucleotide[i-1])+
  annotate("text",x=2,y=160,label=round(res2$corr[i-1],3))

i


############################################ survival curve of UICC on survival
dim(clin.n77)
clin.n77[1:5,1:10]

dim(clin.n144)
clin.n144[1:5,1:10]

clin1 <- rbind(clin.n77[,c("UICC","OS","STATOS")],clin.n144[,c("UICC","OS","STATOS")])
dim(clin1)
head(clin1)

clin1$UICC <- ifelse(clin1$UICC<3,0,1)

# fit survival curve
fit <- survfit(Surv(OS,STATOS)~UICC,data =clin1)
fit #median os 67 vs 25 months
##绘制生存曲线##
ggsurvplot(fit, data = clin1, conf.int = F, # 增加置信区间
           risk.table = TRUE, pval=T) #不同时间点风险人数表，增加P值

res.cox <- coxph(Surv(OS,STATOS)~UICC,data =clin1)
summary(res.cox)
cox.zph(res.cox)


table(clin1$UICC,clin1$STATOS)
sum(diag(table(clin1$UICC,clin1$STATOS)))/sum(table(clin1$UICC,clin1$STATOS))










############################ All important nucleotides present in KEGG maps ###################
library(pathview)

####################### summary on all key metabolites
Nucleotide111 # 8 pan-cancer nucleotides
Nucleotide1111 # 25 NAC related nucleotides
Nucleotide11111 # 24 therapy response related nucleotides

sum111 <- unique(c(Nucleotide111,Nucleotide1111,Nucleotide11111))
length(sum111)

sum1111 <- data.frame(n1=rep(NA,length(sum111)),n2=rep(NA,length(sum111)),n3=rep(NA,length(sum111)))
rownames(sum1111) <- sum111
dim(sum1111)
head(sum1111)
sum1111$n1[match(Nucleotide111,rownames(sum1111))] <- rep(-1,length(match(Nucleotide111,rownames(sum1111))))
sum1111$n2[match(Nucleotide1111,rownames(sum1111))] <- rep(0,length(match(Nucleotide1111,rownames(sum1111))))
sum1111$n3[match(Nucleotide11111,rownames(sum1111))] <- rep(1,length(match(Nucleotide11111,rownames(sum1111))))
sum1111

pathview(cpd.data = sum1111,
         pathway.id = '00230',species = "hsa",out.suffix="summary of important purine metabolites",
         keys.align = "y", kegg.native = T,match.data = F, multi.state = T, same.layer = T,
         low=list(cpd="red"),mid = list(cpd="blue"),high = list(cpd="green"),
         discrete = list(cpd=T),bins = list(cpd=2))

pathview(cpd.data = sum1111,
         pathway.id = '00240',species = "hsa",out.suffix="summary of important pyrimidine metabolites",
         keys.align = "y", kegg.native = T,match.data = F, multi.state = T, same.layer = T,
         low=list(cpd="red"),mid = list(cpd="blue"),high = list(cpd="green"),
         discrete = list(cpd=T),bins = list(cpd=2))



################################ Venn plot for important metabolite
setwd("C:/MyRdata8/Summary")
head(nucleotide)

n51 <- read.csv("summary on single metabolite from multivariate cox analysis5.csv")
dim(n51)
head(n51)
unique(n51$compound_id)

n8 <- read.csv("8 pan-cacner metabolites.csv")
dim(n8)
head(n8)
unique(n8$compound_id)

n25 <- read.csv("25 NAC related metabolites.csv")
dim(n25)
head(n25)
unique(n25$compound_id)

n24 <- read.csv("24 therapy response related metabolites.csv")
dim(n24)
head(n24)
unique(n24$compound_id)

n67 <- read.csv("67 nucleotides for classifier in cancer chemotherapy.csv")
dim(n67)
head(n67)
n67.id <- nucleotide[match(n67$Nucleotide,nucleotide$name),6]
n67.id

n26 <- read.csv("26 nucleotides for classifier in advanced gastric cancer.csv")
head(n26)
dim(n26)
n26.id <- nucleotide[match(n26$Nucleotide,nucleotide$name),6]
n26.id

library(venn)
library(VennDiagram)

venn.list <- list(unique(n51$compound_id),unique(n8$compound_id),unique(n25$compound_id),
                  unique(n24$compound_id),n67.id,n26.id)
names(venn.list) <- c("n51","n8","n25","n24","n67","n26")
venn.list
venn(venn.list)


############################### separate mapping on KEGG
head(nucleotide)
purine.gene <- read.csv("purine gene list.csv")
purine.gene$Name
pyrimidine.gene <- read.csv("pyrimidine gene list.csv")
pyrimidine.gene






n25 <- read.csv("25 NAC related metabolites.csv")
dim(n25)
head(n25)
unique(n25$KEGG.ID)

sum.134 <- setdiff(nucleotide$KEGG.ID,unique(n25$KEGG.ID))
length(sum.134)

sum.159 <- c(rep(0,134),rep(1,25))
sum.159
names(sum.159) <- c(sum.134,unique(n25$KEGG.ID))
sum.159

pathview(gene.data=as.character(purine.gene$Name),cpd.data = sum.159,
         pathway.id = '00230',species = "hsa",out.suffix="purine.gene.cpd.25",
         keys.align = "y", kegg.native = T,key.pos = "topright",na.col = "#E5F5E0",
         limit=list(gene = 1,cpd=1),mid = list(cpd="white"),
         high = list(gene="#E5F5E0",cpd="green"),
         discrete = list(gene=T,cpd=T),
         bins = list(gene=1,cpd=2))

pathview(gene.data=as.character(pyrimidine.gene$Name),cpd.data = sum.159,
         pathway.id = '00240',species = "hsa",out.suffix="pyrimidine.gene.cpd.25",
         keys.align = "y", kegg.native = T,key.pos = "topright",na.col = "#E5F5E0",
         limit=list(gene = 1,cpd=1),mid = list(cpd="white"),
         high = list(gene="#E5F5E0",cpd="green"),
         discrete = list(gene=T,cpd=T),
         bins = list(gene=1,cpd=2))






n24 <- read.csv("24 therapy response related metabolites.csv")
dim(n24)
head(n24)
unique(n24$KEGG.ID)

sum.135 <- setdiff(nucleotide$KEGG.ID,unique(n24$KEGG.ID))
length(sum.135)

sum.159 <- c(rep(0,135),rep(1,24))
sum.159
names(sum.159) <- c(sum.135,unique(n24$KEGG.ID))
sum.159

pathview(gene.data=as.character(purine.gene$Name),cpd.data = sum.159,
         pathway.id = '00230',species = "hsa",out.suffix="purine.gene.cpd.24",
         keys.align = "y", kegg.native = T,key.pos = "topright",na.col = "#E5F5E0",
         limit=list(gene = 1,cpd=1),mid = list(cpd="white"),
         high = list(gene="#E5F5E0",cpd="green"),
         discrete = list(gene=T,cpd=T),
         bins = list(gene=1,cpd=2))

pathview(gene.data=as.character(pyrimidine.gene$Name),cpd.data = sum.159,
         pathway.id = '00240',species = "hsa",out.suffix="pyrimidine.gene.cpd.24",
         keys.align = "y", kegg.native = T,key.pos = "topright",na.col = "#E5F5E0",
         limit=list(gene = 1,cpd=1),mid = list(cpd="white"),
         high = list(gene="#E5F5E0",cpd="green"),
         discrete = list(gene=T,cpd=T),
         bins = list(gene=1,cpd=2))







n67 <- read.csv("67 nucleotides for classifier in cancer chemotherapy.csv")
dim(n67)
head(n67)
n67.id <- nucleotide[match(n67$Nucleotide,nucleotide$name),1]
n67.id

sum.92 <- setdiff(nucleotide$KEGG.ID,n67.id)
length(sum.92)

sum.159 <- c(rep(0,92),rep(1,67))
sum.159
names(sum.159) <- c(sum.92,n67.id)
sum.159

pathview(gene.data=as.character(purine.gene$Name),cpd.data = sum.159,
         pathway.id = '00230',species = "hsa",out.suffix="purine.gene.cpd.67",
         keys.align = "y", kegg.native = T,key.pos = "topright",na.col = "#E5F5E0",
         limit=list(gene = 1,cpd=1),mid = list(cpd="white"),
         high = list(gene="#E5F5E0",cpd="green"),
         discrete = list(gene=T,cpd=T),
         bins = list(gene=1,cpd=2))

pathview(gene.data=as.character(pyrimidine.gene$Name),cpd.data = sum.159,
         pathway.id = '00240',species = "hsa",out.suffix="pyrimidine.gene.cpd.67",
         keys.align = "y", kegg.native = T,key.pos = "topright",na.col = "#E5F5E0",
         limit=list(gene = 1,cpd=1),mid = list(cpd="white"),
         high = list(gene="#E5F5E0",cpd="green"),
         discrete = list(gene=T,cpd=T),
         bins = list(gene=1,cpd=2))






n26 <- read.csv("26 nucleotides for classifier in advanced gastric cancer.csv")
head(n26)
dim(n26)
n26.id <- nucleotide[match(n26$Nucleotide,nucleotide$name),1]
n26.id

sum.133 <- setdiff(nucleotide$KEGG.ID,n26.id)
length(sum.133)

sum.159 <- c(rep(0,133),rep(1,26))
sum.159
names(sum.159) <- c(sum.133,n26.id)
sum.159

pathview(gene.data=as.character(purine.gene$Name),cpd.data = sum.159,
         pathway.id = '00230',species = "hsa",out.suffix="purine.gene.cpd.26",
         keys.align = "y", kegg.native = T,key.pos = "topright",na.col = "#E5F5E0",
         limit=list(gene = 1,cpd=1),mid = list(cpd="white"),
         high = list(gene="#E5F5E0",cpd="green"),
         discrete = list(gene=T,cpd=T),
         bins = list(gene=1,cpd=2))

pathview(gene.data=as.character(pyrimidine.gene$Name),cpd.data = sum.159,
         pathway.id = '00240',species = "hsa",out.suffix="pyrimidine.gene.cpd.26",
         keys.align = "y", kegg.native = T,key.pos = "topright",na.col = "#E5F5E0",
         limit=list(gene = 1,cpd=1),mid = list(cpd="white"),
         high = list(gene="#E5F5E0",cpd="green"),
         discrete = list(gene=T,cpd=T),
         bins = list(gene=1,cpd=2))



######################### nucleotide associated drug pathway #########################

####################################################### SMPDB pathway
setwd("C:/MyRdata8/Summary/SMPDB/smpdb_pathways.csv")

smpdb_pathways <- read.csv("smpdb_pathways.csv")
dim(smpdb_pathways)
head(smpdb_pathways)
factor(smpdb_pathways$Subject)

drug <- smpdb_pathways[smpdb_pathways$Subject=="Drug Action"|smpdb_pathways$Subject=="Drug Metabolism",]
dim(drug)
head(drug)


### all drugs related pathways and metabolites
setwd("C:/MyRdata8/Summary/SMPDB/smpdb_metabolites.csv")

drug1 <- read.csv(paste(drug$SMPDB.ID[1],"_metabolites.csv",sep = ''))
for (i in 2:length(drug$SMPDB.ID)) {
  drug2 <- read.csv(paste(drug$SMPDB.ID[i],"_metabolites.csv",sep = ''))
  drug1 <- rbind(drug1,drug2)
}

dim(drug1)
head(drug1)


### Nucleotide related drug pathway
setwd("C:/MyRdata8/Summary")

nucleotide2 <- read.csv("Nucleotide metabolism1.csv")
dim(nucleotide2)
head(nucleotide2)
nucleotide2 <- nucleotide2[-which(nucleotide2$HMDB.ID==""),]
dim(nucleotide2)
head(nucleotide2)

nucleotide.drug1 <- data.frame()
for (i in 1:131) {
  nucleotide.drug <- drug1[which(drug1$HMDB.ID==nucleotide2$HMDB.ID[i]),]
  nucleotide.drug1 <- rbind(nucleotide.drug1,nucleotide.drug)
}
dim(nucleotide.drug1)
head(nucleotide.drug1)

nucleotide.drug2 <- nucleotide.drug1[,c(2,3,6,7)]
head(nucleotide)
nucleotide.drug2$name <- nucleotide[match(nucleotide.drug2$HMDB.ID,nucleotide$compound_id),5]
dim(nucleotide.drug2)
head(nucleotide.drug2)

drug1 <- c()
for (i in 1:1021) {
  drug2 <- strsplit(nucleotide.drug2$Pathway.Name,split=" ")[[i]][1]
  drug1 <- c(drug1,drug2)
}
head(drug1)
length(drug1)

# nucleotide associated drug
nucleotide.drug2$drug <- drug1
dim(nucleotide.drug2)
head(nucleotide.drug2)

# 105 detectable nucleotide associated drug
sum2 <- read.csv("annotated and detectable 105 nucleotides in all cancers.csv")
sum3 <- sum2[,-1]
head(sum3)
dim(sum3)
sum.105 <- data.frame(HMDB.ID=unique(sum3$compound_id))
dim(sum.105)
head(sum.105)

nucleotide.drug3 <- merge(nucleotide.drug2,sum.105)
dim(nucleotide.drug3)
head(nucleotide.drug3)
write.csv(nucleotide.drug3,"105 nucleotides associated drug pathway.csv")


# 105 detectable nucleotide associated chemotherapy drug
dim(drug)
head(drug)

drug3 <- drug[,c(3,5)]
colnames(drug3)[1] <- "Pathway.Name"

nucleotide.drug4 <- merge(nucleotide.drug3,drug3)
dim(nucleotide.drug4)
head(nucleotide.drug4)

i=100
index <- c()
for (i in 1:nrow(nucleotide.drug4)) {
  if (("cancer" %in% strsplit(nucleotide.drug4$Description[i],' ')[[1]])|
      ("cancers" %in% strsplit(nucleotide.drug4$Description[i],' ')[[1]])|
      ("anticancer" %in% strsplit(nucleotide.drug4$Description[i],' ')[[1]])|
      ("anti-cancer" %in% strsplit(nucleotide.drug4$Description[i],' ')[[1]])|
      ("chemotherapy" %in% strsplit(nucleotide.drug4$Description[i],' ')[[1]])|
      ("chemotherapeutic" %in% strsplit(nucleotide.drug4$Description[i],' ')[[1]])|
      ("chemotherapeutics" %in% strsplit(nucleotide.drug4$Description[i],' ')[[1]])
  ) {
    index <- c(index,i)
  }
}

drug.cancer <- nucleotide.drug4[index,] 
dim(drug.cancer)
head(drug.cancer)
unique(drug.cancer$drug)
unique(drug.cancer$name)

write.csv(drug.cancer,"105 nucleotides associated chemotherapy drug pathway.csv")





################################# nucleotide associated drug
setwd("C:/MyRdata8/Summary")
library(reshape)
nucleotide14 <- read.csv("nucleotide classification.csv")
head(nucleotide14)

ress <- data.frame()
for (i in 1:ncol(nucleotide14)) {
  ress1 <- data.frame(pathway=rep(colnames(nucleotide14)[i],length(nucleotide14[,i])),kegg=nucleotide14[,i])
  ress <- rbind(ress1,ress)
}
dim(ress)
head(ress)
ress2 <- ress[-which(ress$kegg==""),]
dim(ress2)
head(ress2)

ress5 <- read.csv("nucleotide classification1.csv")
dim(ress5)
head(ress5)

ress6 <- merge(ress5,ress2)
dim(ress6)
head(ress6)

nucleotide2 <- read.csv("Nucleotide metabolism1.csv")
dim(nucleotide2)
head(nucleotide2)

ress6$name <- nucleotide2[match(ress6$kegg,nucleotide2$KEGG.ID),5]
dim(ress6)
head(ress6)
write.csv(ress6,"nucleotide-metabolic module.csv")

drug <- read.csv("drugable nucleotide2.csv")
dim(drug)
head(drug)
colnames(drug)[1] <- "name"
head(drug)
unique(drug$name)

ress8 <- merge(ress6,drug,all.y = T)
dim(ress8)
head(ress8)

ress9 <- ress8[,-c(4)]
head(ress9)

dim(nucleotide2)
head(nucleotide2)
ress9$Pathway <- nucleotide2[match(ress9$name,nucleotide2$name1),2]
dim(ress9)

ress10 <- ress9[,c(3,2,1,4,5,6,7,8)]
dim(ress10)
head(ress10)

unique(ress10$name)
write.csv(ress10,"24 nucleotides associated chemotherapy drug pathway.csv")


######################## 11 pan-cancer nucleotides (figure 4) involved pathway/module ####################

pathway <- read.csv("nucleotide-metabolic module.csv")
pathway <- pathway[,-1]
colnames(pathway)[1] <- "module"
dim(pathway)
head(pathway)
intersect(pathway$name,n11$name)

nucleotide2 <- read.csv("Nucleotide metabolism1.csv")
dim(nucleotide2)
head(nucleotide2)

n11$pathway <- nucleotide2[match(n11$name,nucleotide2$name1),2]
n11 <- n11[,c(6,7)]
head(n11)

n11.1 <- data.frame(name=n11$name)
n11.1

pathway1 <- pathway[,c(1,4)]
n11.2 <- merge(n11.1,pathway1)
colnames(n11.2)[2] <- "pathway"
n11.2

n11.3 <- rbind(n11,n11.2)
n11.3$value <- rep("1",nrow(n11.3))
dim(n11.3)
head(n11.3)

unique(n11.3$pathway)
unique(n11.3$name)


library(ggplot2)
n11.3$name <- factor(n11.3$name,levels = rev(unique(n11.3$name)))
n11.3$pathway <- factor(n11.3$pathway,levels = unique(n11.3$pathway))

ggplot(n11.3,aes(x=pathway,y=name))+
  geom_tile(aes(fill=value),color="black")+
  scale_fill_manual(values = c("orange"))+
  xlab("Nucleotide")+ylab("Cancer")+guides(fill=F)+
  theme(axis.text.x = element_text(size = 8,vjust=1,hjust=1,angle=45),
        axis.text.y=element_text(size = 8))


######################## 11 pan-cancer nucleotides (figure 4) involved drug correlation network ####################

drug4 <- read.csv("24 nucleotides associated chemotherapy drug pathway1.csv")
dim(drug4)
head(drug4)
drug5 <- drug4[,c(4:6)]
dim(drug5)
head(drug5)

dim(n11.3)
head(n11.3)

n11.4 <- merge(n11.3,drug5,all.x=T)
head(n11.4)

library(ggplot2)
n11.4$name <- factor(n11.4$name,levels = rev(unique(n11.3$name)))

ggplot(n11.4,aes(x=Drug,y=name))+
  geom_tile(aes(fill=Correlation.type),color="black")+
  scale_fill_manual(values = c("purple","light blue","gold","cyan"),
                    na.translate=T)+
  xlab("Nucleotide")+ylab("Cancer")+guides(fill=F)+
  theme(axis.text.x = element_text(size = 8,vjust=1,hjust=1,angle=45),
        axis.text.y=element_text(size = 8))

head(colors(),300)
######################## 20 nucleotides (figure 5) involved pathway/module ####################
dim(n21)
head(n21)

pathway <- read.csv("nucleotide-metabolic module.csv")
pathway <- pathway[,-1]
colnames(pathway)[1] <- "module"
dim(pathway)
head(pathway)
intersect(pathway$name,n21$Feature)

nucleotide2 <- read.csv("Nucleotide metabolism1.csv")
dim(nucleotide2)
head(nucleotide2)

n21$pathway <- nucleotide2[match(n21$Feature,nucleotide2$name1),2]
n21 <- n21[,c(1,6)]
colnames(n21)[1] <- "name"
head(n21)

n21.1 <- data.frame(name=n21$name)
n21.1

pathway1 <- pathway[,c(1,4)]
n21.2 <- merge(n21.1,pathway1)
colnames(n21.2)[2] <- "pathway"
n21.2

n21.3 <- rbind(n21,n21.2)
n21.3$value <- rep("1",nrow(n21.3))
dim(n21.3)
head(n21.3)

unique(n21.3$pathway)
unique(n21.3$name)


library(ggplot2)
n21.3$name <- factor(n21.3$name,levels = rev(unique(n21.3$name)))
n21.3$pathway <- factor(n21.3$pathway,levels = unique(n21.3$pathway))

ggplot(n21.3,aes(x=pathway,y=name))+
  geom_tile(aes(fill=value),color="black")+
  scale_fill_manual(values = c("orange"))+
  xlab("Nucleotide")+ylab("Cancer")+guides(fill=F)+
  theme(axis.text.x = element_text(size = 7,vjust=1,hjust=1,angle=45),
        axis.text.y=element_text(size = 8))




######################## 20 pan-cancer nucleotides (figure 5A) involved drug correlation network ####################

drug4 <- read.csv("24 nucleotides associated chemotherapy drug pathway1.csv")
dim(drug4)
head(drug4)
drug5 <- drug4[,c(4:6)]
dim(drug5)
head(drug5)

dim(n20)
head(n20)
n20.1 <- n20
colnames(n20.1)[1] <- "name"
head(n20.1)

n20.2 <- merge(n20.1,drug5,all.x=T)
head(n20.2)
n20.3 <- n20.2[,c(1,6,7)]
dim(n20.3)
head(n20.3)

library(ggplot2)
n20.3$name <- factor(n20.3$name,levels = rev(unique(n20.1$name)))
head(n20.3)

ggplot(n20.3,aes(x=Drug,y=name))+
  geom_tile(aes(fill=Correlation.type),color="black")+
  scale_fill_manual(values = c("purple","light blue","gold","cyan"),
                    na.translate=T)+
  xlab("Nucleotide")+ylab("Cancer")+guides(fill=F)+
  theme(axis.text.x = element_text(size = 8,vjust=1,hjust=1,angle=45),
        axis.text.y=element_text(size = 8))

head(colors(),300)
######################## 49 nucleotides (figure 5) involved pathway/module ####################
dim(n49)
head(n49)

pathway <- read.csv("nucleotide-metabolic module.csv")
pathway <- pathway[,-1]
colnames(pathway)[1] <- "module"
dim(pathway)
head(pathway)
intersect(pathway$name,n49$Feature)

nucleotide2 <- read.csv("Nucleotide metabolism1.csv")
dim(nucleotide2)
head(nucleotide2)

n49$pathway <- nucleotide2[match(n49$Feature,nucleotide2$name1),2]
n49 <- n49[,c(1,6)]
colnames(n49)[1] <- "name"
head(n49)
dim(n49)

n49.1 <- data.frame(name=n49$name)
n49.1

pathway1 <- pathway[,c(1,4)]
n49.2 <- merge(n49.1,pathway1)
colnames(n49.2)[2] <- "pathway"
n49.2

n49.3 <- rbind(n49,n49.2)
n49.3$value <- rep("1",nrow(n49.3))
dim(n49.3)
head(n49.3)

unique(n49.3$pathway)
unique(n49.3$name)


library(ggplot2)
n49.3$name <- factor(n49.3$name,levels = rev(unique(n49.3$name)))
n49.3$pathway <- factor(n49.3$pathway,levels = unique(n49.3$pathway))

ggplot(n49.3,aes(x=pathway,y=name))+
  geom_tile(aes(fill=value),color="black")+
  scale_fill_manual(values = c("orange"))+
  xlab("Nucleotide")+ylab("Cancer")+guides(fill=F)+
  theme(axis.text.x = element_text(size = 7,vjust=1,hjust=1,angle=45),
        axis.text.y=element_text(size = 8))




######################## 49 pan-cancer nucleotides (figure 5B) involved drug correlation network ####################

drug4 <- read.csv("24 nucleotides associated chemotherapy drug pathway1.csv")
dim(drug4)
head(drug4)
drug5 <- drug4[,c(4:6)]
dim(drug5)
head(drug5)

dim(n49)
head(n49)
n49.1 <- n49
colnames(n49.1)[1] <- "name"
head(n49.1)

n49.2 <- merge(n49.1,drug5,all.x=T)
head(n49.2)
n49.3 <- n49.2[,c(1,6,7)]
dim(n49.3)
head(n49.3)

library(ggplot2)
n49.3$name <- factor(n49.3$name,levels = rev(unique(n49.1$name)))
head(n49.3)

ggplot(n49.3,aes(x=Drug,y=name))+
  geom_tile(aes(fill=Correlation.type),color="black")+
  scale_fill_manual(values = c("purple","light blue","gold","cyan"),
                    na.translate=T)+
  xlab("Nucleotide")+ylab("Cancer")+guides(fill=F)+
  theme(axis.text.x = element_text(size = 7,vjust=1,hjust=1,angle=45),
        axis.text.y=element_text(size = 7))

######################## 53 nucleotides (figure 6) involved pathway/module ####################
head(n53)
n53 <- data.frame(name=n53)
dim(n53)
head(n53)

setwd("C:/MyRdata8/Summary")
pathway <- read.csv("nucleotide-metabolic module.csv")
pathway <- pathway[,-1]
colnames(pathway)[1] <- "module"
dim(pathway)
head(pathway)
intersect(pathway$name,n53$name)

nucleotide2 <- read.csv("Nucleotide metabolism1.csv")
dim(nucleotide2)
head(nucleotide2)

n53$pathway <- nucleotide2[match(n53$name,nucleotide2$name1),2]
head(n53)
dim(n53)

n53.1 <- data.frame(name=n53$name)
n53.1

head(pathway)
pathway1 <- pathway[,c(1,4)]
head(pathway1)
n53.2 <- merge(n53.1,pathway1)
colnames(n53.2)[2] <- "pathway"
n53.2

n53.3 <- rbind(n53,n53.2)
n53.3$value <- rep("1",nrow(n53.3))
dim(n53.3)
head(n53.3)

unique(n53.3$pathway)
unique(n53.3$name)


library(ggplot2)
n53.3$name <- factor(n53.3$name,levels = rev(unique(n53.3$name)))
n53.3$pathway <- factor(n53.3$pathway,levels = unique(n53.3$pathway)[c(1,2,8,3,9,4,5,6,7,13,14,11,12,16,10,15,17)])

ggplot(n53.3,aes(x=pathway,y=name))+
  geom_tile(aes(fill=value),color="black")+
  scale_fill_manual(values = c("orange"))+
  xlab("Nucleotide")+ylab("Cancer")+guides(fill=F)+
  theme(axis.text.x = element_text(size = 7,vjust=1,hjust=1,angle=45),
        axis.text.y=element_text(size = 7))






library(reshape2)
dim(n53.3)
head(n53.3)
n53.3 <- n53.3[,-3]

n53.4 <- dcast(n53.3,pathway~name)
dim(n53.4)
n53.5 <- melt(n53.4)
dim(n53.5)
head(n53.5)
n53.5$value <- ifelse(n53.5$value>0,1,0)

n53.5$variable <- factor(n53.5$variable,levels = unique(n53.3$name))
n53.5$pathway <- factor(n53.5$pathway,levels = rev(unique(n53.3$pathway)))
n53.5$value <- as.factor(n53.5$value)

ggplot(n53.5,aes(x=variable,y=pathway))+
  geom_tile(aes(fill=value),color="black")+
  scale_fill_manual(values = c("white","orange"))+
  xlab("Nucleotide")+ylab("Cancer")+guides(fill=F)+
  theme(axis.text.x = element_text(size = 8,vjust=1,hjust=1,angle=90),
        axis.text.y=element_text(size = 8))

######################## 53 nucleotides (figure 6) involved drug correlation network ####################
setwd("C:/MyRdata8/Summary")
drug4 <- read.csv("24 nucleotides associated chemotherapy drug pathway1.csv")
dim(drug4)
head(drug4)
drug5 <- drug4[,c(4:6)]
dim(drug5)
head(drug5)

dim(n53)
head(n53)
n53.1 <- data.frame(name=n53$name)
n53.1

n53.2 <- merge(n53.1,drug5,all.x=T)
head(n53.2)
n53.3 <- n53.2
dim(n53.3)
head(n53.3)

library(ggplot2)
n53.3$name <- factor(n53.3$name,levels = rev(unique(n53.1$name)))
head(n53.3)

ggplot(n53.3,aes(x=Drug,y=name))+
  geom_tile(aes(fill=Correlation.type),color="black")+
  scale_fill_manual(values = c("purple","light blue","gold","cyan"),
                    na.translate=T)+
  xlab("Nucleotide")+ylab("Cancer")+guides(fill=F)+
  theme(axis.text.x = element_text(size = 7,vjust=1,hjust=1,angle=45),
        axis.text.y=element_text(size = 7))

##################### key nucleotides by overlaping ######################
setwd("C:/MyRdata8/Summary")

n8 <- read.csv("8 pan-cancer nucleotides.csv")
n8 <- n8$x
n8
#write.csv(n8,"8 pan-cancer nucleotides.csv")

n20 <- read.csv("20 NAC-related nucleotides.csv")
head(n20)
n20.1 <- n20$x
n20.1
#write.csv(n20.1,"20 NAC-related nucleotides.csv")

n49 <- read.csv("49 therapy response-related nucleotides.csv")
head(n49)
n49.1 <- n49$x
n49.1
#write.csv(n49.1,"49 therapy response-related nucleotides.csv")

n53 <- read.csv("53 pan-cancer nucleotides predict chemotherapy response.csv")
head(n53)
n53.1 <- n53$x
n53.1
#write.csv(n53,"53 pan-cancer nucleotides predict chemotherapy response.csv")

n41 <- read.csv("41 nucleotides in RF classifier For GC.csv")
n41.1 <- n41$x
n41.1

intersect(intersect(intersect(intersect(n8,n21.1),n49.1),n53.1),n41.1)



