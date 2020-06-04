## Install packages and load data
install.packages("partykit", dependencies = TRUE)
install.packages("rpart.plot")
install.packages("survminer", repo = 'https://mac.R-project.org')
install.packages("pbkrtest")
install.packages("readxl")
install.packages("ggsurvplot")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.11")
BiocManager::install("survcomp")

library(survival)
library(survminer)
library(ggplot2)
library(rpart)
library(rpart.plot)
library(readxl)
library(tidyverse)
library(dplyr)
library(partykit)
library(survcomp)


# Load in SEER 18 practice data
SEER <- read_excel("SEER18_allvar.xlsx")

# coerce gender into a factor
SEER$Sex<-as.factor(SEER$Sex)

# turn age column into numeric, removing 'years' and '+' from dataset 
SEER$Age=as.numeric(gsub("[a-zA-Z+]", "", SEER$Age)) 

## Drop rows that have blank(s)/unknown/Prostate values for StageA variable
SEER = SEER %>% filter(StageA !="Localized/regional (Prostate cases)") %>% filter(StageA !="Blank(s)") %>% filter(StageA !="Unstaged")
SEER$StageA=as.factor(SEER$StageA)
SEER$StageA=relevel(SEER$StageA,ref="Localized")

## Histology Recode ICDO3/WHO 2008 (remove clear cell carcinoma cases, 9044 code)
SEER$Histology<-factor(SEER$Histology, levels = c("9040","9041","9042","9043","9044"), labels = c("NOS","Monophasic","Monophasic","Biphasic","Clear cell"))
SEER=SEER[SEER$Histology != "Clear cell",]
SEER$Histology=relevel(SEER$Histology, ref = "NOS")

## Grade (Unknown = NaN)
SEER$Grade[SEER$Grade == "Unknown"] <- NaN
SEER$Grade<-factor(SEER$Grade, levels = c("Undifferentiated; anaplastic; Grade IV","Poorly differentiated; Grade III","Moderately differentiated; Grade II","Well differentiated; Grade I"), labels = c("IV","III","II","I"))

## Turn Race into just (White, Black, Other)
for (i in 1:nrow(SEER)){
  if (SEER$Race[i]=="White"){SEER$Race[i]="White"}
  else if (SEER$Race[i]=="Black"){SEER$Race[i]="Black"}
  else{SEER$Race[i]= "Other"}}
SEER$Race = as.factor(SEER$Race)
SEER$Race=relevel(SEER$Race, ref = "White")

## Generate Tumor size column (Size 1 = TS Summary 2016, Size 2 = CS Tumor Size 2004-2015, Size3 = EOD10)
TumorSize=array()
sizes = seq(1,998,1)
for (i in 1:nrow(SEER)){
  if (SEER$Size1[i] %in% sizes){
    TumorSize[i] = as.numeric(SEER$Size1[i])
  }
  else if (SEER$Size2[i] %in% sizes){
    TumorSize[i] = as.numeric(SEER$Size2[i])
  }
  else if (SEER$Size3[i] %in% sizes){
    TumorSize[i] = as.numeric(SEER$Size3[i])
  }
  else{TumorSize[i]=NaN} ## cases 000, 999, Blank(s) are missing
}

## Tumor size (1662 remaining cases after removing 333 NaN and cases >400mm for TS)
SEER$TumorSize=TumorSize 
SEER = SEER %>% filter(!is.na(TumorSize)) %>% filter(TumorSize<400)


## Create train and test sets (80/20 split)
set.seed(12345)
train <- SEER %>% dplyr::sample_frac(0.8)
test  <- dplyr::anti_join(SEER, train, by = "Patient ID")

## create survival object and KM plot for training data (overall, 1 strata)
S = Surv(time=train$Time,event=train$Status) # survival object for test data
fit1 <- survfit(S ~ 1, data = train)
ggsurv<-ggsurvplot(fit1, data = train, pval = TRUE, risk.table = TRUE)
ggarrange(ggsurv$plot, ggsurv$table, heights = c(2, 0.7),
          ncol = 1, nrow = 2, align = "v")

## Rpart Decision Tree with maximum depth (cp=0)
set.seed(20)
tfit_all = rpart(formula = S ~ Age+StageA+Race+Histology+Sex+Grade+TumorSize, data = train, control=rpart.control(minsplit = 20,cp=0), xval=10)
printcp(tfit_all) ## use 1-SD rule to choose cp
plotcp(tfit_all,  minline = FALSE, upper=c("splits"))

## Prune decision tree using to extract min_cp and cp_opt(1+SE rule)
min_cp <- tfit_all$cptable[which.min(tfit_all$cptable[,"xerror"]),"CP"]
min_xerr<-tfit_all$cptable[which.min(tfit_all$cptable[,"xerror"]),"xerror"]
xstd <-tfit_all$cptable[which.min(tfit_all$cptable[,"xerror"]),"xstd"] 
xerror<-min_xerr+xstd
cp_opt=max(tfit_all$cptable[tfit_all$cptable[,"xerror"]<xerror,"CP"])


tfit= prune(tfit_all, cp=cp) ## prune using cp_opt, can also use min_cp
rpart.plot(tfit)
tfit_2 <- as.party(tfit) ## coerce into party object, plot w/ associated KM curves 
plot(tfit_2)


## Conditional tree using partykit (???) much bigger tree
ctree <- ctree(S ~ Age+StageA+Race+Histology+Sex+Grade+TumorSize, data = train, control = ctree_control(testtype = "Bonferroni", alpha = 0.05, minbucket=30))
Groups<- predict(ctree, newdata = test, type="node") # TN 
rpart.plot(ctree)
plot(ctree)


## Rename terminal nodes as groups 1-6, w/6 having highest risk
## Group 1 =  node 3
## group 2 = node 6
## Group 3 = node 7
## Group 4 = node 8
## Group 5 =  node 10
## Group 6 = node 11

## Put train and test subjects into groups, using decision tree and relabel groups 1-6
train$Group = tfit$where 
test$Group <- predict(tfit_2, newdata = test, type = "node" )
train$Group<-factor(train$Group, levels = c("3","6","7","8","10","11"), labels = c("1","2","3","4","5", "6"))
test$Group<-factor(test$Group, levels = c("3","6","7","8","10","11"), labels = c("1","2","3","4","5", "6"))

## KM curve using group on training data (should be excellent fit)
S_train = Surv(train$Time, train$Status) 
fit_group <- survfit(S_train ~ Group, data = train)
ggsurv_group<-ggsurvplot(fit_group, data = train, pval = TRUE, title="KM training cohort", risk.table = TRUE,conf.int = TRUE)
ggarrange(ggsurv_group$plot, ggsurv_group$table, heights = c(2, 0.7),
          ncol = 1, nrow = 2, align = "v")


## Rpart model predictions on test data (individual risk relative to root node)
predtr <- predict(tfit, newdata = train)
predts <- predict(tfit, newdata = test)

## Concordance index for predts and observed risk 
perf <- concordance.index(x = predts, surv.time = test$Time, surv.event = test$Status, method = "noether",na.rm = TRUE)
print(perf[1:5])

## D-index
perf <- D.index(x = predts, surv.time = test$Time, surv.event = test$Status,na.rm = TRUE)
print(perf[1:6])

## ROC curve for binary 5 yr survival period 
perf <- tdrocc(x = predts, surv.time = test$Time,surv.event = test$Status, time = 5 * 365, na.rm = TRUE)
plot(x = 1 - perf$spec, y = perf$sens, type = "l", 
     xlab = "1 - specificity", ylab = "sensitivity", 
     xlim = c(0,1), ylim = c(0, 1), 
     main = "ROC curve\nat 5 years Cox PH")
lines(x = c(0, 1), y = c(0, 1), lty = 3, col = "red")

##Brier score
ddtr <- cbind(time = train$Time, event = train$Status,score = predtr)
ddts <- cbind(time = test$Time, event = test$Status,score = predts)
perf <- sbrier.score2proba(data.tr = data.frame(ddtr), data.ts = data.frame(ddts), method = "cox")
plot(x = perf$time, y = perf$bsc, xlab = "Time (days)",ylab = "Brier score", type = "l")
## null model
ddtr <- cbind(time = train$Time, event = train$Status,score = 1)
ddts <- cbind(time = test$Time, event = test$Status,score = 1)
perfnull <- sbrier.score2proba(data.tr = data.frame(ddtr),data.ts = data.frame(ddts), method = "prodlim")
lines(x = perfnull$time, y = perfnull$bsc, col = "red",lty = 2)
## legend
legend("bottomright", title = "Integrated Brier score", legend = c(sprintf("Rpart model = %.3g", perf$bsc.integrated),
                                                                   sprintf("Null model = %.3g", perfnull$bsc.integrated)),col = c("black", "red"), lty = c(1, 2))


## KM Plot for test data, stratified by group 
S_test = Surv(time=test$Time,event=test$Status) ## survival object for test data 
fit_test <- survfit(S_test ~ Group, data = test)
ggsurv_test<-ggsurvplot(fit_test, data = test, pval = TRUE,title="KM Test cohort (n=332)", risk.table = TRUE, conf.int = TRUE)
ggarrange(ggsurv_test$plot, ggsurv_test$table, heights = c(2, 0.7),
          ncol = 1, nrow = 2, align = "v")

## pairwise comparisons of group curves on test data using log-rank test (rho=0), or Peto Peto (rho = 1)
res<-pairwise_survdiff(S_test ~ Group,data=test, p.adjust.method = "bonferroni", rho =0)






## Cox PH formula (using group), not plotted
fitcox <- coxph(S ~ Group, data = train)
ggforest(fitcox, data = train)
summary(fitcox)
cox.zph(fitcox, transform = "km") # Schoenfeld's test for PH

