## GrimAge Longitudinal Trajectories in the Eighth Decade and Kaplan-Meier Survival Analyses 

## Installing Requisite Packages 

if(!require(ggplot2)){
  install.packages("ggplot2")
}

if(!require(survival)){
  install.packages("survival")
}

if(!require(survminer)){
  install.packages("survminer")
}


if(!require(foreign)){
  install.packages("foreign")
}


library(ggplot2)
library(survival) 
library(survminer) 
library(foreign) 


## Read in Files

setwd('/Cluster_Filespace/Marioni_Group/Rob/GrimAge/')  

grim <- read.csv("GrimAge_LBC36_allwaves.csv")

## Longitudinal DNAm GrimAge Trajectories 

tar = read.csv("target_QC_age_sex_date.csv")
tar1 <- tar[tar$cohort=="LBC36",]
grim36 <- merge(grim, tar1[,c("Basename","ID")], by.x = "SampleID", by.y ="Basename")


pdf("GrimAge_vs_Age.pdf")
p <- ggplot(data = grim36, aes(x = Age, y = DNAmGrimAge, group = as.factor(ID)))
p + geom_line() + stat_smooth(aes(group = 1))  
dev.off()

grim36$grim_diff <- grim36$DNAmGrimAge - grim36$Age

p1 <- ggplot(data = grim36, aes(x = Age, y = grim_diff, group = ID))
p1 + geom_line() + stat_smooth(aes(group = 1))  



## Survival Analyses

tar1 <- tar[tar$cohort=="LBC36" & tar$WAVE==1,]
grim <- read.csv("GrimAge_output.csv")
grim36 <- merge(grim, tar1[,c("Basename","ID")], by.x="SampleID", by.y="Basename")

dead = read.spss("LBC1936_DNA_MethylationBasedEstimatorOfTeloLength_RM_28FEB2019.sav", to.data.frame=T)

dead$event = 0
dead$event[dead$dead=="DEAD"] <- 1
dead$age_event = ifelse(dead$event==0, dead$AgedaysApx_LastCensor/365.25, dead$agedays_death/365.25)

dead$LBC.subject.no = as.character(dead$lbc36no)
dead$age = dead$agedays_w1/365.25

d1 = merge(dead, grim36, by.x="LBC.subject.no", by.y="ID")

d1$tte <- d1$age_event - d1$age

d1$low_hi <- NA
d1$low_hi[d1$DNAmGrimAge <= quantile(d1$DNAmGrimAge, 0.25)] <- 0
d1$low_hi[d1$DNAmGrimAge >= quantile(d1$DNAmGrimAge, 0.75)] <- 1

fit <- survfit(Surv(tte, event) ~ low_hi,
               data = d1)
# Visualize with survminer
ggsurvplot(fit, data = d1)


## GrimAge - CoxPH model and plot

summary(coxph(Surv(tte, event) ~ age + factor(sex) + scale(DNAmGrimAge), data=d1))

d1$low_hi <- NA
d1$low_hi[d1$AgeAccelGrim <= quantile(d1$AgeAccelGrim, 0.25)] <- 0
d1$low_hi[d1$AgeAccelGrim >= quantile(d1$AgeAccelGrim, 0.75)] <- 1

fit <- survfit(Surv(tte, event) ~ low_hi,
               data = d1)
# Visualize with survminer
pdf("Survival_GrimAge.pdf")
ggsurvplot(fit, data = d1, legend.labs = c("Bottom Quartile", "Top Quartile"), xlab = "Years - Follow Up")
dev.off()


## ADM Adj Age - CoxPH model and plot

summary(coxph(Surv(tte, event) ~ age + factor(sex) + scale(DNAmADMAdjAge), data=d1))

d1$low_hi <- NA
d1$low_hi[d1$DNAmADMAdjAge <= quantile(d1$DNAmADMAdjAge, 0.25)] <- 0
d1$low_hi[d1$DNAmADMAdjAge >= quantile(d1$DNAmADMAdjAge, 0.75)] <- 1

fit <- survfit(Surv(tte, event) ~ low_hi,
               data = d1)
# Visualize with survminer
pdf("Survival_ADM_Adj_Age.pdf")
ggsurvplot(fit, data = d1, legend.labs = c("Bottom Quartile", "Top Quartile"), xlab = "Years - Follow Up", xlim = c(0,13))
dev.off()



## B2M Adj Age - CoxPH model and plot

summary(coxph(Surv(tte, event) ~ age + factor(sex) + scale(DNAmB2MAdjAge), data=d1))

d1$low_hi <- NA
d1$low_hi[d1$DNAmB2MAdjAge <= quantile(d1$DNAmB2MAdjAge, 0.25)] <- 0
d1$low_hi[d1$DNAmB2MAdjAge >= quantile(d1$DNAmB2MAdjAge, 0.75)] <- 1

fit <- survfit(Surv(tte, event) ~ low_hi,
               data = d1)
# Visualize with survminer
pdf("Survival_B2M_Adj_Age.pdf")
ggsurvplot(fit, data = d1, legend.labs = c("Bottom Quartile", "Top Quartile"), xlab = "Years - Follow Up")
dev.off()


## Cystatin C - CoxPH model and plot

summary(coxph(Surv(tte, event) ~ age.y + factor(sex.y) + scale(DNAmCystatinCAdjAge), data=d1))

d1$low_hi <- NA
d1$low_hi[d1$DNAmCystatinCAdjAge <= quantile(d1$DNAmCystatinCAdjAge, 0.25)] <- 0
d1$low_hi[d1$DNAmCystatinCAdjAge >= quantile(d1$DNAmCystatinCAdjAge, 0.75)] <- 1

fit <- survfit(Surv(tte, event) ~ low_hi,
               data = d1)
# Visualize with survminer
pdf("Survival_CystatinC_Adj_Age.pdf")
ggsurvplot(fit, data = d1, legend.labs = c("Bottom Quartile", "Top Quartile"), xlab = "Years - Follow Up")
dev.off()


## GDF15 - CoxPH model and plot

summary(coxph(Surv(tte, event) ~ age.y + factor(sex.y) + scale(DNAmGDF15AdjAge), data=d1))

d1$low_hi <- NA
d1$low_hi[d1$DNAmGDF15AdjAge <= quantile(d1$DNAmGDF15AdjAge, 0.25)] <- 0
d1$low_hi[d1$DNAmGDF15AdjAge >= quantile(d1$DNAmGDF15AdjAge, 0.75)] <- 1

fit <- survfit(Surv(tte, event) ~ low_hi,
               data = d1)
# Visualize with survminer
pdf("Survival_GDF15_Adj_Age.pdf")
ggsurvplot(fit, data = d1, legend.labs = c("Bottom Quartile", "Top Quartile"), xlab = "Years - Follow Up")
dev.off()


## Leptin - CoxPH model and plot

summary(coxph(Surv(tte, event) ~ age.y + factor(sex.y) + scale(DNAmLeptinAdjAge), data=d1))

d1$low_hi <- NA
d1$low_hi[d1$DNAmLeptinAdjAge <= quantile(d1$DNAmLeptinAdjAge, 0.25)] <- 0
d1$low_hi[d1$DNAmLeptinAdjAge >= quantile(d1$DNAmLeptinAdjAge, 0.75)] <- 1

fit <- survfit(Surv(tte, event) ~ low_hi,
               data = d1)
# Visualize with survminer
pdf("Survival_Leptin_Adj_Age.pdf")
ggsurvplot(fit, data = d1, legend.labs = c("Bottom Quartile", "Top Quartile"), xlab = "Years - Follow Up")
dev.off()



## PAI1 - CoxPH model and plot

summary(coxph(Surv(tte, event) ~ age.y + factor(sex.y) + scale(DNAmPAI1AdjAge), data=d1))

d1$low_hi <- NA
d1$low_hi[d1$DNAmPAI1AdjAge <= quantile(d1$DNAmPAI1AdjAge, 0.25)] <- 0
d1$low_hi[d1$DNAmPAI1AdjAge >= quantile(d1$DNAmPAI1AdjAge, 0.75)] <- 1

fit <- survfit(Surv(tte, event) ~ low_hi,
               data = d1)
# Visualize with survminer
pdf("Survival_PAI1_Adj_Age.pdf")
ggsurvplot(fit, data = d1, legend.labs = c("Bottom Quartile", "Top Quartile"), xlab = "Years - Follow Up")
dev.off()



## TIMP1 - CoxPH model and plot

summary(coxph(Surv(tte, event) ~ age.y + factor(sex.y) + scale(DNAmTIMP1AdjAge), data=d1))

d1$low_hi <- NA
d1$low_hi[d1$DNAmTIMP1AdjAge <= quantile(d1$DNAmTIMP1AdjAge, 0.25)] <- 0
d1$low_hi[d1$DNAmTIMP1AdjAge >= quantile(d1$DNAmTIMP1AdjAge, 0.75)] <- 1

fit <- survfit(Surv(tte, event) ~ low_hi,
               data = d1)
# Visualize with survminer
pdf("Survival_TIMP1_Adj_Age.pdf")
ggsurvplot(fit, data = d1, legend.labs = c("Bottom Quartile", "Top Quartile"), xlab = "Years - Follow Up")
dev.off()



## PACKYRS - CoxPH model and plot

summary(coxph(Surv(tte, event) ~ age.y + factor(sex.y) + scale(DNAmPACKYRSAdjAge), data=d1))

d1$low_hi <- NA
d1$low_hi[d1$DNAmPACKYRSAdjAge <= quantile(d1$DNAmPACKYRSAdjAge, 0.25)] <- 0
d1$low_hi[d1$DNAmPACKYRSAdjAge >= quantile(d1$DNAmPACKYRSAdjAge, 0.75)] <- 1

fit <- survfit(Surv(tte, event) ~ low_hi,
               data = d1)
# Visualize with survminer
pdf("Survival_PACKYRS_Adj_Age.pdf")
ggsurvplot(fit, data = d1, legend.labs = c("Bottom Quartile", "Top Quartile"), xlab = "Years - Follow Up")
dev.off()
