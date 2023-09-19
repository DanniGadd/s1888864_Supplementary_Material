###################################################################################

### LBC file integration 

###################################################################################

# Mean ages at each wave 

# Mean age W5 with SD
target <- read.csv("/Cluster_Filespace/Marioni_Group/LBC/LBC_methylation/archive/target_QC_age_sex_date.csv")
dat3 <- target %>% filter(WAVE == "5") # 801 people
dat3  <- dat3 %>% filter(cohort == "LBC36") # 801 people 

mean(dat3$age, na.rm = T)
sd(dat3$age, na.rm = T)

# Mean age W4 with SD
target <- read.csv("/Cluster_Filespace/Marioni_Group/LBC/LBC_methylation/archive/target_QC_age_sex_date.csv")
dat3 <- target %>% filter(WAVE == "4") # 801 people
dat3  <- dat3 %>% filter(cohort == "LBC36") # 801 people 

mean(dat3$age, na.rm = T)
sd(dat3$age, na.rm = T)

# Mean age W3 with SD
target <- read.csv("/Cluster_Filespace/Marioni_Group/LBC/LBC_methylation/archive/target_QC_age_sex_date.csv")
dat3 <- target %>% filter(WAVE == "3") # 801 people
dat3  <- dat3 %>% filter(cohort == "LBC36") # 801 people 

mean(dat3$age, na.rm = T)
sd(dat3$age, na.rm = T)

# Mean age W2 with SD
target <- read.csv("/Cluster_Filespace/Marioni_Group/LBC/LBC_methylation/archive/target_QC_age_sex_date.csv")
dat3 <- target %>% filter(WAVE == "2") # 801 people
dat3  <- dat3 %>% filter(cohort == "LBC36") # 801 people 

mean(dat3$age, na.rm = T)
sd(dat3$age, na.rm = T)

# Mean age W1 with SD
target <- read.csv("/Cluster_Filespace/Marioni_Group/LBC/LBC_methylation/archive/target_QC_age_sex_date.csv")
dat3 <- target %>% filter(WAVE == "1") # 801 people
dat3  <- dat3 %>% filter(cohort == "LBC36") # 801 people 

mean(dat3$age, na.rm = T)
sd(dat3$age, na.rm = T)

# Prep analyses for calculation of performance

screen

R

library(tidyverse)
library(foreign)

# Cognitive (g intercept and slope data) - from Simon
# Important to note that you won’t want to use the extracted slope data for folks that don’t have longitudinal data
cog <- read.csv("/Local_Scratch/Danni/GDFBNP/04_LBC/original_data/LBC1936_factors_of_curve_simon.csv")

# Serum markers and ICV, PVS counts, & stroke masks - from Paul 
# The key absences are the latest brain volumes (which I’ve recently received updates for), and dementia ascertainment
data <- read.spss("/Local_Scratch/Danni/GDFBNP/04_LBC/original_data/LBC1936_DNAmPredictorsOfSerumGDF15_and_NTproBNP_DG_23NOV2022.sav", to.data.frame=TRUE)


# Quantify what protein and imaging measurements we have available and from which waves

counts <- data.frame(Variable = 1:2, Present = 1:2, Absent = 1:2)

for(i in 1:length(colnames(data))){
  col_name <- as.character(colnames(data)[[i]])
  present <- table(is.na(data[,col_name]))[1]
  absent <- table(is.na(data[,col_name]))[2]

  counts[i,1] <- col_name
  counts[i,2] <- present
  counts[i,3] <- absent
}

counts[is.na(counts)] <- 0

# Quantify the cognitive data we have 

counts2 <- data.frame(Variable = 1:2, Present = 1:2, Absent = 1:2)

for(i in 1:length(colnames(cog))){
  col_name <- as.character(colnames(cog)[[i]])
  present <- table(is.na(cog[,col_name]))[1]
  absent <- table(is.na(cog[,col_name]))[2]

  counts2[i,1] <- col_name
  counts2[i,2] <- present
  counts2[i,3] <- absent
}

counts2[is.na(counts2)] <- 0
counts2 <- counts2[-1,]

# Save to check plan with simon/riccardo

write.csv(counts, '/Local_Scratch/Danni/GDFBNP/04_LBC/original_data/counts_paul_data.csv', row.names = F)
write.csv(counts2, '/Local_Scratch/Danni/GDFBNP/04_LBC/original_data/counts_simon_data.csv', row.names = F)


### Subset to protein specific data files and check how many have methylation at the specific waves

joint <- left_join(data, cog, by = 'lbc36no')

## GDF

# Get W2 LBC1936 DNAm to index 
target <- read.csv("/Cluster_Filespace/Marioni_Group/LBC/LBC_methylation/target_QC_age_sex_date.csv")
dat3 <- target %>% filter(WAVE == "2") # 801 people
dat3  <- dat3 %>% filter(cohort == "LBC36") # 801 people 
dat3 <- dat3[c(2,4,5)]

length(which(joint$lbc36no %in% dat3$ID)) # all 801 folks have DNAm 

GDF_w2 <- joint[complete.cases(joint$GDF15_W2),] # 798 individuals

length(which(GDF_w2$lbc36no %in% dat3$ID)) # 762 folks have DNAm and GDF at W2 

GDF_w2 <- GDF_w2[which(GDF_w2$lbc36no %in% dat3$ID),]
names(dat3)[2] <- 'lbc36no'
GDF_w2 <- left_join(GDF_w2, dat3, by = 'lbc36no')
# GDF_w2 <- GDF_w2[c(1,2,4,17,72)]
write.csv(GDF_w2, '/Local_Scratch/Danni/GDFBNP/04_LBC/LBC_testing/GDF_w2.csv', row.names = F)

# Get W4 LBC1936 DNAm to index 
target <- read.csv("/Cluster_Filespace/Marioni_Group/LBC/LBC_methylation/target_QC_age_sex_date.csv")
dat3 <- target %>% filter(WAVE == "4") # 589 people
dat3  <- dat3 %>% filter(cohort == "LBC36") # 507 people 
dat3 <- dat3[c(2,4,5)]

length(which(joint$lbc36no %in% dat3$ID)) # all 507 folks have DNAm 

GDF_w4 <- joint[complete.cases(joint$GDF15_W4),] # 331 individuals

length(which(GDF_w4$lbc36no %in% dat3$ID)) # 322 folks have DNAm and GDF at W4

GDF_w4 <- GDF_w4[which(GDF_w4$lbc36no %in% dat3$ID),]
names(dat3)[2] <- 'lbc36no'
GDF_w4 <- left_join(GDF_w4, dat3, by = 'lbc36no')
write.csv(GDF_w4, '/Local_Scratch/Danni/GDFBNP/04_LBC/LBC_testing/GDF_w4.csv', row.names = F)


# Are there common individuals measured at both waves?

length(which(GDF_w2$lbc36no %in% GDF_w4$lbc36no)) # 286 of those in W2 were also done at W4 


## BNP

# Get W3 LBC1936 DNAm to index 
target <- read.csv("/Cluster_Filespace/Marioni_Group/LBC/LBC_methylation/target_QC_age_sex_date.csv")
dat3 <- target %>% filter(WAVE == "3") 
dat3  <- dat3 %>% filter(cohort == "LBC36") # 619 people 
dat3 <- dat3[c(2,4,5)]

length(which(joint$lbc36no %in% dat3$ID)) # all 619 folks have DNAm 

BNP_w3 <- joint[complete.cases(joint$bld_NT_ProBNP_w3),] # 676 individuals

length(which(BNP_w3$lbc36no %in% dat3$ID)) # 616 folks have DNAm and GDF at W4

BNP_w3 <- BNP_w3[which(BNP_w3$lbc36no %in% dat3$ID),]
names(dat3)[2] <- 'lbc36no'
BNP_w3 <- left_join(BNP_w3, dat3, by = 'lbc36no')
write.csv(BNP_w3, '/Local_Scratch/Danni/GDFBNP/04_LBC/LBC_testing/BNP_w3.csv', row.names = F)


# Get W4 LBC1936 DNAm to index 
target <- read.csv("/Cluster_Filespace/Marioni_Group/LBC/LBC_methylation/target_QC_age_sex_date.csv")
dat3 <- target %>% filter(WAVE == "4") 
dat3  <- dat3 %>% filter(cohort == "LBC36") # 507 people 
dat3 <- dat3[c(2,4,5)]

length(which(joint$lbc36no %in% dat3$ID)) # all 507 folks have DNAm 

BNP_w4 <- joint[complete.cases(joint$bld_NT_ProBNP_w4),] # 515 individuals

length(which(BNP_w4$lbc36no %in% dat3$ID)) # 502 folks have DNAm and GDF at W4

BNP_w4 <- BNP_w4[which(BNP_w4$lbc36no %in% dat3$ID),]
names(dat3)[2] <- 'lbc36no'
BNP_w4 <- left_join(BNP_w4, dat3, by = 'lbc36no')
write.csv(BNP_w4, '/Local_Scratch/Danni/GDFBNP/04_LBC/LBC_testing/BNP_w4.csv', row.names = F)

# Are there common individuals measured at both waves?

length(which(BNP_w3$lbc36no %in% BNP_w4$lbc36no)) # 438 of those in W2 were also done at W4 

# # Get W5 LBC1936 DNAm to index 
# target <- read.csv("/Cluster_Filespace/Marioni_Group/LBC/LBC_methylation/target_QC_age_sex_date.csv")
# dat3 <- target %>% filter(WAVE == "5") 
# dat3  <- dat3 %>% filter(cohort == "LBC36") # 0 people - no DNAm



###################################################################################

### Test EpiScores in W4 (complete set)

###################################################################################

cd /Local_Scratch/Danni/GDFBNP/04_LBC/20k_GS_training/

screen

R

library(tidyverse)

W4 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_Analysis_20k/00_Updated_results/Results_scores/projections_LBC_W4_241122.csv")

BNP_w4 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_Analysis_20k/04_LBC/LBC_testing/BNP_w4.csv")
GDF_w4 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_Analysis_20k/04_LBC/LBC_testing/GDF_w4.csv")


### GDF W4
names(W4)[3] <- 'Basename'
join <- left_join(GDF_w4, W4, by = 'Basename')

library(bestNormalize)
for(i in colnames(join)[c(18)]){ 
  join[,i]<- orderNorm(join[,i])$x.t # Rank-Inverse Based Normaliation
}

cor.test(join$GDF15, join$GDF15_W4) # 0.36

join$age <- join$agedays_w4.x / 365.25
join$sex <- join$sex.x

null <- summary(lm(GDF15_W4 ~ age + sex, data=join))$r.squared
full <- summary(lm(GDF15_W4 ~ age + sex + join$GDF15, data=join))$r.squared
print(round(100*(full - null), 3))

# GDF15 450k
# LBC W4: 8.89%

# Add in PRS for GDF in LBC
prs <- read.table("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_PRS/lbc1936_genomewide_gdf.all.score", header = T)

names(prs)[2] <- 'lbc36no'
prs <- prs[-1]
names(prs)[2] <- 'GDFPRS'

join <- left_join(join, prs, by = 'lbc36no')

null <- summary(lm(GDF15_W4 ~ age + sex, data=join))$r.squared
null_prs <- summary(lm(GDF15_W4 ~ age + sex + GDFPRS, data=join))$r.squared
full <- summary(lm(GDF15_W4 ~ age + sex + join$GDF15, data=join))$r.squared
all <- summary(lm(GDF15_W4 ~ age + sex + join$GDF15 + GDFPRS, data=join))$r.squared
print(round(100*(full - null), 3))

# > null
# [1] 0.09482718
# > null_prs
# [1] 0.1240604
# > full
# [1] 0.1837803
# > all
# [1] 0.23148


table(is.na(join$GDF15))
table(is.na(join$GDF15_W4)) # 322 with both available 

library(ggpubr)
plot1 <- ggplot(join, aes(x=GDF15_W4, y=GDF15)) +
geom_point(colour = "tan1", size = 1) +
geom_smooth(method='lm', colour = "tan1") + 
theme(axis.text.x=element_text(size=20),     
      axis.text.y=element_text(size=20),
      axis.title.x=element_text(size=20),
      axis.title.y=element_text(size=20)) + 
stat_cor(aes(label = ..r.label..), method = "pearson", cor.coef.name = "r", size = 7.5) +
xlab("GDF15 protein") + ylab("GDF15 EpiScore") + 
theme_classic() # , label.x = 4.5, label.y = 4.4

# prot_model <- lm(Ig ~ age + sex + GDF15_W4, data = join)
# score_model <- lm(Ig ~ age + sex + GDF15.20k.with.450k.array, data = join)

# join <- join[which(join$enl_peri_space_hl_w2 %in% c('1','0')),]
# prot_model <- glm(enl_peri_space_hl_w2 ~ age + sex + GDF15_W4 + ICV_mm3_wX, 
#   family=binomial(link='logit'), data = join)
# score_model <- glm(enl_peri_space_hl_w2 ~ age + sex + GDF15.20k.with.450k.array + ICV_mm3_wX,
# family=binomial(link='logit'), data = join)


# Compare to GDF15 from grimage

LBC_grim <- read.csv('/Cluster_Filespace/Marioni_Group/LBC/LBC_methylation/LBC_clock_output_3489.csv')

target <- read.csv("/Cluster_Filespace/Marioni_Group/LBC/LBC_methylation/archive/target_QC_age_sex_date.csv")
dat3 <- target[which(target$WAVE %in% '4'),]
dat3  <- dat3 %>% filter(cohort == "LBC36") 
dat3 <- dat3 %>% select('Basename', 'set', 'ID')

which(LBC_grim$Basename %in% dat3$Basename)

join <- left_join(join, LBC_grim, by = 'Basename')

# null <- summary(lm(GDF15_W4 ~ age + sex, data=join))$r.squared
# full <- summary(lm(GDF15_W4 ~ age + sex + join$GDF15, data=join))$r.squared
# print(round(100*(full - null), 3))

null <- summary(lm(GDF15_W4 ~ age + sex, data=join))$r.squared
full <- summary(lm(GDF15_W4 ~ age + sex + join$DNAmGDF15, data=join))$r.squared
print(round(100*(full - null), 3))


### BNP w4
names(W4)[3] <- 'Basename'
join <- left_join(BNP_w4, W4, by = 'Basename')

join$bld_NT_ProBNP_w4 <- as.numeric(join$bld_NT_ProBNP_w4)

library(bestNormalize)
for(i in colnames(join)[c(20)]){ 
  join[,i]<- orderNorm(join[,i])$x.t # Rank-Inverse Based Normaliation
}

cor.test(join$NT.proBNP, join$bld_NT_ProBNP_w4) # 0.25

join$age <- join$agedays_w4.x / 365.25
join$sex <- join$sex.x

null <- summary(lm(bld_NT_ProBNP_w4 ~ age + sex, data=join))$r.squared
full <- summary(lm(bld_NT_ProBNP_w4 ~ age + sex + join$NT.proBNP, data=join))$r.squared
print(round(100*(full - null), 3))

# 8.13

# Add in PRS for BNP in LBC
prs <- read.table("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_PRS/lbc1936_genomewide_bnp.all.score", header = T)

names(prs)[2] <- 'lbc36no'
prs <- prs[-1]
names(prs)[2] <- 'BNPPRS'

join <- left_join(join, prs, by = 'lbc36no')

null <- summary(lm(bld_NT_ProBNP_w4 ~ age + sex, data=join))$r.squared
null_prs <- summary(lm(bld_NT_ProBNP_w4 ~ age + sex + BNPPRS, data=join))$r.squared
full <- summary(lm(bld_NT_ProBNP_w4 ~ age + sex + join$NT.proBNP, data=join))$r.squared
all <- summary(lm(bld_NT_ProBNP_w4 ~ age + sex + join$NT.proBNP + BNPPRS, data=join))$r.squared
print(round(100*(full - null), 3))

# > null
# [1] 0.007537335
# > null_prs
# [1] 0.04060833
# > full
# [1] 0.08884607
# > all
# [1] 0.09824992

table(is.na(join$NT.proBNP))
table(is.na(join$bld_NT_ProBNP_w4)) #  500 (2 missing values)


library(ggpubr)
plot2 <- ggplot(join, aes(x= bld_NT_ProBNP_w4, y= NT.proBNP)) +
geom_point(colour = "firebrick1", size = 1) +
geom_smooth(method='lm', colour = "firebrick1") + 
theme(axis.text.x=element_text(size=20),     
      axis.text.y=element_text(size=20),
      axis.title.x=element_text(size=20),
      axis.title.y=element_text(size=20)) + 
stat_cor(aes(label = ..r.label..), method = "pearson", cor.coef.name = "r", size = 7.5) +
xlab("Nt-proBNP protein") + ylab("Nt-proBNP EpiScore") + 
theme_classic() # , label.x = 4.5, label.y = 4.4


### PLOT TOGETHER
library(patchwork)
pdf('/Local_Scratch/Danni/GDFBNP/00_Updated_results/Results_scores/CORR_PLOT_JOINT_LBC.pdf', width = 12, height = 5)
plot1 + plot2 
dev.off()


###################################################################################

### Use measured proteins as test sets in LBC first 

###################################################################################

cd /Local_Scratch/Danni/GDFBNP/04_LBC/20k_GS_training/

screen

R

library(tidyverse)

W1 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_Analysis_20k/04_LBC/20k_GS_training/projections_LBC_W1_241122.csv")
W2 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_Analysis_20k/04_LBC/20k_GS_training/projections_LBC_W2_241122.csv")
W3 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_Analysis_20k/04_LBC/20k_GS_training/projections_LBC_W3_241122.csv")
W4 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_Analysis_20k/04_LBC/20k_GS_training/projections_LBC_W4_241122.csv")


# > dim(W1)
# [1] 906   3
# > dim(W2)
# [1] 801   3
# > dim(W3)
# [1] 619   3
# > dim(W4)
# [1] 507   3


W2_set1 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_Analysis_20k/04_LBC/20k_GS_training/projections_LBC_W2_set1_241122.csv")
W2_set2 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_Analysis_20k/04_LBC/20k_GS_training/projections_LBC_W2_set2_241122.csv")
W3_set1 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_Analysis_20k/04_LBC/20k_GS_training/projections_LBC_W3_set1_241122.csv")
W3_set2 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_Analysis_20k/04_LBC/20k_GS_training/projections_LBC_W3_set2_241122.csv")


BNP_w3 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_Analysis_20k/04_LBC/LBC_testing/BNP_w3.csv")
BNP_w4 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_Analysis_20k/04_LBC/LBC_testing/BNP_w4.csv")
GDF_w2 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_Analysis_20k/04_LBC/LBC_testing/GDF_w2.csv")
GDF_w4 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_Analysis_20k/04_LBC/LBC_testing/GDF_w4.csv")

# > dim(BNP_w3)
# [1] 616  73
# > dim(BNP_w4)
# [1] 502  73
# > dim(GDF_w2)
# [1] 762  73
# > dim(GDF_w4)
# [1] 322  73


### GDF W2
names(W2)[3] <- 'Basename'
join <- left_join(GDF_w2, W2, by = 'Basename')

library(bestNormalize)
for(i in colnames(join)[c(17,74)]){ 
  join[,i]<- orderNorm(join[,i])$x.t # Rank-Inverse Based Normaliation
}

cor.test(join$GDF15.20k.with.450k.array, join$GDF15_W2) # 0.21

join$age <- join$agedays_w2.x / 365.25
join$sex <- join$sex.x

null <- summary(lm(GDF15_W2 ~ age + sex, data=join))$r.squared
full <- summary(lm(GDF15_W2 ~ age + sex + join$GDF15.20k.with.450k.array, data=join))$r.squared
print(round(100*(full - null), 3))

# GDF15 pQTL 450k
# LBC W2: 2.68% 

table(is.na(join$GDF15.20k.with.450k.array))
table(is.na(join$GDF15_W2)) # 762 with both available 

## try GDF score grimage 

LBC_grim <- read.csv("/Cluster_Filespace/Marioni_Group/LBC/LBC_methylation/LBC_clock_output_3489.csv")
# LBC_grim <- LBC_grim[c(2,78:82,84:85)]
names(LBC_grim)[1] <- "Basename"

join <- left_join(join, LBC_grim, by = 'Basename')

null <- summary(lm(GDF15_W2 ~ age + sex, data=join))$r.squared
full <- summary(lm(GDF15_W2 ~ age + sex + join$DNAmGDF15, data=join))$r.squared
print(round(100*(full - null), 3))

# GDF15 pQTL 450k
# LBC W2: 2.68% 


library(ggpubr)
plot1 <- ggplot(join, aes(x=GDF15_W2, y=GDF15.20k.with.450k.array)) +
geom_point(colour = "tan1", size = 1) +
geom_smooth(method='lm', colour = "tan1") + 
theme(axis.text.x=element_text(size=20),     
      axis.text.y=element_text(size=20),
      axis.title.x=element_text(size=20),
      axis.title.y=element_text(size=20)) + 
stat_cor(aes(label = ..r.label..), method = "pearson", cor.coef.name = "r", size = 7.5) +
xlab("GDF15 protein W2") + ylab("GDF15 EpiScore") # , label.x = 4.5, label.y = 4.4

# Compare protein and episcore associations with stroke mask in w2

prot_model <- lm(Ig ~ age + sex + GDF15_W2, data = join)
score_model <- lm(Ig ~ age + sex + GDF15.20k.with.450k.array, data = join)



### GDF W4
names(W4)[3] <- 'Basename'
join <- left_join(GDF_w4, W4, by = 'Basename')

library(bestNormalize)
for(i in colnames(join)[c(18,74)]){ 
  join[,i]<- orderNorm(join[,i])$x.t # Rank-Inverse Based Normaliation
}

cor.test(join$GDF15.20k.with.450k.array, join$GDF15_W4) # 0.31

join$age <- join$agedays_w4.x / 365.25
join$sex <- join$sex.x

null <- summary(lm(GDF15_W4 ~ age + sex, data=join))$r.squared
full <- summary(lm(GDF15_W4 ~ age + sex + join$GDF15.20k.with.450k.array, data=join))$r.squared
print(round(100*(full - null), 3))

# GDF15 pQTL 450k
# LBC W4: 9.14% 

table(is.na(join$GDF15.20k.with.450k.array))
table(is.na(join$GDF15_W4)) # 322 with both available 

library(ggpubr)
plot2 <- ggplot(join, aes(x=GDF15_W4, y=GDF15.20k.with.450k.array)) +
geom_point(colour = "tan1", size = 1) +
geom_smooth(method='lm', colour = "tan1") + 
theme(axis.text.x=element_text(size=20),     
      axis.text.y=element_text(size=20),
      axis.title.x=element_text(size=20),
      axis.title.y=element_text(size=20)) + 
stat_cor(aes(label = ..r.label..), method = "pearson", cor.coef.name = "r", size = 7.5) +
xlab("GDF15 protein W4") + ylab("GDF15 EpiScore") # , label.x = 4.5, label.y = 4.4

prot_model <- lm(Ig ~ age + sex + GDF15_W4, data = join)
score_model <- lm(Ig ~ age + sex + GDF15.20k.with.450k.array, data = join)

join <- join[which(join$enl_peri_space_hl_w2 %in% c('1','0')),]
prot_model <- glm(enl_peri_space_hl_w2 ~ age + sex + GDF15_W4 + ICV_mm3_wX, 
  family=binomial(link='logit'), data = join)
score_model <- glm(enl_peri_space_hl_w2 ~ age + sex + GDF15.20k.with.450k.array + ICV_mm3_wX,
family=binomial(link='logit'), data = join)



### BNP w4
names(W4)[3] <- 'Basename'
join <- left_join(BNP_w4, W4, by = 'Basename')

join$bld_NT_ProBNP_w4 <- as.numeric(join$bld_NT_ProBNP_w4)

library(bestNormalize)
for(i in colnames(join)[c(20,75)]){ 
  join[,i]<- orderNorm(join[,i])$x.t # Rank-Inverse Based Normaliation
}

cor.test(join$Nt.proBNP.20k.with.450k.array, join$bld_NT_ProBNP_w4) # 0.18

join$age <- join$agedays_w4.x / 365.25
join$sex <- join$sex.x

null <- summary(lm(bld_NT_ProBNP_w4 ~ age + sex, data=join))$r.squared
full <- summary(lm(bld_NT_ProBNP_w4 ~ age + sex + join$Nt.proBNP.20k.with.450k.array, data=join))$r.squared
print(round(100*(full - null), 3))

# BNP pQTL 450k
# LBC W4: 3.795%

table(is.na(join$Nt.proBNP.20k.with.450k.array))
table(is.na(join$bld_NT_ProBNP_w4)) #  500 (2 missing values)


library(ggpubr)
plot3 <- ggplot(join, aes(x=bld_NT_ProBNP_w4, y=Nt.proBNP.20k.with.450k.array)) +
geom_point(colour = "firebrick1", size = 1) +
geom_smooth(method='lm', colour = "firebrick1") + 
theme(axis.text.x=element_text(size=20),     
      axis.text.y=element_text(size=20),
      axis.title.x=element_text(size=20),
      axis.title.y=element_text(size=20)) + 
stat_cor(aes(label = ..r.label..), method = "pearson", cor.coef.name = "r", size = 7.5) +
xlab("Nt-proBNP protein W4") + ylab("Nt-proBNP EpiScore") # , label.x = 4.5, label.y = 4.4




# ### BNP w3
# names(W3)[3] <- 'Basename'
# join <- left_join(BNP_w3, W3, by = 'Basename')

# join$bld_NT_ProBNP_w3 <- as.numeric(join$bld_NT_ProBNP_w3)

# library(bestNormalize)
# for(i in colnames(join)[c(19,75)]){ 
#   join[,i]<- orderNorm(join[,i])$x.t # Rank-Inverse Based Normaliation
# }

# cor.test(join$Nt.proBNP.20k.with.450k.array, join$bld_NT_ProBNP_w3) # 0.13

# join$age <- join$agedays_w3.x / 365.25
# join$sex <- join$sex.x

# null <- summary(lm(bld_NT_ProBNP_w3 ~ age + sex, data=join))$r.squared
# full <- summary(lm(bld_NT_ProBNP_w3 ~ age + sex + join$Nt.proBNP.20k.with.450k.array, data=join))$r.squared
# print(round(100*(full - null), 3))

# # BNP pQTL 450k
# # LBC W4: 1.118%

# table(is.na(join$Nt.proBNP.20k.with.450k.array))
# table(is.na(join$bld_NT_ProBNP_w3)) #  616 with both values 


# library(ggpubr)
# plot4 <- ggplot(join, aes(x=bld_NT_ProBNP_w3, y=Nt.proBNP.20k.with.450k.array)) +
# geom_point(colour = "firebrick1", size = 1) +
# geom_smooth(method='lm', colour = "firebrick1") + 
# theme(axis.text.x=element_text(size=20),     
#       axis.text.y=element_text(size=20),
#       axis.title.x=element_text(size=20),
#       axis.title.y=element_text(size=20)) + 
# stat_cor(aes(label = ..r.label..), method = "pearson", cor.coef.name = "r", size = 7.5) +
# xlab("Nt-proBNP protein W3") + ylab("Nt-proBNP EpiScore") # , label.x = 4.5, label.y = 4.4


# ### PLOT TOGETHER
# library(patchwork)
# pdf('/Local_Scratch/Danni/GDFBNP/04_LBC/LBC_testing/PLOT_JOINT.pdf', width = 15, height = 12)
# plot1 + plot2 + plot3 + plot4
# dev.off()


# ############################################################################################

# ### PROJECT BY SET IN W2 W3 - test performance by set 

# ### GDF15 W2 - set 1

# W2 <- W2_set1
# names(W2)[3] <- 'Basename'
# join <- left_join(GDF_w2, W2, by = 'Basename')

# # Read in target file with matched IDs for methylation data 
# target <- read.csv("/Cluster_Filespace/Marioni_Group/LBC/LBC_methylation/archive/target_QC_age_sex_date.csv")
# dat3 <- target %>% filter(WAVE == "2") # 1342 people
# dat3  <- dat3 %>% filter(cohort == "LBC36") # 906 people 
# dat3 <- dat3 %>% select('Basename', 'set')
# join <- left_join(join, dat3, by = 'Basename')
# join <- filter(join, set %in% c(1))

# library(bestNormalize)
# for(i in colnames(join)[c(17,74)]){ 
#   join[,i]<- orderNorm(join[,i])$x.t # Rank-Inverse Based Normaliation
# }


# cor.test(join$GDF15.20k.with.450k.array, join$GDF15_W2) # 0.21 - to 0.28

# join$age <- join$agedays_w2.x / 365.25
# join$sex <- join$sex.x

# join <- join[complete.cases(join$GDF15.20k.with.450k.array),]
# join <- join[complete.cases(join$GDF15_W2),]

# null <- summary(lm(GDF15_W2 ~ age + sex, data=join))$r.squared
# full <- summary(lm(GDF15_W2 ~ age + sex + join$GDF15.20k.with.450k.array, data=join))$r.squared
# print(round(100*(full - null), 3))




# ### GDF15 W2 - set 2
# W2 <- W2_set2
# names(W2)[3] <- 'Basename'
# join <- left_join(GDF_w2, W2, by = 'Basename')

# # Read in target file with matched IDs for methylation data 
# target <- read.csv("/Cluster_Filespace/Marioni_Group/LBC/LBC_methylation/target_QC_age_sex_date.csv")

# # Filter to wave 1
# dat3 <- target %>% filter(WAVE == "2") # 1342 people

# # Filter to LBC36
# dat3  <- dat3 %>% filter(cohort == "LBC36") # 906 people 
# dat3 <- dat3 %>% select('Basename', 'set')

# join <- left_join(join, dat3, by = 'Basename')
# join <- filter(join, set %in% c(2))

# join <- join[complete.cases(join$GDF15.20k.with.450k.array),]
# join <- join[complete.cases(join$GDF15_W2),]


# library(bestNormalize)
# for(i in colnames(join)[c(17,74)]){ 
#   join[,i]<- orderNorm(join[,i])$x.t # Rank-Inverse Based Normaliation
# }

# cor.test(join$GDF15.20k.with.450k.array, join$GDF15_W2) # 0.21 - to 0.28

# join$age <- join$agedays_w2.x / 365.25
# join$sex <- join$sex.x


# null <- summary(lm(GDF15_W2 ~ age + sex, data=join))$r.squared
# full <- summary(lm(GDF15_W2 ~ age + sex + join$GDF15.20k.with.450k.array, data=join))$r.squared
# print(round(100*(full - null), 3))

# # GDF15 pQTL 450k
# # LBC W2: 2.68% - to 5.674%




# ### BNP w3 - set 1
# W3 <- W3_set1
# names(W3)[3] <- 'Basename'
# join <- left_join(BNP_w3, W3, by = 'Basename')

# join$bld_NT_ProBNP_w3 <- as.numeric(join$bld_NT_ProBNP_w3)

# # Read in target file with matched IDs for methylation data 
# target <- read.csv("/Cluster_Filespace/Marioni_Group/LBC/LBC_methylation/target_QC_age_sex_date.csv")

# # Filter to wave 1
# dat3 <- target %>% filter(WAVE == "3") # 1342 people

# # Filter to LBC36
# dat3  <- dat3 %>% filter(cohort == "LBC36") # 906 people 
# dat3 <- dat3 %>% select('Basename', 'set')

# join <- left_join(join, dat3, by = 'Basename')
# join <- filter(join, set %in% c(1))

# join <- join[complete.cases(join$Nt.proBNP.20k.with.450k.array),]
# join <- join[complete.cases(join$bld_NT_ProBNP_w3),]

# library(bestNormalize)
# for(i in colnames(join)[c(19,75)]){ 
#   join[,i]<- orderNorm(join[,i])$x.t # Rank-Inverse Based Normaliation
# }

# cor.test(join$Nt.proBNP.20k.with.450k.array, join$bld_NT_ProBNP_w3) # 0.13

# join$age <- join$agedays_w3.x / 365.25
# join$sex <- join$sex.x

# null <- summary(lm(bld_NT_ProBNP_w3 ~ age + sex, data=join))$r.squared
# full <- summary(lm(bld_NT_ProBNP_w3 ~ age + sex + join$Nt.proBNP.20k.with.450k.array, data=join))$r.squared
# print(round(100*(full - null), 3))



# ### BNP w3 - set 2
# W3 <- W3_set2
# names(W3)[3] <- 'Basename'
# join <- left_join(BNP_w3, W3, by = 'Basename')

# join$bld_NT_ProBNP_w3 <- as.numeric(join$bld_NT_ProBNP_w3)

# # Read in target file with matched IDs for methylation data 
# target <- read.csv("/Cluster_Filespace/Marioni_Group/LBC/LBC_methylation/target_QC_age_sex_date.csv")

# # Filter to wave 1
# dat3 <- target %>% filter(WAVE == "3") # 1342 people

# # Filter to LBC36
# dat3  <- dat3 %>% filter(cohort == "LBC36") # 906 people 
# dat3 <- dat3 %>% select('Basename', 'set')

# join <- left_join(join, dat3, by = 'Basename')
# join <- filter(join, set %in% c(2))

# join <- join[complete.cases(join$Nt.proBNP.20k.with.450k.array),]
# join <- join[complete.cases(join$bld_NT_ProBNP_w3),]

# library(bestNormalize)
# for(i in colnames(join)[c(19,75)]){ 
#   join[,i]<- orderNorm(join[,i])$x.t # Rank-Inverse Based Normaliation
# }

# cor.test(join$Nt.proBNP.20k.with.450k.array, join$bld_NT_ProBNP_w3) 

# join$age <- join$agedays_w3.x / 365.25
# join$sex <- join$sex.x

# null <- summary(lm(bld_NT_ProBNP_w3 ~ age + sex, data=join))$r.squared
# full <- summary(lm(bld_NT_ProBNP_w3 ~ age + sex + join$Nt.proBNP.20k.with.450k.array, data=join))$r.squared
# print(round(100*(full - null), 3))


# ##################################################################################################

# ### JOIN W2 and W3 by sets and adjust for set when doing regressions


# ### GDF15 W2 

# W2 <- rbind(W2_set1, W2_set2)
# names(W2)[3] <- 'Basename'
# join <- left_join(GDF_w2, W2, by = 'Basename')

# # Read in target file with matched IDs for methylation data 
# target <- read.csv("/Cluster_Filespace/Marioni_Group/LBC/LBC_methylation/target_QC_age_sex_date.csv")

# # Filter to wave 1
# dat3 <- target %>% filter(WAVE == "2") # 1342 people

# # Filter to LBC36
# dat3  <- dat3 %>% filter(cohort == "LBC36") # 906 people 
# dat3 <- dat3 %>% select('Basename', 'set')

# join <- left_join(join, dat3, by = 'Basename')
# join <- filter(join, set %in% c(1,2))

# join$age <- join$agedays_w2.x / 365.25
# join$sex <- join$sex.x

# join <- join[complete.cases(join$GDF15.20k.with.450k.array),]
# join <- join[complete.cases(join$GDF15_W2),]

# library(bestNormalize)
# for(i in colnames(join)[c(17,74)]){ 
#   join[,i]<- orderNorm(join[,i])$x.t # Rank-Inverse Based Normaliation
# }

# cor.test(join$GDF15.20k.with.450k.array, join$GDF15_W2) # 0.21 - to 0.28


# null <- summary(lm(GDF15_W2 ~ age + sex + factor(set), data=join))$r.squared
# full <- summary(lm(GDF15_W2 ~ age + sex + factor(set) + join$GDF15.20k.with.450k.array, data=join))$r.squared
# print(round(100*(full - null), 3))

# # GDF15 pQTL 450k
# # LBC W2: 2.68% - rises to 6.844% by set 1, but 6.124 when sets combined 




# ### BNP w3 - set 2
# W3 <- rbind(W3_set1, W3_set2)
# names(W3)[3] <- 'Basename'
# join <- left_join(BNP_w3, W3, by = 'Basename')

# join$bld_NT_ProBNP_w3 <- as.numeric(join$bld_NT_ProBNP_w3)

# # Read in target file with matched IDs for methylation data 
# target <- read.csv("/Cluster_Filespace/Marioni_Group/LBC/LBC_methylation/target_QC_age_sex_date.csv")

# # Filter to wave 1
# dat3 <- target %>% filter(WAVE == "3") # 1342 people

# # Filter to LBC36
# dat3  <- dat3 %>% filter(cohort == "LBC36") # 906 people 
# dat3 <- dat3 %>% select('Basename', 'set')

# join <- left_join(join, dat3, by = 'Basename')
# join <- filter(join, set %in% c(1,2))


# join <- join[complete.cases(join$Nt.proBNP.20k.with.450k.array),]
# join <- join[complete.cases(join$bld_NT_ProBNP_w3),]

# library(bestNormalize)
# for(i in colnames(join)[c(19,75)]){ 
#   join[,i]<- orderNorm(join[,i])$x.t # Rank-Inverse Based Normaliation
# }

# cor.test(join$Nt.proBNP.20k.with.450k.array, join$bld_NT_ProBNP_w3) # 0.13 - to 0.09 - back to 0.13 when joined together

# join$age <- join$agedays_w3.x / 365.25
# join$sex <- join$sex.x

# null <- summary(lm(bld_NT_ProBNP_w3 ~ age + factor(sex) + factor(set), data=join))$r.squared
# full <- summary(lm(bld_NT_ProBNP_w3 ~ age + factor(sex) + factor(set) + join$Nt.proBNP.20k.with.450k.array, data=join))$r.squared
# print(round(100*(full - null), 3))

# # BNP pQTL 450k
# # LBC W4: 1.118% - to 2.937 - to 1.488 when combined across sets 


