############################################################################################################

### COLLATE COVARIATES FOR HANNAH TO USE

############################################################################################################

### load in health/lifestyle phenotypes - taken from daniels cog EWAS
# original location sourced from: pheno <- read.spss("/Cluster_Filespace/Marioni_Group/Anna/DNAm phenoAge/LBC1936_PheWAS_Ageing_Multi_omics_factor_analysis_RM_06APR2018.sav", to.data.frame=T)
library("foreign")
pheno <- read.spss("/Local_Scratch/Danni/GDFBNP/04_LBC/original_data/LBC1936_PheWAS_Ageing_Multi_omics_factor_analysis_RM_06APR2018.sav", to.data.frame=T)
ph <- pheno[,c("lbc36no", "yrsedu_w1","bmi_w1","smokcat_w1", "alcunitwk_w1", "hibp_w1", "hadsd_w1", "height_w1", "fev_w1", "griprh_w1", "griplh_w1", "depind_w1", "sixmwk_w1")]
ph$lbc36no <- gsub(" ", "", ph$lbc36no)

# remove impossible SIMD value
ph$depind_w1[ph$depind_w1==99999] <- NA

# Add age and sex and convert age to years format 
library(haven)
library(tidyverse)
LBC_phenotypes <- read_sav("/Local_Scratch/Danni/GDFBNP/04_LBC/original_data/LBC1936_Blood_DNAm_And_Brain_DNAm_RM_22APR2020.sav")
LBC_phenotypes <- as.data.frame(LBC_phenotypes)
LBC_phenotypes <- LBC_phenotypes %>% select('lbc36no', 'sex', 'agedays_w1')
LBC_phenotypes$ageyearsw1 <- LBC_phenotypes$agedays_w1 / 365.25
ph <- left_join(LBC_phenotypes, ph, by = 'lbc36no')

# Get W1 LBC1936 DNAm information from the target file - with WBCs and technical variables 
library(tidyverse)
target <- read.csv("/Cluster_Filespace/Marioni_Group/LBC/LBC_methylation/target_QC_age_sex_date.csv")
dat3 <- target %>% filter(WAVE == "1") 
dat3  <- dat3 %>% filter(cohort == "LBC36")
names(dat3)[4] <- 'lbc36no' # 2 missing BMI and 157 missng alc freq
dat3 <- dat3[-c(6,14)]
dat3 <- left_join(ph, dat3, by = 'lbc36no')

# Add epismoker
epismoker <- readRDS("/Cluster_Filespace/Marioni_Group/Danni/lbc_epismoker.rds")
names(epismoker)[1] <- 'Basename'
dat3 <- left_join(dat3, epismoker, by = 'Basename')

# Save out for Hannah 
write.csv(dat3, '/Local_Scratch/Danni/GDFBNP/04_LBC/original_data/LBC_processed_data.csv', row.names = F)

## Take a look at proposed covariate levels
table(dat3$smokcat_w1)
  # never smoked      ex-smoker current smoker
  #          501            465            125
table(is.na(dat3$alcunitwk_w1)) # complete (1091)
table(is.na(dat3$bmi_w1)) # 2 missing individuals
table(dat3$yrsedu_w1)
  # 7   8   9  10  11  12  13  14
  # 1   3  45 584 179 157 119   3

# Proposed basic: age in years + factor(sex)
# Proposed WBC: age in years + factor(sex) + neut + lymph + mono + eosin + baso
# Proposed full: age in years + factor(sex) + neut + lymph + mono + eosin + baso + factor(smokcat_w1) + factor(yrsedu_w1) + bmi_w1 + depind_w1 + alcunitwk_w1 


############################################################################################################

### ADJUST EPISCORES FOR TECHNICAL DNAm COVS

############################################################################################################

library(tidyverse)

# read in W1 episcore projections for LBC1936 
W1 <- read.csv("/Local_Scratch/Danni/GDFBNP/00_Updated_results/Results_scores/projections_LBC_W1_241122.csv")
names(W1)[3] <- 'Basename'

# Read in target file with technical variables for adjustment 
library(tidyverse)
target <- read.csv("/Cluster_Filespace/Marioni_Group/LBC/LBC_methylation/target_QC_age_sex_date.csv")
dat3 <- target %>% filter(WAVE == "1") 
dat3  <- dat3 %>% filter(cohort == "LBC36")

# Join together for regression
W1 <- left_join(W1, dat3, by = 'Basename') # 906 individuals with info 

# # Adjust episcores for technical variables and take residuals forward for SEM analyses
# phen <- W1
# for(i in colnames(phen)[1:2]){ # set columns to the protein EpiScores you want to adjust 
#   phen[,i]<- lm(phen[,i] ~ factor(set) + factor(date) + factor(array), 
#                                     na.action = na.exclude, data = phen)$residuals
# }

# cor.test(phen[,1], W1[,1])
# cor.test(phen[,2], W1[,2])

library(lme4)
# Adjust episcores by technical variables as random effects
phen <- W1
phen$set <- as.factor(phen$set)
phen$array <- as.factor(phen$array)
phen$date <- as.factor(phen$date)


mod <- lmer(phen[,1] ~ (1 | set) + (1 | date) + (1 | array), 
                                    na.action = na.exclude, data = phen)
phen[,1] <- resid(mod)

mod2 <- lmer(phen[,2] ~ (1 | set) + (1 | date) + (1 | array), 
                                    na.action = na.exclude, data = phen)
phen[,2] <- resid(mod2)


# list <- c(3,5)
# for(i in list) {
#     EpiScore[,i] <- scale(resid(lmer(EpiScore[,i] ~ (1 | set) + (1 | date) + (1 | array), na.action=na.exclude)))
# }

# Rank-Inverse Based Normaliation 
library(bestNormalize)
for(i in colnames(phen)[1:2]){ 
  phen[,i]<- orderNorm(phen[,i])$x.t
}

# Save out for Hannah 
phen <- phen[c(1:3)]
names(phen) <- c('GDF15.20k.with.450k.array', 'Nt.proBNP.20k.with.450k.array', 'Basename')
write.csv(phen, '/Local_Scratch/Danni/GDFBNP/04_LBC/original_data/LBC_episcores_adjusted_technical.csv', row.names = F)

############################################################################################################

### CORRECT W2 scores for use in cox PH models 

# read in W2 episcore projections for LBC1936 

set1 <- read.csv("/Local_Scratch/Danni/GDFBNP/00_Updated_results/Results_scores/projections_LBC_W2_set1_241122.csv")
set2 <- read.csv("/Local_Scratch/Danni/GDFBNP/00_Updated_results/Results_scores/projections_LBC_W2_set2_241122.csv")
scores <- rbind(set1, set2)

# scores <- read.csv("/Local_Scratch/Danni/GDFBNP/00_Updated_results/Results_scores/projections_LBC_W2_241122.csv")
names(scores)[3] <- 'Basename'

# Read in target file with technical variables for adjustment 
library(tidyverse)
target <- read.csv("/Cluster_Filespace/Marioni_Group/LBC/LBC_methylation/target_QC_age_sex_date.csv")
dat3 <- target %>% filter(WAVE == "2") 
dat3  <- dat3 %>% filter(cohort == "LBC36")

# Join together for regression
scores <- left_join(scores, dat3, by = 'Basename') # 796 individuals with info 

# # Adjust episcores for technical variables and take residuals forward for SEM analyses
# phen <- W1
# for(i in colnames(phen)[1:2]){ # set columns to the protein EpiScores you want to adjust 
#   phen[,i]<- lm(phen[,i] ~ factor(set) + factor(date) + factor(array), 
#                                     na.action = na.exclude, data = phen)$residuals
# }

# cor.test(phen[,1], W1[,1])
# cor.test(phen[,2], W1[,2])

# Adjust episcores by technical variables as random effects
phen <- scores
phen$set <- as.factor(phen$set)
phen$array <- as.factor(phen$array)
phen$date <- as.factor(phen$date)

library(lme4)

mod <- lmer(phen[,1] ~ (1 | set) + (1 | date) + (1 | array), 
                                    na.action = na.exclude, data = phen)
phen[,1] <- resid(mod)

mod2 <- lmer(phen[,2] ~ (1 | set) + (1 | date) + (1 | array), 
                                    na.action = na.exclude, data = phen)
phen[,2] <- resid(mod2)


# list <- c(3,5)
# for(i in list) {
#     EpiScore[,i] <- scale(resid(lmer(EpiScore[,i] ~ (1 | set) + (1 | date) + (1 | array), na.action=na.exclude)))
# }

# Rank-Inverse Based Normaliation 
library(bestNormalize)
for(i in colnames(phen)[1:2]){ 
  phen[,i]<- orderNorm(phen[,i])$x.t
}

# Save out for use in cox PH
# phen <- phen[c(1:3)]
names(phen) <- c('GDF15_score', 'Nt.proBNP_score', 'Basename')
write.csv(phen, '/Local_Scratch/Danni/GDFBNP/04_LBC/original_data/W2_LBC_episcores_adjusted_technical.csv', row.names = F)
