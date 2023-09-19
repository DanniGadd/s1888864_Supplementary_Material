####################################################################################

### COX MODELLING FOR GS BIOMARKERS 

####################################################################################

screen

R

# Load in packages 
library(survival)
library(kinship2)
library(coxme)
library(readxl)
library(tidyverse)
library(gdata)

# Read in the prepped pedigree file to cluster and create kinship matrix
ped <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Proxies_stradl/KORA_train_recieved_020221/Cox/pedigree_formatted.csv")
kin <- with(ped, pedigree(volid, father, mother, sex, famid=famid))
kin_model <- kinship(kin) 

# Load function to Extract Lmekin Results	
extract_coxme_table <- function (mod){
  beta <- mod$coefficients #$fixed is not needed
  nvar <- length(beta)
  nfrail <- nrow(mod$var) - nvar
  se <- sqrt(diag(mod$var)[nfrail + 1:nvar])
  z<- round(beta/se, 2)
  p<- 1 - pchisq((beta/se)^2, 1)
  table=data.frame(cbind(beta,se,z,p))
  return(table)
}


####################################################################################

#### PREP PHENOTYPE FILE WITH AGE ALIVE AND AGE DEATH INFO 

####################################################################################

## DEMOGRAPHICS

## Read in GS 20k age/sex base file and add prevalent data to it 
all <- read.csv("/Cluster_Filespace/Marioni_Group/GS/GS_dataset/agemonths.csv")
names(all)[2] <- "age"
names(all)[1] <- "Sample_Name"
dim(all) # 24088   3
prevalent <- read.csv("/Cluster_Filespace/Marioni_Group/GS/GS_dataset/PCQ/disease.csv")
names(prevalent)[1] <- "Sample_Name"
merged_prev <- merge(all, prevalent, by = "Sample_Name", all = T)
dim(merged_prev) # 24092   115

## SURVIVAL INFO

# Read in deaths data and count how many individuals have died since oct 2020 (when GP data was last sampled)
age_dead <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/GS_COVID/Diabetes_update_codes_yipeng/2022-03-10_age_at_death.csv")
length(which(age_dead$dod_ym > 202010)) # 227 individuals since oct 2020 

## extract year/month of death as separate variables ##
age_dead$y_dead <- as.numeric(substr(age_dead$dod_ym, 1, 4))
age_dead$m_dead <- as.numeric(substr(age_dead$dod_ym, 5, 6))

# Assign the dead individuals as a separate subset
age_dead_include <- age_dead[which(age_dead$dod_ym <= 202010),] # 1350
age_dead_exclude <- age_dead[which(age_dead$dod_ym > 202010),] # 227

# Calculate a more exact estimate (by year and month) for age of death in the included 1350 individuals that died 
age_dead_include$y_diff <- age_dead_include$y_dead - age_dead_include$yob
age_dead_include$m_diff <- (age_dead_include$m_dead - age_dead_include$mob)/12
age_dead_include$diff <- age_dead_include$y_diff + age_dead_include$m_diff

# Work out those who are alive (i.e. not in the list of 1350 dead people from above)
age_alive <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/GS_COVID/Diabetes_update_codes_yipeng/2022-03-10_age_alive.csv")
age_alive <- age_alive[c(1:3,6,8)]
age_alive = age_alive[!age_alive$id %in% age_dead_include$id,] # 22738 individuals who are not in the dead people we include - this covers the remainder of GS who did not die

# Ensure that all individuals in the 22738 sample are coded as alive and all individuals in the dead file are coded as such 
age_alive$dead <- 0
age_dead_include$dead <- 1

# Find age the 'alive' people were in oct 2020
age_alive$y_diff <- 2020 - age_alive$yob
age_alive$m_diff <- (10 - age_alive$mob)/12
age_alive$diff <- age_alive$y_diff + age_alive$m_diff

# The included dead people will have their age at death taken forward (i.e. pre 2020) as the 'aged' column in cox loops - this is now calculated
# The excluded dead people will be classed as alive and will have their age at 2020 taken forward as the 'aged' column - we have just calculated this as part of the wider group of alive individuals

# Subset to just the cols needed for joining dead and alive data
age_alive <- age_alive[c(1:4,8)] 
age_dead_include <- age_dead_include[c(1,2,3,7,13)]

# > dim(age_alive)
# [1] 22738     5
# > dim(age_dead_include)
# [1] 1350    5

# Bind the rows of the alive and dead files together for the whole GS sample 
names(age_dead_include) <- c("id", "yob", "mob", "dead", "aged")
names(age_alive) <- c("id", "yob", "mob", "dead", "aged")
age = rbind(age_alive, age_dead_include)
dim(age) # 24088     5
table(age$dead)

#     0     1
# 22738  1350

names(age)[1] <- "Sample_Name"

## Add survival info to the 20k base file 
d1 <- left_join(merged_prev, age, by = "Sample_Name")
dim(d1) # [1] 24092   119

# Create a subset with DOBMOB to use to filter cases by in the diabetes code processing below
d2 <- d1[c(1,2,116,117)]

# # Calculate prevalent AD (it was split into males and females)
# d1$AD <- 0 
# d1[d1$Sample_Name %in% d1[(d1$alzheimers_M %in% 1 | d1$alzheimers_F %in% 1),"Sample_Name"],"AD"] <- 1

####################################################################################

## Read in phenotypes to covary for in the fully-adjusted models 

# details of how Rob preprocessed original traits used in the epiclocks paper:
# "Outliers were defined as those values which were beyond 3.5 standard deviations from the mean for a given trait. 
# These outliers were removed prior to analyses. Body mass index was log-transformed. 
# To reduce skewness in the distribution of alcohol consumption and smoking pack years, a log(units +1) or log(pack years +1) transformation was performed."

# Load in the original d1 files

# Alcohol (units) and (usual)
alcohol <- read.csv("/Cluster_Filespace/Marioni_Group/GS/GS_dataset/PCQ/alcohol.csv")
names(alcohol)[1] <- "Sample_Name"

# BMI at baseline 
BMI <- read.csv("/Cluster_Filespace/Marioni_Group/GS/GS_dataset/clinical/body.csv")
names(BMI)[1] <- "Sample_Name"

# SIMD
simd <- read.csv("/Cluster_Filespace/Marioni_Group/GS/GS_dataset/clinical/SIMD.csv")
names(simd)[1] <- "Sample_Name"

# EA
ea <- read.csv("/Cluster_Filespace/Marioni_Group/GS/GS_dataset/PCQ/education.csv")
names(ea)[1] <- "Sample_Name"

# Join up the covariate data to the new d1 file in 20k 
d1 <- left_join(d1, alcohol, by = "Sample_Name")
d1 <- left_join(d1, BMI, by = "Sample_Name")
d1 <- left_join(d1, simd, by = "Sample_Name")
d1 <- left_join(d1, ea, by = "Sample_Name")

# Read in epismoekr 
w1 <- readRDS("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/DNAm_preps/wave1_epismoker.rds")
w3 <- readRDS("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/DNAm_preps/wave3_epismoker.rds")
w4 <- readRDS("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/DNAm_preps/wave4_epismoker.rds")
bind <- rbind(w1, w3)
bind <- rbind(bind, w4) # 18779 individuals with DNAm epismoker calculated 
bind$Sample_Sentrix_ID <- row.names(bind)
bind <- bind[-1]

# Join Sample_Name info 
target <- readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/GS20k/GS20k_Targets.rds")
bind <- left_join(bind, target, by = "Sample_Sentrix_ID")
bind <- bind[c(1,3)]

# Join to the data 
d1 <- left_join(d1, bind, by = "Sample_Name")

## Remove outlying values from key covairates used in models (>3.5 SDs from mean in each direction)
list <- c("units", "usual", "bmi", "rank", "years")
for(i in list){ 
  
  cutoff1 = mean(d1[,i], na.rm = T) + 3.5*sd(d1[,i], na.rm = T)
  cutoff2 = mean(d1[,i], na.rm = T) - 3.5*sd(d1[,i], na.rm = T)
  
  d1[,i][which(d1[,i] > cutoff1 | d1[,i] < cutoff2)] <- NA 
} 

# ### MAKE COVARIATES MISSINGNESS TABLE 

# my.list <- c("pack_years", "units", "usual", "bmi", "rank", "years")

# names <- list("Pack Years Smoked", "Units Alcohol", "Usual Alcohol", "BMI", "SIMD", "Educational Attainment")

# output <- matrix(nrow = 1*length(my.list), ncol = 3)
# output <- as.data.frame(output)
# j=c(1:length(my.list))

# for(i in 1:length(my.list)){
#   trait <- as.character(my.list[[i]])
#   df <- d1[,c(trait)]
#   df <- na.omit(df)
#   table <- length(df)
#   missing <- (dim(d1)[1]) - table
#   # output metrics 
#   output[j[i],1] <- names[i]
#   output[j[i],2] <- table
#   output[j[i],3] <- missing
#   names(output) <- c("Trait", "Present (n)", "Missing (n)")
# }

# # order by cases low to high
# sort <- output[order(output[,3]),]  

# write.csv(sort, file = "/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/covariate_missingness_table_GS_20k_090821.csv", row.names = F)

####################################################################################

# Add maximal covariate info from protein file
prot <- read.delim("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/file_input_030821/GS20K_GDF15_NT_proBNP.PHE")
names(prot)[2] <- "Sample_Name"
prot <- prot[-c(3:6)]
d1 <- left_join(d1, prot, by = "Sample_Name")

# Add in protein data as predictors - pre prepped in N used for analyses
GDF15 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_All_initial_files/GDF15_data_cox.csv")
NT <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_All_initial_files/BNP_data_cox.csv")

GDF15 <- GDF15[c(2,4)]
NT <- NT[c(2,6)]

d1 <- left_join(d1, GDF15, by = 'Sample_Name')
d1 <- left_join(d1, NT, by = 'Sample_Name')

table(is.na(d1$gdf15))
table(is.na(d1$nt.probnp))

# Set the clock for variables of interest - the usual clock for running all proteins - just to rank transformed first 
clock <- names(d1)[c(160,161)] 
write.csv(d1, "/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_All_initial_files/d1.csv", row.names = F)

####################################################################################

# Read in the prepped cases files for each trait to be assessed 

# Diabetes <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_All_initial_files/Diabetes_combined_consent.csv")

# Stroke <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_All_initial_files/Stroke_combined_consent.csv")

# IHD <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_All_initial_files/IHD_combined_consent.csv")

# DEM <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_All_initial_files/all_dementia_consent.csv")

Diabetes <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_Cox_proteins_20k/Trait_preps/Diabetes_combined.csv")

Stroke <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_Cox_proteins_20k/Trait_preps/Stroke_combined.csv")

IHD <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_Cox_proteins_20k/Trait_preps/IHD_combined.csv")

DEM <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_Cox_proteins_20k/Trait_preps/all_dementia.csv")

####################################################################################

## HAZARD MODELS - PREDICTING TIME-TO-ONSET OF DISEASES FROM STUDY BASELINE (2006)

####################################################################################

## Run Incidence of Alzheimer's Disease Analysis (aged equal to or over 65)

d1_DEM <- d1

mat_hazard_ad <- matrix(nrow=length(clock),ncol=9)
output_hazard_DEM<- as.data.frame(mat_hazard_ad)
for(j in 1:length(clock)){ 
  tryCatch({ 
  dat1= d1_DEM
  tmp1 = DEM[which(DEM$id %in% dat1$Sample_Name),] # 135 individuals 
  
  ## Obtain Age of Onset 
  affected = dat1[which(dat1$Sample_Name %in% tmp1$id),] 
  age_onset = DEM[,c("first", "id")]
  affected = merge(age_onset, affected, by.x = "id", by.y = "Sample_Name")
  affected$Event = 1
  affected$yoe = substring(affected$first, 1, 4)
  affected$moe = substring(affected$first, 5,6)
  affected$month_event1 = (as.numeric(affected$moe) - as.numeric(affected$mob))/12
  affected$age_event1 = as.numeric(affected$yoe) - as.numeric(affected$yob)
  affected$age_event = affected$age_event1 + affected$month_event1
  affected$first = NULL
  affected$yoe = NULL 
  affected$moe = NULL
  affected$month_event1 = NULL 
  affected$age_event1 = NULL
  
  healthy = dat1[-which(dat1$Sample_Name %in% DEM$id),]
  healthy$Event = 0
  healthy$age_event = 0 
  affected$id.y <- NULL
  healthy$id <- NULL
  names(affected)[names(affected)=="id"] <- "Sample_Name"
  cox = rbind(affected, healthy)
  
  cox$age_death = 0
  cox$age_death = ifelse(cox$dead %in% 1, cox$aged, 0)
  cox$age_at_event = ifelse(cox$Event %in% 1, cox$age_event, (ifelse(cox$dead %in% 1 & cox$Event %in% 0, cox$age_death, cox$aged)))
  cox$tte = cox$age_at_event - cox$age.x
  cox$tte = as.numeric(cox$tte)
  cox$tte <- ifelse(cox$tte < -1, "NA", cox$tte)
  cox$tte = ifelse(cox$tte < 0, 0, cox$tte)
  cox$Event = as.numeric(cox$Event)
  cox$tte<-as.numeric(cox$tte)

  cox <- cox[complete.cases(cox[,clock[[j]]]),]
  
  cox = cox[cox$age_at_event >=65,]
  mod = coxme(Surv(cox$tte, cox$Event) ~ scale(cox[,clock[[j]]]) + factor(cox$sex.x) + cox$age.x + (1|cox$Sample_Name), varlist = kin_model*2)
  print(clock[[j]])
  print("DEM")
  output_hazard_DEM[j,1] <- as.character(clock[[j]])
  output_hazard_DEM[j,2] <- as.character("Dementia")
  output_hazard_DEM[j,3:5]<-round(exp(cbind(coef(mod), confint(mod)))[1,1:3],2)
  output_hazard_DEM[j,6] <- extract_coxme_table(mod)[1,4]
  output_hazard_DEM[j,7] <- mod$n[1]
  output_hazard_DEM[j,8] <- mod$n[2]-mod$n[1]
  output_hazard_DEM[j,9] <-cox.zph(mod)[1][[1]][3]

  # Get time to event info for included cases, with max tte reported
  cox$Event <- ifelse(cox$tte < 0, "NA", cox$Event) # add this in to change event to NA as well as tte
  check <- cox %>% filter(Event == "1")
  check <- check %>% filter(!tte == 'NA')

  mean_tte <- mean(check$tte, na.rm = T) %>% round(digits = 1)
  sd_tte <- sd(check$tte, na.rm = T) %>% round(digits = 1) 
  mean_sd <- paste0(mean_tte, " ", "(", sd_tte, ")") 
  output_hazard_DEM[j,10] <- mean_sd

  max <- max(check$tte, na.rm = T)
  output_hazard_DEM[j,11] <- max

  }, error = function(e) cat("skipped"))
} 


## Create List of Remaining Dataframes

my.list = list(Stroke,Diabetes,IHD)
names = list("Stroke","Diabetes","IHD")
names(my.list) <- names 

l=lapply(my.list, "[", c(1:2))

names(d1)[names(d1) == "heart_disease_Y"] <- "IHD"
names(d1)[names(d1) == "stroke_Y"] <- "Stroke"
names(d1)[names(d1) == "diabetes_Y"] <- "Diabetes"

mat_hazard <- matrix(nrow=200*length(my.list),ncol=9)
output_hazard <- as.data.frame(mat_hazard)
k=c(0,200,400,600,800,1000,1200,1400)

## Loop of Survival Models - Longitudinal Associations
for(j in 1:length(clock)){
  for(i in 1:length(l)){ 
   tryCatch({ 
    tmp <- l[[i]]
    
    ## Exclude Indiviudals who Reported Disease at Study Baseline  
    dat1= d1[-which(d1[,names[[i]]] %in% 1),]
    tmp1 = tmp[which(tmp$id %in% dat1$Sample_Name),]
    
    ## Obtain Age of Onset 
    affected = dat1[which(dat1$Sample_Name %in% tmp1$id),] 
    age_onset = tmp[,c("first", "id")]
    affected = merge(age_onset, affected, by.x = "id", by.y = "Sample_Name")
    affected$Event = 1
    affected$yoe = substring(affected$first, 1, 4)
    affected$moe = substring(affected$first, 5,6)
    affected$month_event1 = (as.numeric(affected$moe) - as.numeric(affected$mob))/12
    affected$age_event1 = as.numeric(affected$yoe) - as.numeric(affected$yob)
    affected$age_event = affected$age_event1 + affected$month_event1
    affected$first = NULL
    affected$yoe = NULL 
    affected$moe = NULL
    affected$month_event1 = NULL 
    affected$age_event1 = NULL
    
    healthy = dat1[-which(dat1$Sample_Name %in% tmp$id),]
    healthy$Event = 0
    healthy$age_event = 0 
    affected$id.y <- NULL
    healthy$id <- NULL
    names(affected)[names(affected)=="id"] <- "Sample_Name"
    cox = rbind(affected, healthy)
    
    cox$age_death = 0
    cox$age_death = ifelse(cox$dead %in% 1, cox$aged, 0)
    cox$age_at_event = ifelse(cox$Event %in% 1, cox$age_event, (ifelse(cox$dead %in% 1 & cox$Event %in% 0, cox$age_death, cox$aged)))
    cox$tte = cox$age_at_event - cox$age.x
    cox$tte = as.numeric(cox$tte)
    cox$tte <- ifelse(cox$tte < -1, "NA", cox$tte)
    cox$tte = ifelse(cox$tte < 0, 0, cox$tte)
    cox$Event = as.numeric(cox$Event)
    cox$tte<-as.numeric(cox$tte)
       
    cox <- cox[complete.cases(cox[,clock[[j]]]),] 

    mod = coxme(Surv(cox$tte, cox$Event) ~ scale(cox[,clock[[j]]]) + cox$age.x + factor(cox$sex.x) + (1|cox$Sample_Name), varlist = kin_model*2)
    print(names[[i]])
    print(clock[[j]])
    output_hazard[j+k[[i]],1] <- as.character(clock[[j]])
    output_hazard[j+k[[i]],2] <- as.character(names[[i]])
    output_hazard[j+k[[i]],3:5] <- round(exp(cbind(coef(mod), confint(mod)))[1,1:3],2)
    output_hazard[j+k[[i]],6] <- extract_coxme_table(mod)[1,4]
    output_hazard[j+k[[i]],7] <- mod$n[1]
    output_hazard[j+k[[i]],8] <- mod$n[2]-mod$n[1]
    output_hazard[j+k[[i]],9] <-cox.zph(mod)[1][[1]][3]

    cox$Event <- ifelse(cox$tte < 0, "NA", cox$Event) # add this in to change event to NA as well as tte
    check <- cox %>% filter(Event == "1")
    check <- check %>% filter(!tte == 'NA')

    mean_tte <- mean(check$tte, na.rm = T) %>% round(digits = 1)
    sd_tte <- sd(check$tte, na.rm = T) %>% round(digits = 1) 
    mean_sd <- paste0(mean_tte, " ", "(", sd_tte, ")") 
    output_hazard[j+k[[i]],10]  <- mean_sd

    max <- max(check$tte, na.rm = T)
    output_hazard[j+k[[i]],11] <- max

   }, error = function(e) cat("skipped"))
  } 
} 

# Collate results and save 
comb <- rbind(output_hazard_DEM, output_hazard)
comb <- na.omit(comb)
names(comb) <- c("Predictor", "Outcome", "Hazard Ratio", "LCI", "UCI", "P.Value", "No. of Cases", "No. of Controls", "cox.zph", 'tte_mean_sd', 'max_tte')
write.csv(comb, "/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_All_initial_files/01_Cox_20k_results/basic.csv", row.names = F)

######################################################################################

### FULLY ADJUSTED MODELS 

d1_DEM <- d1

mat_hazard_ad <- matrix(nrow=length(clock),ncol=9)
output_hazard_DEM<- as.data.frame(mat_hazard_ad)
for(j in 1:length(clock)){ 
  tryCatch({ 
  dat1= d1_DEM
  tmp1 = DEM[which(DEM$id %in% dat1$Sample_Name),] # 135 individuals 
  
  ## Obtain Age of Onset 
  affected = dat1[which(dat1$Sample_Name %in% tmp1$id),] 
  age_onset = DEM[,c("first", "id")]
  affected = merge(age_onset, affected, by.x = "id", by.y = "Sample_Name")
  affected$Event = 1
  affected$yoe = substring(affected$first, 1, 4)
  affected$moe = substring(affected$first, 5,6)
  affected$month_event1 = (as.numeric(affected$moe) - as.numeric(affected$mob))/12
  affected$age_event1 = as.numeric(affected$yoe) - as.numeric(affected$yob)
  affected$age_event = affected$age_event1 + affected$month_event1
  affected$first = NULL
  affected$yoe = NULL 
  affected$moe = NULL
  affected$month_event1 = NULL 
  affected$age_event1 = NULL
  
  healthy = dat1[-which(dat1$Sample_Name %in% DEM$id),]
  healthy$Event = 0
  healthy$age_event = 0 
  affected$id.y <- NULL
  healthy$id <- NULL
  names(affected)[names(affected)=="id"] <- "Sample_Name"
  cox = rbind(affected, healthy)
  
  cox$age_death = 0
  cox$age_death = ifelse(cox$dead %in% 1, cox$aged, 0)
  cox$age_at_event = ifelse(cox$Event %in% 1, cox$age_event, (ifelse(cox$dead %in% 1 & cox$Event %in% 0, cox$age_death, cox$aged)))
  cox$tte = cox$age_at_event - cox$age.x
  cox$tte = as.numeric(cox$tte)
  cox$tte <- ifelse(cox$tte < -1, "NA", cox$tte)
  cox$tte = ifelse(cox$tte < 0, 0, cox$tte)
  cox$Event = as.numeric(cox$Event)
  cox$tte<-as.numeric(cox$tte)

  cox <- cox[complete.cases(cox[,clock[[j]]]),]
  
  cox = cox[cox$age_at_event >=65,]
  mod = coxme(Surv(cox$tte, cox$Event) ~ scale(cox[,clock[[j]]]) + factor(cox$sex.x) + cox$age.x + cox$units + cox$smokingScore + cox$rank + cox$years + cox$bmi + (1|cox$Sample_Name), varlist = kin_model*2)
  print(clock[[j]])
  print("DEM")
  output_hazard_DEM[j,1] <- as.character(clock[[j]])
  output_hazard_DEM[j,2] <- as.character("Dementia")
  output_hazard_DEM[j,3:5]<-round(exp(cbind(coef(mod), confint(mod)))[1,1:3],2)
  output_hazard_DEM[j,6] <- extract_coxme_table(mod)[1,4]
  output_hazard_DEM[j,7] <- mod$n[1]
  output_hazard_DEM[j,8] <- mod$n[2]-mod$n[1]
  output_hazard_DEM[j,9] <-cox.zph(mod)[1][[1]][3]

  cox$Event <- ifelse(cox$tte < 0, "NA", cox$Event) # add this in to change event to NA as well as tte
  check <- cox %>% filter(Event == "1")
  check <- check %>% filter(!tte == 'NA')

    # remove missing covariates 
    join <- check
    join <- join[!is.na(join$units), ]
    join <- join[!is.na(join$smokingScore), ]
    join <- join[!is.na(join$rank), ]
    join <- join[!is.na(join$years), ]
    join <- join[!is.na(join$bmi), ]

    check <- join

  mean_tte <- mean(check$tte, na.rm = T) %>% round(digits = 1)
  sd_tte <- sd(check$tte, na.rm = T) %>% round(digits = 1) 
  mean_sd <- paste0(mean_tte, " ", "(", sd_tte, ")") 
  output_hazard_DEM[j,10] <- mean_sd

  max <- max(check$tte, na.rm = T)
  output_hazard_DEM[j,11] <- max

  }, error = function(e) cat("skipped"))
} 


## Create List of Remaining Dataframes

my.list = list(Stroke,Diabetes,IHD)
names = list("Stroke","Diabetes","IHD")
names(my.list) <- names 

l=lapply(my.list, "[", c(1:2))

names(d1)[names(d1) == "heart_disease_Y"] <- "IHD"
names(d1)[names(d1) == "stroke_Y"] <- "Stroke"
names(d1)[names(d1) == "diabetes_Y"] <- "Diabetes"

mat_hazard <- matrix(nrow=200*length(my.list),ncol=9)
output_hazard <- as.data.frame(mat_hazard)
k=c(0,200,400,600,800,1000,1200,1400)

## Loop of Survival Models - Longitudinal Associations
for(j in 1:length(clock)){
  for(i in 1:length(l)){ 
   tryCatch({ 
    tmp <- l[[i]]
    
    ## Exclude Indiviudals who Reported Disease at Study Baseline  
    dat1= d1[-which(d1[,names[[i]]] %in% 1),]
    tmp1 = tmp[which(tmp$id %in% dat1$Sample_Name),]
    
    ## Obtain Age of Onset 
    affected = dat1[which(dat1$Sample_Name %in% tmp1$id),] 
    age_onset = tmp[,c("first", "id")]
    affected = merge(age_onset, affected, by.x = "id", by.y = "Sample_Name")
    affected$Event = 1
    affected$yoe = substring(affected$first, 1, 4)
    affected$moe = substring(affected$first, 5,6)
    affected$month_event1 = (as.numeric(affected$moe) - as.numeric(affected$mob))/12
    affected$age_event1 = as.numeric(affected$yoe) - as.numeric(affected$yob)
    affected$age_event = affected$age_event1 + affected$month_event1
    affected$first = NULL
    affected$yoe = NULL 
    affected$moe = NULL
    affected$month_event1 = NULL 
    affected$age_event1 = NULL
    
    healthy = dat1[-which(dat1$Sample_Name %in% tmp$id),]
    healthy$Event = 0
    healthy$age_event = 0 
    affected$id.y <- NULL
    healthy$id <- NULL
    names(affected)[names(affected)=="id"] <- "Sample_Name"
    cox = rbind(affected, healthy)
    
    cox$age_death = 0
    cox$age_death = ifelse(cox$dead %in% 1, cox$aged, 0)
    cox$age_at_event = ifelse(cox$Event %in% 1, cox$age_event, (ifelse(cox$dead %in% 1 & cox$Event %in% 0, cox$age_death, cox$aged)))
    cox$tte = cox$age_at_event - cox$age.x
    cox$tte = as.numeric(cox$tte)
    cox$tte <- ifelse(cox$tte < -1, "NA", cox$tte)
    cox$tte = ifelse(cox$tte < 0, 0, cox$tte)
    cox$Event = as.numeric(cox$Event)
    cox$tte<-as.numeric(cox$tte)
       
    cox <- cox[complete.cases(cox[,clock[[j]]]),]

    mod = coxme(Surv(cox$tte, cox$Event) ~ scale(cox[,clock[[j]]]) + factor(cox$sex.x) + cox$age.x + cox$units + cox$smokingScore + cox$rank + cox$years + cox$bmi + (1|cox$Sample_Name), varlist = kin_model*2)
    print(names[[i]])
    print(clock[[j]])
    output_hazard[j+k[[i]],1] <- as.character(clock[[j]])
    output_hazard[j+k[[i]],2] <- as.character(names[[i]])
    output_hazard[j+k[[i]],3:5] <- round(exp(cbind(coef(mod), confint(mod)))[1,1:3],2)
    output_hazard[j+k[[i]],6] <- extract_coxme_table(mod)[1,4]
    output_hazard[j+k[[i]],7] <- mod$n[1]
    output_hazard[j+k[[i]],8] <- mod$n[2]-mod$n[1]
    output_hazard[j+k[[i]],9] <-cox.zph(mod)[1][[1]][3]

    cox$Event <- ifelse(cox$tte < 0, "NA", cox$Event) # add this in to change event to NA as well as tte
    check <- cox %>% filter(Event == "1")
    check <- check %>% filter(!tte == 'NA')

    # remove missing covariates 
    join <- check
    join <- join[!is.na(join$units), ]
    join <- join[!is.na(join$smokingScore), ]
    join <- join[!is.na(join$rank), ]
    join <- join[!is.na(join$years), ]
    join <- join[!is.na(join$bmi), ]

    check <- join
    mean_tte <- mean(check$tte, na.rm = T) %>% round(digits = 1)
    sd_tte <- sd(check$tte, na.rm = T) %>% round(digits = 1) 
    mean_sd <- paste0(mean_tte, " ", "(", sd_tte, ")") 
    output_hazard[j+k[[i]],10]  <- mean_sd

    max <- max(check$tte, na.rm = T)
    output_hazard[j+k[[i]],11] <- max

   }, error = function(e) cat("skipped"))
  } 
} 


# Collate results and save 
comb <- rbind(output_hazard_DEM, output_hazard)
names(comb) <- c("Predictor", "Outcome", "Hazard Ratio", "LCI", "UCI", "P.Value", "No. of Cases", "No. of Controls", "cox.zph", 'tte_mean_sd')
comb <- na.omit(comb)
# Order by p value 
# comb <- comb[order(comb$`P.Value`),]
# Write ordered basic results 
write.csv(comb, "/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_All_initial_files/01_Cox_20k_results/full.csv", row.names = F)


# sessionInfo()

# R version 4.0.3 (2020-10-10)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: openSUSE Leap 15.1

# Matrix products: default
# BLAS:   /usr/local/lib64/R/lib/libRblas.so
# LAPACK: /usr/local/lib64/R/lib/libRlapack.so

# locale:
#  [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C
#  [3] LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8
#  [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8
#  [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C
#  [9] LC_ADDRESS=C               LC_TELEPHONE=C
# [11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C

# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base

# other attached packages:
#  [1] gdata_2.18.0        coxme_2.2-16        bdsmatrix_1.3-4
#  [4] kinship2_1.8.5      quadprog_1.5-8      Matrix_1.4-1
#  [7] survival_3.3-1      bestNormalize_1.8.2 forcats_0.5.1
# [10] stringr_1.5.0       dplyr_1.0.10        purrr_1.0.0
# [13] readr_2.1.3         tidyr_1.2.1         tibble_3.1.8
# [16] ggplot2_3.4.0       tidyverse_1.3.1     readxl_1.4.0

# loaded via a namespace (and not attached):
#  [1] httr_1.4.4          jsonlite_1.8.4      splines_4.0.3
#  [4] foreach_1.5.2       gtools_3.9.4        prodlim_2019.11.13
#  [7] modelr_0.1.8        assertthat_0.2.1    doRNG_1.8.2
# [10] cellranger_1.1.0    globals_0.16.2      ipred_0.9-13
# [13] pillar_1.8.1        backports_1.4.1     lattice_0.20-45
# [16] glue_1.6.2          digest_0.6.31       rvest_1.0.2
# [19] hardhat_1.2.0       colorspace_2.1-0    recipes_1.0.3
# [22] timeDate_4021.107   pkgconfig_2.0.3     broom_0.8.0
# [25] listenv_0.9.0       haven_2.5.0         scales_1.2.1
# [28] gower_1.0.0         lava_1.7.0          tzdb_0.3.0
# [31] timechange_0.1.1    generics_0.1.3      usethis_2.1.6
# [34] ellipsis_0.3.2      withr_2.5.0         nnet_7.3-17
# [37] cli_3.5.0           magrittr_2.0.3      crayon_1.5.2
# [40] fs_1.6.1            fansi_1.0.4         future_1.30.0
# [43] parallelly_1.33.0   nlme_3.1-157        doParallel_1.0.17
# [46] MASS_7.3-56         xml2_1.3.2          class_7.3-20
# [49] tools_4.0.3         hms_1.1.2           lifecycle_1.0.3
# [52] munsell_0.5.0       reprex_2.0.1        rngtools_1.5.2
# [55] compiler_4.0.3      rlang_1.0.6         grid_4.0.3
# [58] iterators_1.0.14    rstudioapi_0.14     gtable_0.3.1
# [61] codetools_0.2-18    DBI_1.1.3           R6_2.5.1
# [64] lubridate_1.9.0     future.apply_1.10.0 utf8_1.2.3
# [67] nortest_1.0-4       butcher_0.1.5       stringi_1.7.8
# [70] parallel_4.0.3      Rcpp_1.0.9          vctrs_0.5.1
# [73] rpart_4.1.16        dbplyr_2.1.1        tidyselect_1.2.0





