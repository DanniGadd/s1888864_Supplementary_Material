
####################################################################################

# 20k - Run coxph models and extract local coxzph for each protein in 25 associations

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

# Calculate prevalent AD (it was split into males and females)
d1$AD <- 0 
d1[d1$Sample_Name %in% d1[(d1$alzheimers_M %in% 1 | d1$alzheimers_F %in% 1),"Sample_Name"],"AD"] <- 1

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

# Add in protein data as predictors
prot <- read.delim("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/file_input_030821/GS20K_GDF15_NT_proBNP.PHE")
names(prot)[2] <- "Sample_Name"
d1 <- left_join(d1, prot, by = "Sample_Name")


####################################################################################

# Read in the prepped cases files for each trait to be assessed 

# dem <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_Cox_proteins_20k/Trait_preps/all_dementia.csv")

AD <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_Cox_proteins_20k/Trait_preps/AD.csv")

VD <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_Cox_proteins_20k/Trait_preps/VD.csv")

Diabetes <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_Cox_proteins_20k/Trait_preps/Diabetes_combined.csv")

Stroke <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_Cox_proteins_20k/Trait_preps/Stroke_combined.csv")

IHD <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_Cox_proteins_20k/Trait_preps/IHD_combined.csv")

DEM <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_Cox_proteins_20k/Trait_preps/all_dementia.csv")


####################################################################################

# Set the clock for variables of interest - the usual clock for running all proteins - just to rank transformed first 
clock <- names(d1)[c(139,141)] 

# Work out max sample for each protein used in the associations 
table(is.na(d1$gdf15_rnk))
table(is.na(d1$nt.probnp_rnk))

# > table(is.na(d1$gdf15_rnk))
# FALSE  TRUE
# 18414  5678
# > table(is.na(d1$nt.probnp_rnk))
# FALSE  TRUE
# 17863  6229

# Subset to just those with DNAm available
d1 <- d1[complete.cases(d1$smokingScore),]

table(is.na(d1$gdf15_rnk))
table(is.na(d1$nt.probnp_rnk))


# FALSE  TRUE
# 17489   924
# > table(is.na(d1$nt.probnp_rnk))

# FALSE  TRUE
# 16963  1450

# write.csv(d1, "/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_Cox_proteins_20k/d1.csv", row.names = F)



####################################################################################

# Run coxph models and extract local coxzph for each protein

####################################################################################


# Load in basic model IHD cox tables for GDF and BNP


# read in the cox tables for each trait

DEM_GDF <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_Cox_proteins_20k/Cox_proteins_210522/Cox_table_extraction_basic/Cox_DEM_gdf15_rnk.csv")
DEM_GDF <- DEM_GDF %>% select("Sample_Name", "age.x", "sex.x", "dead", "aged", "Event", "age_death", "age_at_event", "tte",
  "bmi", "units", "smokingScore", "rank", "years", "gdf15_rnk", "nt.probnp_rnk")

DEM_BNP <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_Cox_proteins_20k/Cox_proteins_210522/Cox_table_extraction_basic/Cox_DEM_nt.probnp_rnk.csv")
DEM_BNP <- DEM_BNP %>% select("Sample_Name", "age.x", "sex.x", "dead", "aged", "Event", "age_death", "age_at_event", "tte",
  "bmi", "units", "smokingScore", "rank", "years", "gdf15_rnk", "nt.probnp_rnk")

DI_GDF <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_Cox_proteins_20k/Cox_proteins_210522/Cox_table_extraction_basic/Cox_Diabetes_gdf15_rnk.csv")
DI_GDF <- DI_GDF %>% select("Sample_Name", "age.x", "sex.x", "dead", "aged", "Event", "age_death", "age_at_event", "tte",
  "bmi", "units", "smokingScore", "rank", "years", "gdf15_rnk", "nt.probnp_rnk")

DI_BNP <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_Cox_proteins_20k/Cox_proteins_210522/Cox_table_extraction_basic/Cox_Diabetes_nt.probnp_rnk.csv")
DI_BNP <- DI_BNP %>% select("Sample_Name", "age.x", "sex.x", "dead", "aged", "Event", "age_death", "age_at_event", "tte",
  "bmi", "units", "smokingScore", "rank", "years", "gdf15_rnk", "nt.probnp_rnk")

ST_GDF <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_Cox_proteins_20k/Cox_proteins_210522/Cox_table_extraction_basic/Cox_Stroke_gdf15_rnk.csv")
ST_GDF <- ST_GDF %>% select("Sample_Name", "age.x", "sex.x", "dead", "aged", "Event", "age_death", "age_at_event", "tte",
  "bmi", "units", "smokingScore", "rank", "years", "gdf15_rnk", "nt.probnp_rnk")

ST_BNP <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_Cox_proteins_20k/Cox_proteins_210522/Cox_table_extraction_basic/Cox_Stroke_nt.probnp_rnk.csv")
ST_BNP <- ST_BNP %>% select("Sample_Name", "age.x", "sex.x", "dead", "aged", "Event", "age_death", "age_at_event", "tte",
  "bmi", "units", "smokingScore", "rank", "years", "gdf15_rnk", "nt.probnp_rnk")

IHD_GDF <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_Cox_proteins_20k/Cox_proteins_210522/Cox_table_extraction_basic/Cox_IHD_gdf15_rnk.csv")
IHD_GDF <- IHD_GDF %>% select("Sample_Name", "age.x", "sex.x", "dead", "aged", "Event", "age_death", "age_at_event", "tte",
  "bmi", "units", "smokingScore", "rank", "years", "gdf15_rnk", "nt.probnp_rnk")

IHD_BNP <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_Cox_proteins_20k/Cox_proteins_210522/Cox_table_extraction_basic/Cox_IHD_nt.probnp_rnk.csv")
IHD_BNP <- IHD_BNP %>% select("Sample_Name", "age.x", "sex.x", "dead", "aged", "Event", "age_death", "age_at_event", "tte",
  "bmi", "units", "smokingScore", "rank", "years", "gdf15_rnk", "nt.probnp_rnk")



# Run coxph 

assoc <- data.frame(Predictor = c("gdf15_rnk", "gdf15_rnk", "gdf15_rnk", "gdf15_rnk", "nt.probnp_rnk", "nt.probnp_rnk", "nt.probnp_rnk", "nt.probnp_rnk"), 
	Outcome = c("DI_GDF", "ST_GDF", "IHD_GDF", "DEM_GDF", "ST_BNP", "DI_BNP", "IHD_BNP", "DEM_BNP"))
names1 <- assoc$Predictor
names2 <- assoc$Outcome

# make a repository for the data 
results <- data.frame(predictor = names1, trait = names2, basic_cox.zph = 1:8, basic_global = 1:8)

# make list of possible data tables 
list <- list(DEM_GDF, DEM_BNP, DI_GDF, DI_BNP, ST_GDF, ST_BNP, IHD_GDF, IHD_BNP)
names(list) <- c("DEM_GDF", "DEM_BNP", "DI_GDF", "DI_BNP", "ST_GDF", "ST_BNP", "IHD_GDF", "IHD_BNP")

for (i in 1:8){
	outcome <- as.character(assoc[i,2]) # get the disease outcome of interest 
	data <- list[[outcome]] # extract the cox info for the disease of interest from the list 
	# data <- left_join(data, d1, by = "Sample_Name") # join in the phneotypic information to those in the disease set 
	# Getting the protein name 
	protein <- as.character(assoc[i,1])
	cox <- data[c("tte", "Event", protein, "sex.x", "age.x", "Sample_Name")]
	names(cox)[3] <- "prot"

	# cox$filter <- ifelse(cox$time_to_event <= 5, 0, 1)
	# cox$new_tte <- ifelse(cox$time_to_event > 5, 5, cox$time_to_event)
	# cox$new_event <- ifelse(cox$time_to_event > 5, 0, cox$Event)

	# #cox  <- cox[cox$time_to_event < 10,]
	# for(j in 1:16){
	# 	cox$filter <- ifelse(cox$time_to_event <= j, 0, 1)
	# 	cox$new_tte <- ifelse(cox$time_to_event > j, j, cox$time_to_event)
	# 	cox$new_event <- ifelse(cox$time_to_event > j, 0, cox$Event)
	# }

	# Run basic model
	mod1 = coxph(Surv(tte, Event) ~ scale(cox$prot) + cox$age.x + factor(cox$sex.x), data = cox)
	# mod1 = coxph(Surv(time_to_event, Event) ~ scale(cox$prot)*cox$filter + cox$age.x + factor(cox$sex.x), data = cox)
	# Extract results 
	all <- cox.zph(mod1) 

	# #cox  <- cox[cox$time_to_event < 10,]
	# # Run basic model
	# mod1 = coxph(Surv(new_tte, new_event) ~ scale(cox$prot) + cox$Age.x + factor(cox$Female.x)  + factor(cox$Set), data = cox)
	# #mod1 = coxph(Surv(time_to_event, Event) ~ scale(cox$prot)*cox$filter + cox$Age.x + factor(cox$Female.x)  + factor(cox$Set), data = cox)
	# # Extract results 
	# all <- cox.zph(mod1) 

	p <- all$table[,"p"]
	local <- p[1]
	global <- p[4]
	# Save into results table
	results[i,3] <- local
	results[i,4] <- global
	# # Run fully ajd model
	# mod1 = coxph(Surv(time_to_event, Event) ~ scale(cox$prot) + factor(cox$Female.x) + cox$Age.x + cox$units + factor(cox$usual) + cox$smokingScore + cox$simd + cox$EA + cox$bmi + factor(cox$Set), data = cox)
	# # Extract results 
	# all1 <- cox.zph(mod1) 
	# p1 <- all1$table[,"p"]
	# local1 <- p1[1]
	# global1 <- p1[11]
	# # Save into results table
	# results[i,5] <- local1
	# results[i,6] <- global1
}


write.csv(results, "/Local_Scratch/Danni/GDFBNP/03_Cox/Cox_sensitivities/20k_basic_coxzph.csv", row.names = F)




##########

### sensitivity analysis - for failures of coxph in 20k

# Run a loop for each association to extract 1-12 options and plot 
library(survival)
library(coxme)
library(tidyverse)
library(survivalAnalysis)

# make a repository for the data 
results <- data.frame(predictor = "X", trait = "X", TTFU = 1:300, HR = 1:300, LC = 1:300, UC = 1:300, P = 1:300, control = 1:300, case = 1:300, local_cox.zph = 1:300, global_cox.zph = 1:300)
results$predictor <- as.character(results$predictor)
results$trait <- as.character(results$trait)
timer <- seq(0,300, by = 12)

# make list of possible data tables 
list <- list(IHD_GDF, IHD_BNP)
names(list) <- c("IHD", "IHD")
list2 <- c("IHD", "IHD")
prots <- c("gdf15_rnk", "nt.probnp_rnk")

# Run loop
for (i in 1:2){
	k <- timer[[i]]
	outcome <- as.character(list2[i]) # get the disease outcome of interest 
	data <- list[[outcome]] # extract the cox table info for the disease of interest from the list 
	# Getting the protein name 
	protein <- as.character(prots[i])
	cox <- data[c("tte", "Event", protein, "sex.x", "age.x", "Sample_Name")]
	names(cox)[3] <- "prot"

	for(j in 1:12){
		cox$new_tte <- ifelse(cox$tte < j, cox$tte, j)
		cox$new_ev <- ifelse(cox$tte < j, cox$Event, 0)
		# Run basic model
		mod = coxph(Surv(new_tte, new_ev) ~ scale(prot) + factor(sex.x) + age.x, data = cox)
		results[j+k,1] <- as.character(protein)
  		results[j+k,2] <- as.character(outcome)
  		results[j+k,3] <- j
		results[j+k,4:6]<-round(exp(cbind(coef(mod), confint(mod)))[1,1:3],2)
		results[j+k,7] <- cox_as_data_frame(mod)[1,10]
		results[j+k,8] <- mod$n[1] - mod$nevent[1]
		results[j+k,9] <- mod$nevent[1]
		table <- cox.zph(mod)
		p1 <- table$table[,"p"]
		results[j+k,10] <-p1[1]
		results[j+k,11] <-p1[4]
	}
}

results <- results[1:24,]

write.csv(results, "/Local_Scratch/Danni/GDFBNP/03_Cox/Cox_sensitivities/20k_loops.csv", row.names = F)


#########################################################################################################################

### DO FOR W4 subset 



# read in the cox tables for each trait

DEM_GDF <- read.csv("/Local_Scratch/Danni/GDFBNP/03_Cox/Cox_w1w3_w4_220722/Tables/Cox_DEM_gdf15_rnk_EPIC.csv")
DEM_GDF <- DEM_GDF %>% select("Sample_Name", "age.x", "sex.x", "dead", "aged", "Event", "age_death", "age_at_event", "tte",
  "bmi", "units", "smokingScore", "rank", "years", "gdf15_rnk", "nt.probnp_rnk", "GDF_score", "BNP_score")

DEM_BNP <- read.csv("/Local_Scratch/Danni/GDFBNP/03_Cox/Cox_w1w3_w4_220722/Tables/Cox_DEM_BNP_score_EPIC.csv")
DEM_BNP <- DEM_BNP %>% select("Sample_Name", "age.x", "sex.x", "dead", "aged", "Event", "age_death", "age_at_event", "tte",
  "bmi", "units", "smokingScore", "rank", "years", "gdf15_rnk", "nt.probnp_rnk", "GDF_score", "BNP_score")

DI_GDF <- read.csv("/Local_Scratch/Danni/GDFBNP/03_Cox/Cox_w1w3_w4_220722/Tables/Cox_Diabetes_gdf15_rnk_EPIC.csv")
DI_GDF <- DI_GDF %>% select("Sample_Name", "age.x", "sex.x", "dead", "aged", "Event", "age_death", "age_at_event", "tte",
  "bmi", "units", "smokingScore", "rank", "years", "gdf15_rnk", "nt.probnp_rnk", "GDF_score", "BNP_score")

DI_BNP <- read.csv("/Local_Scratch/Danni/GDFBNP/03_Cox/Cox_w1w3_w4_220722/Tables/Cox_Diabetes_BNP_score_EPIC.csv")
DI_BNP <- DI_BNP %>% select("Sample_Name", "age.x", "sex.x", "dead", "aged", "Event", "age_death", "age_at_event", "tte",
  "bmi", "units", "smokingScore", "rank", "years", "gdf15_rnk", "nt.probnp_rnk", "GDF_score", "BNP_score")

ST_GDF <- read.csv("/Local_Scratch/Danni/GDFBNP/03_Cox/Cox_w1w3_w4_220722/Tables/Cox_Stroke_gdf15_rnk_EPIC.csv")
ST_GDF <- ST_GDF %>% select("Sample_Name", "age.x", "sex.x", "dead", "aged", "Event", "age_death", "age_at_event", "tte",
  "bmi", "units", "smokingScore", "rank", "years", "gdf15_rnk", "nt.probnp_rnk", "GDF_score", "BNP_score")

ST_BNP <- read.csv("/Local_Scratch/Danni/GDFBNP/03_Cox/Cox_w1w3_w4_220722/Tables/Cox_Stroke_BNP_score_EPIC.csv")
ST_BNP <- ST_BNP %>% select("Sample_Name", "age.x", "sex.x", "dead", "aged", "Event", "age_death", "age_at_event", "tte",
  "bmi", "units", "smokingScore", "rank", "years", "gdf15_rnk", "nt.probnp_rnk", "GDF_score", "BNP_score")

IHD_GDF <- read.csv("/Local_Scratch/Danni/GDFBNP/03_Cox/Cox_w1w3_w4_220722/Tables/Cox_IHD_gdf15_rnk_EPIC.csv")
IHD_GDF <- IHD_GDF %>% select("Sample_Name", "age.x", "sex.x", "dead", "aged", "Event", "age_death", "age_at_event", "tte",
  "bmi", "units", "smokingScore", "rank", "years", "gdf15_rnk", "nt.probnp_rnk", "GDF_score", "BNP_score")

IHD_BNP <- read.csv("/Local_Scratch/Danni/GDFBNP/03_Cox/Cox_w1w3_w4_220722/Tables/Cox_IHD_BNP_score_EPIC.csv")
IHD_BNP <- IHD_BNP %>% select("Sample_Name", "age.x", "sex.x", "dead", "aged", "Event", "age_death", "age_at_event", "tte",
  "bmi", "units", "smokingScore", "rank", "years", "gdf15_rnk", "nt.probnp_rnk", "GDF_score", "BNP_score")




# Run coxph 

assoc <- read.csv("/Local_Scratch/Danni/GDFBNP/03_Cox/Cox_w1w3_w4_220722/Results/basic_EPIC.csv")
names1 <- assoc$Predictor

assoc$Outcome <- c("DI_GDF", "ST_GDF", "DI_GDF", "ST_GDF", "DI_BNP", "DI_BNP", "DEM_GDF", "DEM_GDF", "IHD_GDF", "ST_BNP", "ST_BNP", "IHD_GDF", "IHD_BNP", "DEM_BNP", "IHD_BNP", "DEM_BNP")
names2 <- assoc$Outcome

# make a repository for the data 
results <- data.frame(predictor = names1, trait = names2, basic_cox.zph = 1:16, basic_global = 1:16)

# make list of possible data tables 
list <- list(DEM_GDF, DEM_BNP, DI_GDF, DI_BNP, ST_GDF, ST_BNP, IHD_GDF, IHD_BNP)
names(list) <- c("DEM_GDF", "DEM_BNP", "DI_GDF", "DI_BNP", "ST_GDF", "ST_BNP", "IHD_GDF", "IHD_BNP")


for (i in 1:16){
	outcome <- as.character(assoc[i,2]) # get the disease outcome of interest 
	data <- list[[outcome]] # extract the cox info for the disease of interest from the list 
	# data <- left_join(data, d1, by = "Sample_Name") # join in the phneotypic information to those in the disease set 
	# Getting the protein name 
	protein <- as.character(assoc[i,1])
	cox <- data[c("tte", "Event", protein, "sex.x", "age.x", "Sample_Name")]
	names(cox)[3] <- "prot"

	# cox$filter <- ifelse(cox$time_to_event <= 5, 0, 1)
	# cox$new_tte <- ifelse(cox$time_to_event > 5, 5, cox$time_to_event)
	# cox$new_event <- ifelse(cox$time_to_event > 5, 0, cox$Event)

	# #cox  <- cox[cox$time_to_event < 10,]
	# for(j in 1:16){
	# 	cox$filter <- ifelse(cox$time_to_event <= j, 0, 1)
	# 	cox$new_tte <- ifelse(cox$time_to_event > j, j, cox$time_to_event)
	# 	cox$new_event <- ifelse(cox$time_to_event > j, 0, cox$Event)
	# }

	# Run basic model
	mod1 = coxph(Surv(tte, Event) ~ scale(cox$prot) + cox$age.x + factor(cox$sex.x), data = cox)
	# mod1 = coxph(Surv(time_to_event, Event) ~ scale(cox$prot)*cox$filter + cox$age.x + factor(cox$sex.x), data = cox)
	# Extract results 
	all <- cox.zph(mod1) 

	# #cox  <- cox[cox$time_to_event < 10,]
	# # Run basic model
	# mod1 = coxph(Surv(new_tte, new_event) ~ scale(cox$prot) + cox$Age.x + factor(cox$Female.x)  + factor(cox$Set), data = cox)
	# #mod1 = coxph(Surv(time_to_event, Event) ~ scale(cox$prot)*cox$filter + cox$Age.x + factor(cox$Female.x)  + factor(cox$Set), data = cox)
	# # Extract results 
	# all <- cox.zph(mod1) 

	p <- all$table[,"p"]
	local <- p[1]
	global <- p[4]
	# Save into results table
	results[i,3] <- local
	results[i,4] <- global
	# # Run fully ajd model
	# mod1 = coxph(Surv(time_to_event, Event) ~ scale(cox$prot) + factor(cox$Female.x) + cox$Age.x + cox$units + factor(cox$usual) + cox$smokingScore + cox$simd + cox$EA + cox$bmi + factor(cox$Set), data = cox)
	# # Extract results 
	# all1 <- cox.zph(mod1) 
	# p1 <- all1$table[,"p"]
	# local1 <- p1[1]
	# global1 <- p1[11]
	# # Save into results table
	# results[i,5] <- local1
	# results[i,6] <- global1
}


write.csv(results, "/Local_Scratch/Danni/GDFBNP/03_Cox/Cox_sensitivities/W4_basic_coxzph.csv", row.names = F)


###################################################

### sensitivity analysis 

# Load in basic model IHD cox tables for GDF and BNP

# Run a loop for each association to extract 1-12 options and plot 
library(survival)
library(coxme)
library(tidyverse)
library(survivalAnalysis)

# make a repository for the data 
results <- data.frame(predictor = "X", trait = "X", TTFU = 1:300, HR = 1:300, LC = 1:300, UC = 1:300, P = 1:300, control = 1:300, case = 1:300, local_cox.zph = 1:300, global_cox.zph = 1:300)
results$predictor <- as.character(results$predictor)
results$trait <- as.character(results$trait)
timer <- seq(0,300, by = 12)

# make list of possible data tables 
list <- list(IHD_BNP, ST_BNP)
names(list) <- c("IHD", "ST")
list2 <- c("IHD", "ST")
prots <- c("nt.probnp_rnk", "nt.probnp_rnk")

# Run loop
for (i in 1:2){
	k <- timer[[i]]
	outcome <- as.character(list2[i]) # get the disease outcome of interest 
	data <- list[[outcome]] # extract the cox table info for the disease of interest from the list 
	# Getting the protein name 
	protein <- as.character(prots[i])
	cox <- data[c("tte", "Event", protein, "sex.x", "age.x", "Sample_Name", "units", "smokingScore", "rank", "years", "bmi")]
	names(cox)[3] <- "prot"

	for(j in 1:12){
		cox$new_tte <- ifelse(cox$tte < j, cox$tte, j)
		cox$new_ev <- ifelse(cox$tte < j, cox$Event, 0)
		# Run basic model
		mod = coxph(Surv(new_tte, new_ev) ~ scale(prot) + factor(sex.x) + age.x, data = cox)
		results[j+k,1] <- as.character(protein)
  		results[j+k,2] <- as.character(outcome)
  		results[j+k,3] <- j
		results[j+k,4:6]<-round(exp(cbind(coef(mod), confint(mod)))[1,1:3],2)
		results[j+k,7] <- cox_as_data_frame(mod)[1,10]
		results[j+k,8] <- mod$n[1] - mod$nevent[1]
		results[j+k,9] <- mod$nevent[1]
		table <- cox.zph(mod)
		p1 <- table$table[,"p"]
		results[j+k,10] <-p1[1]
		results[j+k,11] <-p1[4]
	}
}

results <- results[1:24,]

write.csv(results, "/Local_Scratch/Danni/GDFBNP/03_Cox/Cox_sensitivities/W4_loops.csv", row.names = F)









































