######################################################################################################

### TRAITS IN GS - DATA LINKAGE

######################################################################################################

# •	Codes are taken from here supplementary info of this paper 1.  
# •	Read codes were pulled in from consenting GPs, ICD10 codes were pulled in from smr (hospitalisation data) and morbidity records.
# •	There's a broad dementia flag (dementia) as well as subtypes (ad, ftd, vd, and dlb), based on codes in the supplementary table above.
# •	For each flag, there's a dt1 variable, which corresponds to the first occurrence of the related code
# •	Table is here:
# Z:\Generation_Scotland_data_Sep2021\Dementia_case_update_24Jan2022.csv

cd /Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/

screen

R

library(readxl)
library(tidyverse)

# Dementia 

# Read in dementia codes from daniel 
dem <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Dementia_GS/00_Data_exploration/Dementia_case_update_24Jan2022.csv")

# Summarise subgroup counts
table(dem$dementia) # 583
table(dem$ad) # 188
table(dem$vd) # 158
table(dem$dlb) # 1
table(dem$ftd) # 7

# Check that all GP codes have given consent for GP sharing 
consent <- read_excel("/Cluster_Filespace/Marioni_Group/Danni/Dementia_GS/00_Data_exploration/2020-09-02 GP consent ids (1).xlsx")
consent <- as.data.frame(consent)
consent$Consent <- "YES"
consent <- consent[c(1,5)]
dem <- left_join(dem, consent, by = "id")
gp <- dem[which(dem$source == "gp"),]
table(gp$Consent) # yes, all GP cases consented

# Process dementis trait - remove dlb and ftd
all <- dem[which(dem$dementia == "1"),]
# all <- all[-which(all$dlb == '1' | all$ftd == '1'),]

all <- all[c(1,10)]
names(all) <- c("id", "first")
all <- all[order(all$first),] # Check that the file is ordered by date (lowest to highest)
all <- all[which(all$first > 199001),] # Remove codes pre 199001
all2 <- all[-which(duplicated(all$id)),] # Remove any duplicate ID codes by taking the earliest date for that person (if ordered this should take first instance)
dim(all2) # 318
write.csv(all2, "/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_Cox_proteins_20k/Trait_preps/all_dementia.csv", row.names = F)

## Stroke 
Stroke <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/00_Rerun/Cox/Input_tables/Run_by_Danni/Disease_files/Stroke.csv")
Stroke_2 <- read_excel("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/00_Rerun/Cox/Input_tables/Run_by_Danni/Disease_files_GP/Stroke.xlsx", col_types = c("text", "text", "text", "date", "text"))
Stroke_2 <- Stroke_2[-grep("Injury|injury", Stroke_2$description),]
Stroke_2 <- Stroke_2[-grep("Mitochondrial", Stroke_2$description),]
Stroke_2 <- Stroke_2[-grep("Personal", Stroke_2$description),]
Stroke_2 <- Stroke_2[,c(5,4)]
names(Stroke_2) <- c("id", "first")
Stroke_2$first <- paste0(substring(Stroke_2$first, 1, 4), substring(Stroke_2$first, 6,7))
Stroke <- Stroke[,c(1,2)]
Stroke = rbind(Stroke, Stroke_2)
Stroke <- Stroke[order(Stroke$first),]
Stroke <- Stroke[which(Stroke$first > 199001),]
Stroke = Stroke[-which(duplicated(Stroke$id)),]
write.csv(Stroke, "/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_Cox_proteins_20k/Trait_preps/Stroke_combined.csv", row.names = F)


## Heart Disease 
IHD <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/00_Rerun/Cox/Input_tables/Run_by_Danni/Disease_files/Heart.csv")
IHD_2 <- read_excel("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/00_Rerun/Cox/Input_tables/Run_by_Danni/Disease_files_GP/IHD.xlsx", col_types = c("text", "text", "text", "date", "text"))
IHD_2 <- IHD_2[,c(5,4)]
names(IHD_2) <- c("id", "first")
IHD_2$first <- paste0(substring(IHD_2$first, 1, 4), substring(IHD_2$first, 6,7))
IHD <- IHD[,c(1,2)]
IHD = rbind(IHD, IHD_2)
IHD <- IHD[order(IHD$first),]
IHD <- IHD[which(IHD$first > 199001),]
IHD = IHD[-which(duplicated(IHD$id)),]
write.csv(IHD, "/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_Cox_proteins_20k/Trait_preps/IHD_combined.csv", row.names = F)


# Diabetes
Diabetes <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/GS_COVID/Daniel_diab_comparison/2022-15-03_diabetes_smr_update.csv")
Diabetes <- Diabetes[-grep("MODY|1|Other|Unknown", Diabetes$tname),] # remove type 1 and other cases
Diabetes <- na.omit(Diabetes)
names(Diabetes)[4] <- "first"
Diabetes_2 <- read_excel("/Cluster_Filespace/Marioni_Group/Danni/GS_COVID/Diabetes_update_codes_yipeng/Diabetes.xlsx", col_types = c("text", "text", "text", "date", "text")) 
Diabetes_2 <- as.data.frame(Diabetes_2)
Diabetes_2 <- Diabetes_2[-grep("Juvenile|juvenile|1|one|One|Steroid|steriod|glycosuria", Diabetes_2$description),]
Diabetes_2 <- Diabetes_2[-grep("Insulin dependent|Insulin-dependent|Insulin treated|Haemochromatosis", Diabetes_2$description),]
Diabetes_2 <- Diabetes_2[-which(Diabetes_2$description %in% c("[X]Insulin and oral hypoglycaemic [antidiabetic] drugs causing adverse effects in therapeutic use", "[V]Dietary surveillance and counselling")),]
Diabetes_2 <- Diabetes_2[-which(Diabetes_2$description %in% c("Diabetes mellitus induced by steroids")),]
Diabetes_2 <- Diabetes_2[,c(5,2,3,4)]
names(Diabetes_2) <- c("id", "code", "diag", "first")
Diabetes_2$first <- paste0(substring(Diabetes_2$first, 1, 4), substring(Diabetes_2$first, 6,7))
names(Diabetes) <- names(Diabetes_2)
Diabetes = rbind(Diabetes, Diabetes_2)
Diabetes <- Diabetes[order(Diabetes$first),]
Diabetes <- Diabetes[-which(Diabetes$first <= 199001),]
Diabetes <- Diabetes[-which(Diabetes$first == "NANA"),]
Diabetes <- Diabetes[-which(duplicated(Diabetes$id)),] # 1652
Diabetes <- Diabetes[c(1,4)]
write.csv(Diabetes, "/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_Cox_proteins_20k/Trait_preps/Diabetes_combined.csv", row.names = F)




# # Dementia 

# # Read in dementia codes from daniel 
# dem <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/Dementia_case_update_18Jan2023.csv")

# # Summarise subgroup counts
# table(dem$dementia) # 752
# table(dem$ad) # 268
# table(dem$vd) # 214
# table(dem$dlb) # 3
# table(dem$ftd) # 7

# # Consent
# consent <- read_excel("/Cluster_Filespace/Marioni_Group/Danni/Dementia_GS/00_Data_exploration/2020-09-02 GP consent ids (1).xlsx")
# consent <- as.data.frame(consent)
# consent$Consent <- "YES"
# consent <- consent[c(1,5)]
# dem <- left_join(dem, consent, by = "id")
# gp <- dem[which(dem$source == "gp"),]
# table(gp$Consent) # 89 of 249 entries

# # Process dementis trait - remove dlb and ftd
# all <- dem[which(dem$dementia == "1"),]
# # all <- all[-which(all$dlb == '1' | all$ftd == '1'),]

# # Summarise subgroup counts after exclusions
# table(all$dementia) # 742
# table(all$ad) # 268
# table(all$vd) # 214
# table(all$dlb) # 0
# table(all$ftd) # 0

# all <- all[c(2,11)]
# names(all) <- c("id", "first")
# all <- all[order(all$first),] # Check that the file is ordered by date (lowest to highest)
# all <- all[which(all$first > 199001),] # Remove codes pre 199001
# all2 <- all[-which(duplicated(all$id)),] # Remove any duplicate ID codes by taking the earliest date for that person (if ordered this should take first instance)
# dim(all2) # 373
# write.csv(all2, "/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_Cox_proteins_20k/Trait_preps/all_dementia.csv", row.names = F)


# ######################################################################################################

# ### DEMENTIA GS 

# ######################################################################################################

# # •	Codes are taken from here supplementary info of this paper 1.  
# # •	Read codes were pulled in from consenting GPs, ICD10 codes were pulled in from smr (hospitalisation data) and morbidity records.
# # •	There's a broad dementia flag (dementia) as well as subtypes (ad, ftd, vd, and dlb), based on codes in the supplementary table above.
# # •	For each flag, there's a dt1 variable, which corresponds to the first occurrence of the related code
# # •	Table is here:
# # Z:\Generation_Scotland_data_Sep2021\Dementia_case_update_24Jan2022.csv

# screen

# R

# library(readxl)
# library(tidyverse)

# # Read in dementia codes from daniel 
# dem <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Dementia_GS/00_Data_exploration/Dementia_case_update_24Jan2022.csv")

# # Will need to ask daniel about code and description inclusion to remove cases that do not fall in range 

# # Summarise subgroup counts
# table(dem$dementia) # 583
# table(dem$ad) # 188
# table(dem$vd) # 158
# table(dem$dlb) # 1
# table(dem$ftd) # 7

# # Check that all GP codes have given consent for GP sharing 
# consent <- read_excel("/Cluster_Filespace/Marioni_Group/Danni/Dementia_GS/00_Data_exploration/2020-09-02 GP consent ids (1).xlsx")
# consent <- as.data.frame(consent)
# consent$Consent <- "YES"
# consent <- consent[c(1,5)]
# dem <- left_join(dem, consent, by = "id")
# gp <- dem[which(dem$source == "gp"),]
# table(gp$Consent) # yes, all GP cases consented

# ######################################################################################################

# ## Check prevalent vs incident cases for each group 

# death <- read_excel("/Cluster_Filespace/Marioni_Group/Danni/Dementia_GS/00_Data_exploration/2021-04-19 dead Feb 2021.xlsx")
# death <- as.data.frame(death)
# death$died <- 1

# # Join to main info dataset 
# affected <- left_join(dem, death, by = "id") # 276 

# # Add mobyob info for each person 
# yob <- read.csv("/Cluster_Filespace/Marioni_Group/GS/GS_dataset/agesex_yobmob.csv")
# affected <- left_join(affected, yob, by = "id")

# # Work out which were prevalent and which were incident cases for each trait of interest with enough cases 
# # Dementia (all), dementia (AD) and dementia (vascular)

# # Dementia (all)
# affected$Event = 1
# affected$yoe = substring(affected$dt1_dementia, 1, 4)
# affected$moe = substring(affected$dt1_dementia, 5,6)
# affected$month_event1 = (as.numeric(affected$moe) - as.numeric(affected$mob))/12
# affected$age_event1 = as.numeric(affected$yoe) - as.numeric(affected$yob)
# affected$age_event = affected$age_event1 + affected$month_event1

# affected$age_death = 0
# affected$age_death = ifelse(affected$died %in% 1, affected$aged, 0)
# affected$age_at_event = ifelse(affected$Event %in% 1, affected$age_event, (ifelse(affected$died %in% 1 & affected$Event %in% 0, affected$age_death, affected$aged)))
# affected$tte = affected$age_at_event - affected$age
# affected$tte = as.numeric(affected$tte)
# affected$Event <- as.numeric(affected$Event)
# affected$incident_tte <- ifelse(affected$tte < -1, "NA", affected$tte)
# #affected$incident_tte = as.numeric(affected$incident_tte)
# #affected$tte = ifelse(cox$tte < 0, 0, cox$tte)
# #cox$Event = as.numeric(cox$Event)
# #cox$tte<-as.numeric(cox$tte)

# # How many were diagnosed pre baseline
# length(which(affected$incident_tte == "NA")) # 18 
# length(which(!affected$incident_tte == "NA")) # 544
# dim(affected) # 586 - of 586 dementia cases, 18 prevalent and 544 incident 

# # Filter to just those aged over 65 when diagnosed
# affected <- affected[which(affected$age_at_event >= 65),] # of 586 cases, 531 were diagnosed either at or over age 65

# # Check number of incident/prevalent after this filtering step > 65 years old at diagnosis
# length(which(affected$incident_tte == "NA")) # 13
# length(which(!affected$incident_tte == "NA")) # 518

# # So, of 531 cases >= 65 at diagnosis, we have 13 prevalent and 518 incident cases 

# # Check to see how many incident cases are AD and vascular subtypes 
# AD <- affected[which(affected$ad == "1"),]
# VD <- affected[which(affected$vd == "1"),]

# dim(AD)
# # [1] 176  32
# dim(VD)
# # [1] 153  32

# AD <- AD[unique(AD$id),] # 130 
# VD <- VD[unique(VD$id),] # 111

# AD$incident_tte <- as.numeric(AD$incident_tte)

# IDs <- unique(AD$id)

# mean(AD$incident_tte, na.rm = T) # mean of 8.7 years 

# ######################################################################################################

# ## PREP COX FILES 

# # Load episcores from wave 1/3 projections in the original episcores paper (taken from script and files at: Y:\Protein_DNAm_Proxies\Work_and_code_post_KORA\Gitlab_code_files\03_Projections\Projections_with_weights)
# episcores <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Dementia_GS/01_EpiScore_projections_GS_W4/EpiScore_projections_GS_9537.csv", check.names = F)

# # Load W4 projections calculated for this analaysis 
# target <- readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/wave4/w4-samplesheet.rds")
# W4 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Dementia_GS/01_EpiScore_projections_GS_W4/EpiScore_projections_W4_V2.csv", check.names = F)
# names(W4)[1] <- "ID"

# # Rank Inverse Based Normalisation of the data
# for(i in names(episcores)[2:110]){ 
#  episcores[,i] <- qnorm((rank(episcores[,i],na.last="keep")-0.5)/sum(!is.na(episcores[,i])))
# }

# for(i in names(W4)[2:110]){ 
#  episcores[,i] <- qnorm((rank(episcores[,i],na.last="keep")-0.5)/sum(!is.na(episcores[,i])))
# }


# # Join episcore projections to full set of 20k DNAm
# episcores <- rbind(episcores, W4)

# # Read in age information (mortality doc 2019)
# age_alive <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Proxies_stradl/KORA_train_recieved_020221/Cox/age_may19.csv")

# # Work out those who are alive (i.e. not set to 1 but set to 0)
# age_alive = age_alive[age_alive$dead %in% 0,]

# # Read in dead info
# age_dead = read.csv("/Cluster_Filespace/Marioni_Group/Danni/Proxies_stradl/KORA_train_recieved_020221/Cox/age_at_death.csv")

# # Assign the age in the age alive file to "aged" to match with aged at death file 
# names(age_alive)[7] <- "aged"

# # Bind the rows of the alive and dead files together 
# age = rbind(age_alive, age_dead)

# # Join up target files 
# tar1 <- readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/GS10k_Targets.rds")
# tar2 <- readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/wave4/w4-samplesheet.rds")

# tar1 <- tar1[c(1,2)]
# tar2$Sample_Sentrix_ID <- rownames(tar2)
# tar2 <- tar2[c("Sample_Name", "Sample_Sentrix_ID")]
# rownames(tar2) <- NULL
# tar2$Set <- "W4"
# target <- rbind(tar1, tar2)

# # Join target files to episcores
# names(episcores)[1] <- "Sample_Sentrix_ID"
# episcores <- left_join(episcores, target, by = "Sample_Sentrix_ID")

# # Join episcores information to main file with demographics
# names(age)[1] <- "Sample_Name"
# age$Sample_Name <- as.character(age$Sample_Name)
# d1 <- left_join(episcores, age, by = "Sample_Name")

# # Add age sex info in 
# agesex <- read.csv("/Cluster_Filespace/Marioni_Group/GS/GS_dataset/agesex_yobmob.csv")
# agesex <- agesex[c(1:3)]
# names(agesex)[1] <- "Sample_Name"
# agesex$Sample_Name <- as.character(agesex$Sample_Name)
# names(agesex)[2] <- "Age"
# d1 <- left_join(d1, agesex, by = "Sample_Name")

# # Isolate AD information 
# dem <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Dementia_GS/00_Data_exploration/Dementia_case_update_24Jan2022.csv")
# AD <- dem[c("id", "ad", "dt1_ad")]
# AD <- AD[which(AD$ad == "1"),] # 188
# AD <- AD[-2]
# AD <- AD[-which(duplicated(AD$id)),] # 136
# names(AD) <- c("id", "first")

# ######################################################################################################

# # Run the cox analyses in GS for AD, VD and joint dementia against episcores in 20k at >= 65 age at diagnosis

# library(survival)
# library(kinship2)
# library(coxme)
# library(readxl)
# library(tidyverse)
# library(gdata)

# # Read in the prepped file to cluster 
# ped <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Proxies_stradl/KORA_train_recieved_020221/Cox/pedigree_formatted.csv")
# # ped <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Dementia_GS/00_Data_exploration/2022-01-17_pedigree.csv")

# kin <- with(ped, pedigree(volid, father, mother, sex, famid=famid))
# kin_model <- kinship(kin) 

# # Function to Extract Lmekin Results	
# extract_coxme_table <- function (mod){
#   beta <- mod$coefficients #$fixed is not needed
#   nvar <- length(beta)
#   nfrail <- nrow(mod$var) - nvar
#   se <- sqrt(diag(mod$var)[nfrail + 1:nvar])
#   z<- round(beta/se, 2)
#   p<- signif(1 - pchisq((beta/se)^2, 1), 2)
#   table=data.frame(cbind(beta,se,z,p))
#   return(table)
# }

# # the usual clock for running all proteins
# clock <- names(d1)[2:110] 

# ## Run Incidence of Alzheimer's Disease Analysis - Restricted to Participants > 60 

# d1_AD <- d1

# mat_hazard_ad <- matrix(nrow=length(clock),ncol=9)
# output_hazard_AD<- as.data.frame(mat_hazard_ad)
# for(j in 1:length(clock)){ 
#   tryCatch({ 
#   dat1= d1_AD
#   tmp1 = AD[which(AD$id %in% dat1$Sample_Name),] # 110 individuals had AD and were over 60 at study baseline 
  
#   ## Obtain Age of Onset 
#   affected = dat1[which(dat1$Sample_Name %in% tmp1$id),] 
#   age_onset = AD[,c("first", "id")]
#   affected = merge(age_onset, affected, by.x = "id", by.y = "Sample_Name")
#   affected$Event = 1
#   affected$yoe = substring(affected$first, 1, 4)
#   affected$moe = substring(affected$first, 5,6)
#   affected$month_event1 = (as.numeric(affected$moe) - as.numeric(affected$mob))/12
#   affected$age_event1 = as.numeric(affected$yoe) - as.numeric(affected$yob)
#   affected$age_event = affected$age_event1 + affected$month_event1
#   affected$first = NULL
#   affected$yoe = NULL 
#   affected$moe = NULL
#   affected$month_event1 = NULL 
#   affected$age_event1 = NULL
  
#   healthy = dat1[-which(dat1$Sample_Name %in% AD$id),]
#   healthy$Event = 0
#   healthy$age_event = 0 
#   affected$id.y <- NULL
#   healthy$id <- NULL
#   names(affected)[names(affected)=="id"] <- "Sample_Name"
#   cox = rbind(affected, healthy)
  
#   cox$age_death = 0
#   cox$age_death = ifelse(cox$dead %in% 1, cox$aged, 0)
#   cox$age_at_event = ifelse(cox$Event %in% 1, cox$age_event, (ifelse(cox$dead %in% 1 & cox$Event %in% 0, cox$age_death, cox$aged)))
#   cox$tte = cox$age_at_event - cox$Age
#   cox$tte = as.numeric(cox$tte)
#   cox$tte <- ifelse(cox$tte < -1, "NA", cox$tte)
#   cox$tte = ifelse(cox$tte < 0, 0, cox$tte)
#   cox$Event = as.numeric(cox$Event)
#   cox$tte<-as.numeric(cox$tte)
  
#   cox = cox[cox$age_at_event >=65,]
#   mod = coxme(Surv(cox$tte, cox$Event) ~ scale(cox[,clock[[j]]]) + factor(cox$sex) + cox$Age + (1|cox$Sample_Name), varlist = kin_model*2)
#   print(clock[[j]])
#   print("AD")
#   output_hazard_AD[j,1] <- as.character(clock[[j]])
#   output_hazard_AD[j,2] <- as.character("Alzheimer's Disease")
#   output_hazard_AD[j,3:5]<-round(exp(cbind(coef(mod), confint(mod)))[1,1:3],2)
#   output_hazard_AD[j,6] <- extract_coxme_table(mod)[1,4]
#   output_hazard_AD[j,7] <- mod$n[1]
#   output_hazard_AD[j,8] <- mod$n[2]-mod$n[1]
#   # output_hazard_AD[j,9] <-cox.zph(mod)[1][[1]][3]
#   }, error = function(e) cat("skipped"))
# } 

# write.csv(output_hazard_AD, "/Cluster_Filespace/Marioni_Group/Danni/Dementia_GS/00_Data_exploration/AD_results_basic_no_set_070221.csv")


# ## Run a predictor of dementia subtypes using the episcores as input features in GS 20k 

