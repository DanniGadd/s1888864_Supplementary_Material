################################################################################################################################

### PREP ALL PheWAS PHENOTYPES IN STRADL

################################################################################################################################

cd /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/Depression_case_controls_200521/prep_files/

screen

R

library(tidyverse)
library(readxl)

# This script will format and combine all variables required for PheWAS analyses 
# These variables (N=1065 total) will be loaded into lmekin models for association tests
# Each variable will be tested for 4235 protein levels to identify protein markers of brain health outcomes

################################################################################################################################

### PROTEIN DATA

# Protein data (raw, without processing)
prot <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Data_for_inputs/GS+soma+QC+normalized.csv")

# Annotation linker file for SeqIds
link <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Data_for_inputs/annotation.csv", check.names = F)

# Update naming so were working in SeqIds 
names <- colnames(link)
names(prot)[c(33:4267)] <- names

# Load target file 
target <- read.delim("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Data_for_inputs/STRADL_DNAm_target_REM_17April2020.txt")

## TRANSFORM PROTEIN DATA AND JOIN FOR REGRESSIONS

## Log Transform 
for(i in colnames(prot)[33:4267]){ 
  prot[,i]<- log(prot[,i])
}


## Rank-Inverse Based Normaliation
library(bestNormalize)
for(i in colnames(prot)[33:4267]){ 
  prot[,i]<- orderNorm(prot[,i])$x.t
}

################################################################################################################################

### ADD PHENOTYPES

# Load demographics for all people in STRADL 
demo <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Data_for_inputs/demographicsV2.csv")
names(demo)[1] <- "SampleId"

# Join phenotype info to protein dataset 
prot <- left_join(prot, demo, by = "SampleId")

# Join in the GS id linkage so that further phenotypes can be joined in 
IDs <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Data_for_inputs/ST id linkage.csv")
names(IDs)[2] <- "GS_id"
names(IDs)[1] <- "SampleId"
IDs <- IDs[c(1,2)]
prot <- left_join(prot, IDs, by = "SampleId")

# Join in APOE
APOE <- read.csv("/Cluster_Filespace/Marioni_Group/GS/GS_dataset/APOE_snps/apoe_GS.csv")
names(APOE)[1] <- "GS_id"
prot <- left_join(prot, APOE, by = "GS_id")

# Check how many have each of the apoe phenotypes 
test <- prot
test$apoe %>% unique() #  "e3e4" "e3e3" "e2e3" NA     "e4e4" "e2e4" "e2e2"
table(is.na(test$apoe)) # 15 missing NA values 
outcomes <- test$apoe %>% as.data.frame() 
names(outcomes)[1] <- "X"
count(outcomes, X)

prot <- prot %>% mutate(APOE = case_when(
  apoe == "e4e4" ~ 2,
  apoe == "e3e4" ~ 2,
  apoe == "e3e3" ~ 1,
  apoe == "e2e2" ~ 0,
  apoe == "e2e3" ~ 0))

table(prot$APOE)
#   0   1   2
# 126 633 269

# Now join in the processed cognitive data from the script daniel shared with me - composite gf and g scores and other scores with outliers > 3.5 sd from mean removed 
comp <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS/01_Phenotype_collation/cog1_270121.csv")
names(comp)[1] <- "SampleId"
prot <- left_join(prot, comp, by = "SampleId")

# Add depression status into the dataset (read in file generated in depression covariate check, that is prepped with combined case/controls) 
dep <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/Depression_case_controls_200521/prep_files/SCID_joint_to_GP_hospital_combined_SCID_only_270121.csv")
dep <- dep[c(1,2)]
dep$DiagnosisGiven <- as.character(dep$DiagnosisGiven)
names(dep)[1] <- "SampleId"
prot <- left_join(prot, dep, by = "SampleId")
table(is.na(prot$DiagnosisGiven))
table(prot$DiagnosisGiven)

# Add BMI as a covariate 
demo_table <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/protAA/demo_table.csv")
BMI <- demo_table[c(14,54)]
prot = left_join(prot, BMI, by = "GS_id")
library(imputeTS)
prot$BMI = na_mean(prot$BMI) 


### ADD IN IMAGING DATA 

# Read in variables
brain_age <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Imaging_data/STRADL_Brain_Measures_BrainAge.csv")
age <- brain_age[c(1,5)]
names(age)[1] <- "ID"

vol <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Imaging_data/STRADL_Brain_Measures_Volumes.csv")

WM <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS/01_Phenotype_collation/Updated_WM/STRADL_Brain_Measures_gFA-gMD.csv")

WMHV <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_WMHV_110821/STRADL_WMHV_Measures_Complete_07.2021.csv")

faz <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Imaging_data/STRADL_Brain_Measures_Fazekas.csv")


### Join
vol <- full_join(vol, WM, by = "ID", all = TRUE)
vol <- full_join(vol, WMHV, by = "ID", all = TRUE)
vol <- full_join(vol, age, by = "ID", all = TRUE)
vol <- full_join(vol, faz, by = "ID", all = TRUE)
write.csv(vol, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS/01_Phenotype_collation/imaging_joint.csv", row.names = F)

# Recode points > 3.5 sd from mean for each variable (this has already been done for cognitive variables in the cognitive prep script for g and gf)
low = mean(vol$Global_GM_Volume, na.rm=T) - 3.5*sd(vol$Global_GM_Volume, na.rm=T)
high = mean(vol$Global_GM_Volume, na.rm=T) + 3.5*sd(vol$Global_GM_Volume, na.rm=T)
table(vol$Global_GM_Volume < low | vol$Global_GM_Volume > high)
vol$Global_GM_Volume[vol$Global_GM_Volume < low | vol$Global_GM_Volume > high] <- NA # no outliers GGM

low = mean(vol$WBV_No_Ventricles, na.rm=T) - 3.5*sd(vol$WBV_No_Ventricles, na.rm=T)
high = mean(vol$WBV_No_Ventricles, na.rm=T) + 3.5*sd(vol$WBV_No_Ventricles, na.rm=T)
table(vol$WBV_No_Ventricles < low | vol$WBV_No_Ventricles > high)
vol$WBV_No_Ventricles[vol$WBV_No_Ventricles < low | vol$WBV_No_Ventricles > high] <- NA # one outlier WBV

low = mean(vol$gFA, na.rm=T) - 3.5*sd(vol$gFA, na.rm=T)
high = mean(vol$gFA, na.rm=T) + 3.5*sd(vol$gFA, na.rm=T)
table(vol$gFA < low | vol$gFA > high)
vol$gFA[vol$gFA < low | vol$gFA > high] <- NA # one outlier gFA

low = mean(vol$gMD, na.rm=T) - 3.5*sd(vol$gMD, na.rm=T)
high = mean(vol$gMD, na.rm=T) + 3.5*sd(vol$gMD, na.rm=T)
table(vol$gMD < low | vol$gMD > high)
vol$gMD[vol$gMD < low | vol$gMD > high] <- NA # one outlier gMD

low = mean(vol$Brain_age, na.rm=T) - 3.5*sd(vol$Brain_age, na.rm=T)
high = mean(vol$Brain_age, na.rm=T) + 3.5*sd(vol$Brain_age, na.rm=T)
table(vol$Brain_age < low | vol$Brain_age > high)
vol$Brain_age[vol$Brain_age < low | vol$Brain_age > high] <- NA # no outliers brain age 

low = mean(vol$Fazekas_Score_Total, na.rm=T) - 3.5*sd(vol$Fazekas_Score_Total, na.rm=T)
high = mean(vol$Fazekas_Score_Total, na.rm=T) + 3.5*sd(vol$Fazekas_Score_Total, na.rm=T)
table(vol$Fazekas_Score_Total < low | vol$Fazekas_Score_Total > high)
vol$Fazekas_Score_Total[vol$Fazekas_Score_Total < low | vol$Fazekas_Score_Total > high] <- NA # 4 outliers

# +1 and log transform the WMHV variable
vol$WMH_Volume_Total <- as.numeric(vol$WMH_Volume_Total)
vol$WMH_Volume_Total <- vol$WMH_Volume_Total + 1
vol$WMH_Volume_Total <- log(vol$WMH_Volume_Total)

low = mean(vol$WMH_Volume_Total, na.rm=T) - 3.5*sd(vol$WMH_Volume_Total, na.rm=T)
high = mean(vol$WMH_Volume_Total, na.rm=T) + 3.5*sd(vol$WMH_Volume_Total, na.rm=T)
table(vol$WMH_Volume_Total < low | vol$WMH_Volume_Total > high)
# vol$gMD[vol$WMH_Volume_Total < low | vol$WMH_Volume_Total > high] <- NA # 6 outlying values for WMHV
vol$WMH_Volume_Total[vol$WMH_Volume_Total < low | vol$WMH_Volume_Total > high] <- NA # 6 outlying values for WMHV

### Join the imaging data into the protein dataset for regressions 
# We want imaging data for all of the individuals with protein data, so will leftjoin 
names(vol)[1] <- "SampleId"
prot <- left_join(prot, vol, by = "SampleId")

# Calculate brain age acceleration score 
prot$brain_accel <- resid(lm(Brain_age ~ st_age, na.action = na.exclude, data = prot))

# add in lag group and study site 
ewas <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/soma_demo1_file_for_rosie.csv")
var <- ewas[c("st","study_site", "lag_group")]
names(var)[1] <- "SampleId"

# Add them into the main protein file 
prot <- left_join(prot, var, by = "SampleId")

# Save out prot file that is prepped wth proteins and phenotypes
write.csv(prot, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS/01_Phenotype_collation/prot_file_220221.csv", row.names = F)

###################

