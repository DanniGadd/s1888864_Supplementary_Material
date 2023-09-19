################################################################################################################################

### PREP DEPRESSION CASES AND CONTROLS STRADL  

################################################################################################################################

cd /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/Depression_case_controls_200521/prep_files/

screen

R

library(tidyverse)
library(readxl)

################################################################################################################################

### PROTEIN FILE - LOAD 

# Protein data (raw, without processing)
prot <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Data_for_inputs/GS+soma+QC+normalized.csv")

# Annotation linker file for SeqIds
link <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Data_for_inputs/annotation.csv", check.names = F)

# Update naming so were working in SeqIds 
names <- colnames(link)
names(prot)[c(33:4267)] <- names

# Load target file 
target <- read.delim("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Data_for_inputs/STRADL_DNAm_target_REM_17April2020.txt")

# Check for missing protein data 
proteins <- prot[c(13,33:4267)]
table(is.na(proteins))

## Log Transform 
for(i in colnames(prot)[33:4267]){ 
  prot[,i]<- log(prot[,i])
}

## Rank-Inverse Based Normaliation
library(bestNormalize)
for(i in colnames(prot)[33:4267]){ 
  prot[,i]<- orderNorm(prot[,i])$x.t
}

# Load demographics for all people in STRADL 
demo <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Data_for_inputs/demographicsV2.csv")
names(demo)[1] <- "SampleId"

# > dim(demo)
# [1] 1069   11

# Join phenotype info to protein dataset 
prot <- left_join(prot, demo, by = "SampleId")

# Scale age data prior to running model 
# prot[,4269:4270] <- apply(prot[,4269:4270], 2, scale)

# Join in the GS id linkage so that further phenotypes can be joined in 
IDs <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Data_for_inputs/ST id linkage.csv")
names(IDs)[2] <- "GS_id"
names(IDs)[1] <- "SampleId"
IDs <- IDs[c(1,2)]

prot <- left_join(prot, IDs, by = "SampleId")

# Join in APOE and cognitive phenotypes 
APOE <- read.csv("/Cluster_Filespace/Marioni_Group/GS/GS_dataset/APOE_snps/apoe_GS.csv")
names(APOE)[1] <- "GS_id"

prot <- left_join(prot, APOE, by = "GS_id")

# Cognitive scores - read in and join to dataset 
cog <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Data_for_inputs/st_cognitive.csv")
names(cog)[1] <- "SampleId"

prot <- left_join(prot, cog, by = "SampleId")

# Now join in the processed cognitive data from the code daniel shared with me - composite gf and g scores 
comp <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Data_for_inputs/cog1_240321.csv")
names(comp)[1] <- "SampleId"
comp <- comp[c(1,5,6,7)]

prot <- left_join(prot, comp, by = "SampleId")

################################################################################################################################

### DEPRESSION FILE - LOAD AND CALCULATE CASES/CONTROLS

# Read in the depression SCID data diagnoses from Miruna 
SCID <- read_excel("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/Depression_case_controls_200521/SCID.xlsx")
SCID <- as.data.frame(SCID)

diagnosed <- SCID[which(SCID$DiagnosisGiven %in% "1"),] # 366
dim(diagnosed)

# Read in depression codes from GP hospital extraction Cliff 
codes <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/Depression_case_controls_200521/prep_files/Depression_mod_severe_combined.csv")

# Subset to those with consent 
consent <- read_excel("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/Depression_case_controls_200521/prep_files/2020-09-02 GP consent ids (1).xlsx")

codes <- codes[which(codes$id %in% consent$id),]

# I need to work out the individuals that had recieved MDD diagnoses leading up to or at STRADL sampling
affected = prot[which(prot$GS_id %in% codes$id),] # 51 people from STRADL 1,065 with depression
age_onset = codes[,c("first", "id")]
affected = merge(age_onset, affected, by.x = "id", by.y = "GS_id")
affected$Event = 1
affected$yoe = substring(affected$first, 1, 4)
affected$moe = substring(affected$first, 5,6)
mob <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/Depression_case_controls_200521/prep_files/agesex_yobmob.csv")
affected <- left_join(affected, mob, by = "id")
affected$month_event1 = (as.numeric(affected$moe) - as.numeric(affected$mob))/12
affected$age_event1 = as.numeric(affected$yoe) - as.numeric(affected$yob.y)
affected$age_event = affected$age_event1 + affected$month_event1
affected$diff = affected$st_age - affected$age_event
affected$linkage <- ifelse(affected$diff > 0, "1", "0")
cases <- affected[which(affected$linkage == "1"),]
dim(cases) 
case_ID <- cases[,which(colnames(cases) %in% c("id", "SampleId", "linkage"))]
names(case_ID) <- c("GS", "ID", "linkage")

# Combine SCID and MDD codes for depression 
target <- read.delim("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Data_for_inputs/STRADL_DNAm_target_REM_17April2020.txt")
IDs <- target[c(1:2)]
SCID2 <- left_join(SCID, case_ID, by ="ID")
SCID2$combined <- ifelse(SCID2$DiagnosisGiven == 1 | SCID2$linkage == "1", 1, 0)
SCID2[is.na(SCID2)] <- "0"
names(SCID2)[1] <- "SampleId"

# Save a copy such that combined depression can be used as a covariate in the EWAS too, and rest of the pheWAS 
write.csv(SCID2, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/Depression_case_controls_200521/prep_files/SCID_joint_to_GP_hospital_combined_SCID_only_270121.csv", row.names = F)

# Overall, there are no additional cases derived from the lookup of incident cases 
# We can therefore use the SCID information as is 