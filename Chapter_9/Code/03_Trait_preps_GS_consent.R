######################################################################################################

### TRAITS IN GS - DATA LINKAGE (CONSENT GP, PLUS FULL ICD)

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
dem <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/Dementia_case_update_18Jan2023.csv")

# Consent
consent <- read_excel("/Cluster_Filespace/Marioni_Group/Danni/Dementia_GS/00_Data_exploration/2020-09-02 GP consent ids (1).xlsx")
consent <- as.data.frame(consent)
consent$Consent <- "YES"
consent <- consent[c(1,5)]
dem <- left_join(dem, consent, by = "id")
dem <- dem[which(dem$dementia == "1"),]
gp <- dem[which(dem$source == "gp"),]
table(gp$Consent) # 89 of 249 entries
gp <- gp[which(gp$Consent %in% 'YES'),]
smr <- dem[which(dem$source %in% 'smr01'),]

all <- rbind(gp, smr)
# all <- all[-which(all$dlb == '1' | all$ftd == '1'),]

names(all)[2] <- c("id")
names(all)[11] <- c("first")
all <- all[order(all$first),] # Check that the file is ordered by date (lowest to highest)
all <- all[which(all$first > 199001),] # Remove codes pre 199001
all2 <- all[-which(duplicated(all$id)),] # Remove any duplicate ID codes by taking the earliest date for that person (if ordered this should take first instance)
dim(all2) # 290

# Summarise subgroup counts after exclusions
table(all2$dementia) # 290
table(all2$ad) # 99
table(all2$vd) # 68
table(all2$dlb) # 1
table(all2$ftd) # 3

all2 <- all2[c(2,11)]
write.csv(all2, "/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_All_initial_files/all_dementia_consent.csv", row.names = F)




## Stroke 
Stroke <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/00_Rerun/Cox/Input_tables/Run_by_Danni/Disease_files/Stroke.csv")
Stroke_2 <- read_excel("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/00_Rerun/Cox/Input_tables/Run_by_Danni/Disease_files_GP/Stroke.xlsx", col_types = c("text", "text", "text", "date", "text"))

Stroke_2 <- Stroke_2[which(Stroke_2$id %in% consent$id),]

Stroke_2 <- Stroke_2[-grep("Injury|injury", Stroke_2$description),]
Stroke_2 <- Stroke_2[-grep("Mitochondrial", Stroke_2$description),]
# Stroke_2 <- Stroke_2[-grep("Personal", Stroke_2$description),]
Stroke_2 <- Stroke_2[,c(5,4)]
names(Stroke_2) <- c("id", "first")
Stroke_2$first <- paste0(substring(Stroke_2$first, 1, 4), substring(Stroke_2$first, 6,7))
Stroke <- Stroke[,c(1,2)]
Stroke = rbind(Stroke, Stroke_2)
Stroke <- Stroke[order(Stroke$first),]
Stroke <- Stroke[which(Stroke$first > 199001),]
Stroke = Stroke[-which(duplicated(Stroke$id)),] # 986 drop to 696
write.csv(Stroke, "/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_All_initial_files/Stroke_combined_consent.csv", row.names = F)


## Heart Disease 
IHD <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/00_Rerun/Cox/Input_tables/Run_by_Danni/Disease_files/Heart.csv")
IHD_2 <- read_excel("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/00_Rerun/Cox/Input_tables/Run_by_Danni/Disease_files_GP/IHD.xlsx", col_types = c("text", "text", "text", "date", "text"))

IHD_2 <- IHD_2[which(IHD_2$id %in% consent$id),]

IHD_2 <- IHD_2[,c(5,4)]
names(IHD_2) <- c("id", "first")
IHD_2$first <- paste0(substring(IHD_2$first, 1, 4), substring(IHD_2$first, 6,7))
IHD <- IHD[,c(1,2)]
IHD = rbind(IHD, IHD_2)
IHD <- IHD[order(IHD$first),]
IHD <- IHD[which(IHD$first > 199001),]
IHD = IHD[-which(duplicated(IHD$id)),] # 1715 to 1549
write.csv(IHD, "/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_All_initial_files/IHD_combined_consent.csv", row.names = F)


# Diabetes
Diabetes <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/GS_COVID/Daniel_diab_comparison/2022-15-03_diabetes_smr_update.csv")
Diabetes <- Diabetes[-grep("MODY|1|Other|Unknown", Diabetes$tname),] # remove type 1 and other cases
Diabetes <- na.omit(Diabetes)
names(Diabetes)[4] <- "first"
Diabetes_2 <- read_excel("/Cluster_Filespace/Marioni_Group/Danni/GS_COVID/Diabetes_update_codes_yipeng/Diabetes.xlsx", col_types = c("text", "text", "text", "date", "text")) 

Diabetes_2 <- Diabetes_2[which(Diabetes_2$id %in% consent$id),]


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
Diabetes <- Diabetes[-which(duplicated(Diabetes$id)),] # 1652 to 1355
Diabetes <- Diabetes[c(1,4)]
write.csv(Diabetes, "/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_All_initial_files/Diabetes_combined_consent.csv", row.names = F)


