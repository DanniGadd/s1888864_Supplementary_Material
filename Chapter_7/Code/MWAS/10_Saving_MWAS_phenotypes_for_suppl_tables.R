### Save out phenotypes file from MWAS for full 778 with protein data 

## Add in GS ids instead of DNAm_ids 
stradl_ids <- read.table("/Cluster_Filespace/Marioni_Group/Rob/Somalogic/STRADL_DNAm_target_REM_17April2020.txt",header =T)
id.778 <- read.csv("/Cluster_Filespace/Marioni_Group/Rob/Somalogic/IDs_778.csv")

## Subset to the people in the ID file 
stradl_ids <- stradl_ids[which(stradl_ids$GS_id %in% id.778$x),]
ids = stradl_ids$DNAm_id


## Prepare phenotypes file for covariates that we will use to pre-correct data 
stradl_ids <- stradl_ids[,c("GS_id", "age_stradl", "sex", "wave", "Stradl_id")]

# Add batch 
w2_batch <- read.csv("/Cluster_Filespace/Marioni_Group/Rob/Somalogic/wave2_batch.csv") 
w2_batch$Batch[w2_batch$Batch == 1]<- "w1_1"
w2_batch$Batch[w2_batch$Batch == 2]<- "w1_2"
w2_batch$Batch[w2_batch$Batch == 3]<- "w1_3"
w2_batch$Batch[w2_batch$Batch == 4]<- "w1_4"
w2_batch$Batch[w2_batch$Batch == 5]<- "w1_5"
w2_batch$Batch[w2_batch$Batch == 6]<- "w1_6"
w3_batch <- read.csv("/Cluster_Filespace/Marioni_Group/Rob/Somalogic/wave_3_original_samplesheet.csv") 
w2_batch <- w2_batch[,c("Sample_Name","Batch")]
names(w2_batch) <- c("GS_id", "Batch")
w3_batch <- w3_batch[,c("Sample_Name","Batch_all")]
names(w3_batch) <- c("GS_id", "Batch")
w3_batch$GS_id <- gsub("ST", "", w3_batch$GS_id)
batch = rbind(w2_batch, w3_batch)
batch = batch[-which(duplicated(batch$GS_id)),]

phenotypes= merge(stradl_ids, batch, by = "GS_id")

# Join in SCID depression file
dep <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/Depression_case_controls_200521/prep_files/SCID_joint_to_GP_hospital_combined_SCID_only_270121.csv")
dep <- dep[c(1,2)]
names(dep)[2] <- "combined"
dep$combined <- as.character(dep$combined)
names(dep)[1] <- "Stradl_id"
# which(phenotypes$Stradl_id %in% dep$Stradl_id) 
phenotypes = left_join(phenotypes, dep, by = "Stradl_id")

# # Subset to 774 individuals with complete depression information
# sample <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/774_IDs_meth_order.csv")
# phenotypes <- phenotypes[which(phenotypes$Stradl_id %in% sample$Stradl_id),]

# Add WBC info into phenotypes
w2cells <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/w2_cells_joint.csv")
w3cells <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/w3_cells_joint.csv")
merged <- rbind(w2cells, w3cells)
names(merged)[1] <- "GS_id"

phenotypes = left_join(phenotypes, merged, by = "GS_id")

# Join in epismoker info for those in STRADL to linker file, then to cov file 
epi <- readRDS("/Cluster_Filespace/Marioni_Group/Danni/stradl_epismoker.rds")
names(epi)[1] <- "DNAm_id"
target <- read.table("/Cluster_Filespace/Marioni_Group/Rob/Somalogic/STRADL_DNAm_target_REM_17April2020.txt",header =T)
target <- target[c(2,3)]
names(target)[2] <- "DNAm_id"
library(tidyverse)
epi <- left_join(epi, target, by = "DNAm_id")

phenotypes = left_join(phenotypes, epi, by = "GS_id")

# Add BMI
BMI <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS/01_Phenotype_collation/prot_file_270121.csv")
BMI <- BMI[c("GS_id", "BMI")]

phenotypes = left_join(phenotypes, BMI, by = "GS_id")

# Write out phenotype file 
write.csv(phenotypes, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/demo_suppl_table/MWAS_phenotypes_file.csv", row.names = F)

