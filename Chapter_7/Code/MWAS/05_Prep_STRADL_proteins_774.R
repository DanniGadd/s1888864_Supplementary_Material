############################################################################################

### STRADL protein preparation

############################################################################################

cd /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Protein_preps/

screen

R

## read in the somalogic data ##
setwd("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Prep/")

soma <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/GS+soma+QC+normalized.csv")

## read in Archie's master linker file ##
link <- read.csv("ST id linkage.csv")

## read in the wave 2 target file ##
w2_tar <- read.csv("ALL-wave2.csv") 
a <- grep("S", w2_tar$Sample_Name)
w2_tar <- w2_tar[a,1:4]
w2_tar$id <- gsub("S", "", w2_tar$Sample_Name)

## merge target file with link file using GS id ##
w2_tar_update <- merge(link, w2_tar, by="id")

## filter to those who passed somalogic qc ##
b <- which(w2_tar_update$st %in% soma$SampleId)
w2_tar_soma <- w2_tar_update[b,]

## filter to stradl id (same as proteomics id), gs id, and DNAm id ##
w2_target <- w2_tar_soma[,c("st","id","Sample_Sentrix_ID", "sex", "age")]
names(w2_target) <- c("Stradl_id","GS_id","DNAm_id", "sex", "age_stradl")

## read in wave 3 target file ##
w3_tar <- read.csv("ST.csv")
w3_tar$id <- gsub("ST", "", w3_tar$Sample_Name)

## merge with link file ##
w3_tar_update <- merge(link, w3_tar, by="id")

## filter to those who passed somalogic qc ##
b <- which(w3_tar_update$st %in% soma$SampleId)
w3_tar_soma <- w3_tar_update[b,]

## filter to stradl id (same as proteomics id), gs id, and DNAm id ##
w3_target <- w3_tar_soma[,c("st","id","Sample_Sentrix_ID", "sex", "age")]
names(w3_target) <- c("Stradl_id","GS_id","DNAm_id", "sex", "age_stradl")

## harmonise into a single file ##
w3_target$wave <- "w3"
w2_target$wave <- "w2"

stradl_DNAm_target <- rbind(w2_target, w3_target)

# #write.table(stradl_DNAm_target, file="U:/Datastore/IGMM/marioni-lab/STRADL/STRADL_DNAm_target_REM_17April2020.txt", quote=F, col.names=T, row.names=F, sep="\t")
# # i will write a copy of this to my prep folder so i know how it was generated for the 844 individuals:
# write.table(stradl_DNAm_target, file="/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Prep/STRADL_DNAm_target_REM_17April2020_regen_140121.txt", quote=F, col.names=T, row.names=F, sep="\t")


# soma1 <- soma[which(soma$SampleId %in% stradl_DNAm_target$Stradl_id),]
# write.csv(soma1, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Prep/soma_844_file.csv", row.names = F)

####################################################################################

## Read in w2 and w3 cell count info and add into target

####################################################################################

## WAVE 2 
stradl_w2 <- stradl_DNAm_target[stradl_DNAm_target$wave %in% "w2",]
w2_cell <- read.csv("wave2_cell.csv")
w2_cell$Sample_Name <- gsub(".*S", "", w2_cell$Sample_Name)
w2_cell$Sample_Name <- gsub(".*R", "", w2_cell$Sample_Name)
w2_cell <- w2_cell[-which(duplicated(w2_cell$Sample_Name)),]
stradl_w2 <- merge(stradl_w2, w2_cell, by.x="GS_id", by.y="Sample_Name", all.x=T)
for(i in 7:12){ 
stradl_w2[,i][stradl_w2[,i] %in% NA] <- mean(stradl_w2[,i],na.rm = T)
} 

write.csv(stradl_w2, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/stradl_w2_cells.csv", row.names = F)


## WAVE 3 
stradl_w3 <- stradl_DNAm_target[stradl_DNAm_target$wave %in% "w3",]
w3_cell <- read.csv("w3_cell.csv")
w3_cell$Sample_Name <- gsub(".*ST", "", w3_cell$Sample_Name)
w3_cell <- w3_cell[-which(duplicated(w3_cell$Sample_Name)),]
stradl_w3 <- merge(stradl_w3, w3_cell, by.x = "GS_id", by.y="Sample_Name", all.x=T)

write.csv(stradl_w3, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/stradl_w3_cells.csv", row.names = F)


####################################################################################

## Read in w2 and w3 batch info and add into target

####################################################################################

## WAVE 2 

w2_batch <- read.csv("wave2_batch.csv")
stradl_w2 <- merge(stradl_w2, w2_batch, by.x="GS_id", by.y="Sample_Name", all.x=T)

## WAVE 3 

w3_batch <- read.csv("wave_3_original_samplesheet.csv")
w3_batch$Sample_Name <- gsub(".*ST", "", w3_batch$Sample_Name)
w3_batch = w3_batch[-which(duplicated(w3_batch$Sample_Name)),]
stradl_w3 <- merge(stradl_w3, w3_batch[,c("Sample_Name", "Batch_all")], by.x="GS_id",by.y="Sample_Name",all.x=T)
stradl_w3$Batch_all[stradl_w3$Batch_all == 1] <- 32
stradl_w3$Batch_all[stradl_w3$Batch_all == 3] <- 33
stradl_w3$Batch_all[stradl_w3$Batch_all == 4] <- 34
stradl_w3$Batch_all[stradl_w3$Batch_all == 5] <- 35
stradl_w3$Batch_all[stradl_w3$Batch_all == 6] <- 36
names(stradl_w3)[13] <- "Batch"
stradl <- rbind(stradl_w2, stradl_w3)

nrow(stradl) # this file has the WBCs and the batch for the 844 STRADL people

# Write out batch info 
batch <- stradl[c(1:3,13)]

write.csv(batch, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/stradl_batch.csv", row.names = F)

####################################################################################

## Add covariates into target file 

####################################################################################

soma_demo <- merge(link[,c("st","id", "age","sex")], soma, by.y = "SampleId",by.x="st") # merging in proteins with demographics 
demo <- read.csv("demographicsV2.csv")

soma_demo1 <- merge(soma_demo, demo[,c("st","st_age", "sample_dt")])
table(soma_demo1$age == soma_demo1$st_age)
soma_demo1$st_age <- NULL
library(stringi)
soma_demo1$study_site <- stri_sub(soma_demo1$st, 5,6)
soma_demo1$month_of_sample <- substring(soma_demo1$PlateRunDate, 4,5)
soma_demo1$month_of_draw <- substring(soma_demo1$sample_dt, 4,5)
soma_demo1$year_of_sample <- substring(soma_demo1$PlateRunDate, 7,10)
soma_demo1$year_of_draw <- substring(soma_demo1$sample_dt, 7,10)
soma_demo1$year_lag <- as.numeric(soma_demo1$year_of_sample) - as.numeric(soma_demo1$year_of_draw)
soma_demo1$month_lag <- as.numeric(soma_demo1$month_of_sample) - as.numeric(soma_demo1$month_of_draw)
soma_demo1$month_lag <- soma_demo1$month_lag/12
soma_demo1$lag_time <- soma_demo1$year_lag + soma_demo1$month_lag

eigen <- read.table("GS20K_ALL_MAF5_PCA.eigenvec")
head(eigen)
eigen <- eigen[which(eigen$V2 %in% soma_demo1$id),]
eigen$V1 <- NULL
names(eigen)[1] <- "id"
names(eigen)[2:21] <- paste0("PC", 1:20)
soma_demo1 <- merge(eigen, soma_demo1, by= "id")

# Assign lag group 
soma_demo1$lag_group <- 0
soma_demo1[soma_demo1$lag_time <= 2, "lag_group"] <- 1
soma_demo1[soma_demo1$lag_time > 2 & soma_demo1$lag_time <= 2.9, "lag_group"] <- 2
soma_demo1[soma_demo1$lag_time > 2.9 & soma_demo1$lag_time <= 3.5, "lag_group"] <- 3
soma_demo1[soma_demo1$lag_time > 3.5, "lag_group"] <- 4

# Save out copy for rosie 
write.csv(soma_demo1, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/soma_demo1_file_for_rosie.csv", row.names = F)

# Subset to just the 774 with methylation and depression status complete 
order <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/774_IDs_meth_order.csv")
soma_demo2 <- soma_demo1[soma_demo1$id %in% order$GS_id,]
ids = order$GS_id
soma_demo2 <- soma_demo2[match(ids,soma_demo2$id), ]
phenotypes <- soma_demo2

## Log Transform 
for(i in colnames(phenotypes)[56:4290]){ 
  phenotypes[,i]<- log(phenotypes[,i])
}

## Regress Proteins onto Covariates 
phenotypes_residualised <- phenotypes
for(i in colnames(phenotypes_residualised)[56:4290]){ 
  phenotypes_residualised[,i]<- lm(phenotypes_residualised[,i] ~ age + factor(sex) + factor(study_site) + factor(lag_group) + 
                                     PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20, 
                                    na.action = na.exclude, data = phenotypes_residualised)$residuals
}

## Rank-Inverse Based Normaliation on Pre-corrected Protein Phenotypes
library(bestNormalize)
for(i in colnames(phenotypes_residualised)[56:4290]){ 
  phenotypes_residualised[,i]<- orderNorm(phenotypes_residualised[,i])$x.t
}

## Scale somalogic data 
phenotypes_residualised[,56:4290] <- apply(phenotypes_residualised[,56:4290], 2, scale)

## Save out somalogic file whih has no covariates 
phe <- phenotypes_residualised[,c(1,22,56:4290)]	
write.csv(phe, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Protein_preps/Phenotype_file_774.csv", row.names = F)


