################################################################################################################################

### PheWAS analyses - Age and Sex

################################################################################################################################


screen

R

library(tidyverse)
library(readxl)

## Load prepped joint PheWAS phenotypes file 
prot <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS/01_Phenotype_collation/prot_file_220221.csv", check.names = F)

## Load lmekin requirements
library(survival)
library(kinship2)
library(coxme)
library(readxl)
library(tidyverse)
library(gdata)

## Code that is already processed and read in as the ped file below
# ped = read.csv("/Cluster_Filespace/Marioni_Group/Danni/Proxies_stradl/KORA_train_recieved_020221/Cox/pedigree.csv")
# ped$father <- as.numeric(ped$father)
# ped$mother <- as.numeric(ped$mother)
# ped$father[ped$father==0] <- NA
# ped$mother[ped$mother==0] <- NA
# table(ped$sex)
# ped$sex <- as.numeric(ped$sex)
# ped$sex[ped$sex==2] <- 0
# ped$sex <- ped$sex+1

# Read in the prepped file to cluster 
ped <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Proxies_stradl/KORA_train_recieved_020221/Cox/pedigree_formatted.csv")

# Create kinship matrix for GS
kin <- with(ped, pedigree(volid, father, mother, sex, famid=famid)) # Pedigree list with 26 total subjects in 5 families
kin_model <- kinship(kin) 

# Function to Extract Lmekin Results  
extract_coxme_table <- function (mod){
  #beta <- mod$coefficients #$fixed is not needed
  beta <- fixef(mod)
  nvar <- length(beta)
  nfrail <- nrow(mod$var) - nvar
  se <- sqrt(diag(mod$var)[nfrail + 1:nvar])
  z<- beta/se
  p<- 1 - pchisq((beta/se)^2, 1)
  table=data.frame(cbind(beta,se,z,p))
  return(table)
}

## Set marker names to loop through in models as proteins (x)
markers <- prot[c(33:4267)]
marker_names <- colnames(markers)

# Rename depression variable
names(prot)[4297] <- "combined"



#####################################################################################################################################

### AGE AND SEX PHEWAS

#####################################################################################################################################

### BATCHES

length <- 4235
results <- data.frame(SeqId = 1:length, Age_beta = 1:length, Age_SE = 1:length, Age_P = 1:length)
results2 <- data.frame(SeqId = 1:length, Sex_beta = 1:length, Sex_SE = 1:length, Sex_P = 1:length)

for(i in 1:500){
  prot_name <- as.character(colnames(markers[i]))
  
  mod <- lmekin(scale(prot[,prot_name]) ~ scale(prot$st_age) + as.factor(prot$sex) + as.factor(prot$combined) + as.factor(study_site) + as.factor(lag_group) + (1|prot$GS_id), data = prot, varlist = kin_model*2)

  print(i)
  print(prot_name)
  results[i,1] <- prot_name
  results[i,2] <- extract_coxme_table(mod)[2,1]
  results[i,3] <- extract_coxme_table(mod)[2,2]
  results[i,4] <- extract_coxme_table(mod)[2,4]

  results2[i,1] <- prot_name
  results2[i,2] <- extract_coxme_table(mod)[3,1]
  results2[i,3] <- extract_coxme_table(mod)[3,2]
  results2[i,4] <- extract_coxme_table(mod)[3,4]
}

# save a copy of results tables 
write.csv(results, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/Age_sex/batches/age_b1.csv", row.names = F)
write.csv(results2, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/Age_sex/batches/sex_b1.csv", row.names = F)


length <- 4235
results <- data.frame(SeqId = 1:length, Age_beta = 1:length, Age_SE = 1:length, Age_P = 1:length)
results2 <- data.frame(SeqId = 1:length, Sex_beta = 1:length, Sex_SE = 1:length, Sex_P = 1:length)

for(i in 501:1000){
  prot_name <- as.character(colnames(markers[i]))
  
  mod <- lmekin(scale(prot[,prot_name]) ~ scale(prot$st_age) + as.factor(prot$sex) + as.factor(prot$combined) + as.factor(study_site) + as.factor(lag_group) + (1|prot$GS_id), data = prot, varlist = kin_model*2)

  print(i)
  print(prot_name)
  results[i,1] <- prot_name
  results[i,2] <- extract_coxme_table(mod)[2,1]
  results[i,3] <- extract_coxme_table(mod)[2,2]
  results[i,4] <- extract_coxme_table(mod)[2,4]

  results2[i,1] <- prot_name
  results2[i,2] <- extract_coxme_table(mod)[3,1]
  results2[i,3] <- extract_coxme_table(mod)[3,2]
  results2[i,4] <- extract_coxme_table(mod)[3,4]
}

# save a copy of results tables 
write.csv(results, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/Age_sex/batches/age_b2.csv", row.names = F)
write.csv(results2, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/Age_sex/batches/sex_b2.csv", row.names = F)


length <- 4235
results <- data.frame(SeqId = 1:length, Age_beta = 1:length, Age_SE = 1:length, Age_P = 1:length)
results2 <- data.frame(SeqId = 1:length, Sex_beta = 1:length, Sex_SE = 1:length, Sex_P = 1:length)

for(i in 1001:1500){
  prot_name <- as.character(colnames(markers[i]))
  
  mod <- lmekin(scale(prot[,prot_name]) ~ scale(prot$st_age) + as.factor(prot$sex) + as.factor(prot$combined) + as.factor(study_site) + as.factor(lag_group) + (1|prot$GS_id), data = prot, varlist = kin_model*2)

  print(i)
  print(prot_name)
  results[i,1] <- prot_name
  results[i,2] <- extract_coxme_table(mod)[2,1]
  results[i,3] <- extract_coxme_table(mod)[2,2]
  results[i,4] <- extract_coxme_table(mod)[2,4]

  results2[i,1] <- prot_name
  results2[i,2] <- extract_coxme_table(mod)[3,1]
  results2[i,3] <- extract_coxme_table(mod)[3,2]
  results2[i,4] <- extract_coxme_table(mod)[3,4]
}

# save a copy of results tables 
write.csv(results, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/Age_sex/batches/age_b3.csv", row.names = F)
write.csv(results2, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/Age_sex/batches/sex_b3.csv", row.names = F)


length <- 4235
results <- data.frame(SeqId = 1:length, Age_beta = 1:length, Age_SE = 1:length, Age_P = 1:length)
results2 <- data.frame(SeqId = 1:length, Sex_beta = 1:length, Sex_SE = 1:length, Sex_P = 1:length)

for(i in 1501:2000){
  prot_name <- as.character(colnames(markers[i]))
  
  mod <- lmekin(scale(prot[,prot_name]) ~ scale(prot$st_age) + as.factor(prot$sex) + as.factor(prot$combined) + as.factor(study_site) + as.factor(lag_group) + (1|prot$GS_id), data = prot, varlist = kin_model*2)

  print(i)
  print(prot_name)
  results[i,1] <- prot_name
  results[i,2] <- extract_coxme_table(mod)[2,1]
  results[i,3] <- extract_coxme_table(mod)[2,2]
  results[i,4] <- extract_coxme_table(mod)[2,4]

  results2[i,1] <- prot_name
  results2[i,2] <- extract_coxme_table(mod)[3,1]
  results2[i,3] <- extract_coxme_table(mod)[3,2]
  results2[i,4] <- extract_coxme_table(mod)[3,4]
}

# batch 4
# save a copy of results tables 
write.csv(results, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/Age_sex/batches/age_b4.csv", row.names = F)
write.csv(results2, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/Age_sex/batches/sex_b4.csv", row.names = F)



length <- 4235
results <- data.frame(SeqId = 1:length, Age_beta = 1:length, Age_SE = 1:length, Age_P = 1:length)
results2 <- data.frame(SeqId = 1:length, Sex_beta = 1:length, Sex_SE = 1:length, Sex_P = 1:length)

for(i in 2001:2500){
  prot_name <- as.character(colnames(markers[i]))
  
  mod <- lmekin(scale(prot[,prot_name]) ~ scale(prot$st_age) + as.factor(prot$sex) + as.factor(prot$combined) + as.factor(study_site) + as.factor(lag_group) + (1|prot$GS_id), data = prot, varlist = kin_model*2)

  print(i)
  print(prot_name)
  results[i,1] <- prot_name
  results[i,2] <- extract_coxme_table(mod)[2,1]
  results[i,3] <- extract_coxme_table(mod)[2,2]
  results[i,4] <- extract_coxme_table(mod)[2,4]

  results2[i,1] <- prot_name
  results2[i,2] <- extract_coxme_table(mod)[3,1]
  results2[i,3] <- extract_coxme_table(mod)[3,2]
  results2[i,4] <- extract_coxme_table(mod)[3,4]
}

# save a copy of results tables 
write.csv(results, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/Age_sex/batches/age_b5.csv", row.names = F)
write.csv(results2, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/Age_sex/batches/sex_b5.csv", row.names = F)


length <- 4235
results <- data.frame(SeqId = 1:length, Age_beta = 1:length, Age_SE = 1:length, Age_P = 1:length)
results2 <- data.frame(SeqId = 1:length, Sex_beta = 1:length, Sex_SE = 1:length, Sex_P = 1:length)

for(i in 2501:3000){
  prot_name <- as.character(colnames(markers[i]))
  
  mod <- lmekin(scale(prot[,prot_name]) ~ scale(prot$st_age) + as.factor(prot$sex) + as.factor(prot$combined) + as.factor(study_site) + as.factor(lag_group) + (1|prot$GS_id), data = prot, varlist = kin_model*2)

  print(i)
  print(prot_name)
  results[i,1] <- prot_name
  results[i,2] <- extract_coxme_table(mod)[2,1]
  results[i,3] <- extract_coxme_table(mod)[2,2]
  results[i,4] <- extract_coxme_table(mod)[2,4]

  results2[i,1] <- prot_name
  results2[i,2] <- extract_coxme_table(mod)[3,1]
  results2[i,3] <- extract_coxme_table(mod)[3,2]
  results2[i,4] <- extract_coxme_table(mod)[3,4]
}

# save a copy of results tables 
write.csv(results, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/Age_sex/batches/age_b6.csv", row.names = F)
write.csv(results2, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/Age_sex/batches/sex_b6.csv", row.names = F)



length <- 4235
results <- data.frame(SeqId = 1:length, Age_beta = 1:length, Age_SE = 1:length, Age_P = 1:length)
results2 <- data.frame(SeqId = 1:length, Sex_beta = 1:length, Sex_SE = 1:length, Sex_P = 1:length)

for(i in 3001:3500){
  prot_name <- as.character(colnames(markers[i]))
  
  mod <- lmekin(scale(prot[,prot_name]) ~ scale(prot$st_age) + as.factor(prot$sex) + as.factor(prot$combined) + as.factor(study_site) + as.factor(lag_group) + (1|prot$GS_id), data = prot, varlist = kin_model*2)

  print(i)
  print(prot_name)
  results[i,1] <- prot_name
  results[i,2] <- extract_coxme_table(mod)[2,1]
  results[i,3] <- extract_coxme_table(mod)[2,2]
  results[i,4] <- extract_coxme_table(mod)[2,4]

  results2[i,1] <- prot_name
  results2[i,2] <- extract_coxme_table(mod)[3,1]
  results2[i,3] <- extract_coxme_table(mod)[3,2]
  results2[i,4] <- extract_coxme_table(mod)[3,4]
}

# save a copy of results tables 
write.csv(results, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/Age_sex/batches/age_b7.csv", row.names = F)
write.csv(results2, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/Age_sex/batches/sex_b7.csv", row.names = F)



length <- 4235
results <- data.frame(SeqId = 1:length, Age_beta = 1:length, Age_SE = 1:length, Age_P = 1:length)
results2 <- data.frame(SeqId = 1:length, Sex_beta = 1:length, Sex_SE = 1:length, Sex_P = 1:length)

for(i in 3501:4000){
  prot_name <- as.character(colnames(markers[i]))
  
  mod <- lmekin(scale(prot[,prot_name]) ~ scale(prot$st_age) + as.factor(prot$sex) + as.factor(prot$combined) + as.factor(study_site) + as.factor(lag_group) + (1|prot$GS_id), data = prot, varlist = kin_model*2)

  print(i)
  print(prot_name)
  results[i,1] <- prot_name
  results[i,2] <- extract_coxme_table(mod)[2,1]
  results[i,3] <- extract_coxme_table(mod)[2,2]
  results[i,4] <- extract_coxme_table(mod)[2,4]

  results2[i,1] <- prot_name
  results2[i,2] <- extract_coxme_table(mod)[3,1]
  results2[i,3] <- extract_coxme_table(mod)[3,2]
  results2[i,4] <- extract_coxme_table(mod)[3,4]
}

# save a copy of results tables 
write.csv(results, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/Age_sex/batches/age_b8.csv", row.names = F)
write.csv(results2, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/Age_sex/batches/sex_b8.csv", row.names = F)


length <- 4235
results <- data.frame(SeqId = 1:length, Age_beta = 1:length, Age_SE = 1:length, Age_P = 1:length)
results2 <- data.frame(SeqId = 1:length, Sex_beta = 1:length, Sex_SE = 1:length, Sex_P = 1:length)

for(i in 4001:4235){
  prot_name <- as.character(colnames(markers[i]))
  
  mod <- lmekin(scale(prot[,prot_name]) ~ scale(prot$st_age) + as.factor(prot$sex) + as.factor(prot$combined) + as.factor(study_site) + as.factor(lag_group) + (1|prot$GS_id), data = prot, varlist = kin_model*2)

  print(i)
  print(prot_name)
  results[i,1] <- prot_name
  results[i,2] <- extract_coxme_table(mod)[2,1]
  results[i,3] <- extract_coxme_table(mod)[2,2]
  results[i,4] <- extract_coxme_table(mod)[2,4]

  results2[i,1] <- prot_name
  results2[i,2] <- extract_coxme_table(mod)[3,1]
  results2[i,3] <- extract_coxme_table(mod)[3,2]
  results2[i,4] <- extract_coxme_table(mod)[3,4]
}

# save a copy of results tables 
write.csv(results, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/Age_sex/batches/age_b9.csv", row.names = F)
write.csv(results2, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/Age_sex/batches/sex_b9.csv", row.names = F)


#########################################################################################################################

### PROCESSING AGE SEX RESULTS

#########################################################################################################################

screen

R

library(tidyverse)

# Sort the results for each model and combine the batches

# Age 

e1 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/Age_sex/batches/age_b1.csv")
e2 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/Age_sex/batches/age_b2.csv")
e3 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/Age_sex/batches/age_b3.csv")
e4 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/Age_sex/batches/age_b4.csv")
e5 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/Age_sex/batches/age_b5.csv")
e6 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/Age_sex/batches/age_b6.csv")
e7 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/Age_sex/batches/age_b7.csv")
e8 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/Age_sex/batches/age_b8.csv")
e9 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/Age_sex/batches/age_b9.csv")

e1 <- e1[1:500,]
e2 <- e2[501:1000,]
e3 <- e3[1001:1500,]
e4 <- e4[1501:2000,]
e5 <- e5[2001:2500,]
e6 <- e6[2501:3000,]
e7 <- e7[3001:3500,]
e8 <- e8[3501:4000,]
e9 <- e9[4001:4235,]

e <- rbind(e1, e2)
e <- rbind(e, e3)
e <- rbind(e, e4)
e <- rbind(e, e5)
e <- rbind(e, e6)
e <- rbind(e, e7)
e <- rbind(e, e8)
e <- rbind(e, e9)

results <- e

# Sex 

e1 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/Age_sex/batches/sex_b1.csv")
e2 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/Age_sex/batches/sex_b2.csv")
e3 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/Age_sex/batches/sex_b3.csv")
e4 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/Age_sex/batches/sex_b4.csv")
e5 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/Age_sex/batches/sex_b5.csv")
e6 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/Age_sex/batches/sex_b6.csv")
e7 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/Age_sex/batches/sex_b7.csv")
e8 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/Age_sex/batches/sex_b8.csv")
e9 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/Age_sex/batches/sex_b9.csv")

e1 <- e1[1:500,]
e2 <- e2[501:1000,]
e3 <- e3[1001:1500,]
e4 <- e4[1501:2000,]
e5 <- e5[2001:2500,]
e6 <- e6[2501:3000,]
e7 <- e7[3001:3500,]
e8 <- e8[3501:4000,]
e9 <- e9[4001:4235,]

e <- rbind(e1, e2)
e <- rbind(e, e3)
e <- rbind(e, e4)
e <- rbind(e, e5)
e <- rbind(e, e6)
e <- rbind(e, e7)
e <- rbind(e, e8)
e <- rbind(e, e9)

results2 <- e

# Rank
results <- results[order(results$Age_P),]
results2 <- results2[order(results2$Sex_P),]

# Merge in protein info 
anno <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Annotations_protein_file/annotation_formatted_for_paper.csv")

results <- left_join(results, anno, by = "SeqId")

results2 <- left_join(results2, anno, by = "SeqId")

# record of previous adjustment:

# # Manually calculate p vals 
# results$Age_Pcalc <- pchisq((results$Age_beta/results$Age_SE)^2, 1, lower.tail=F)
# results2$Sex_Pcalc <- pchisq((results2$Sex_beta/results2$Sex_SE)^2, 1, lower.tail=F)

# # # Plot p vals extracted vs calculated 
# # plot(results$Age_P, results$Age_Pcalc) # these look good

# # Do p value adjustment 
# results$Age_P_adjust <- p.adjust(results$Age_Pcalc, method = "BH")
# results2$Sex_P_adjust <- p.adjust(results2$Sex_Pcalc, method = "BH")

# # Add identifier for those that pass/fail significance in each 
# results$Age_Status <- ifelse(results$Age_P_adjust < 0.05, "pass", "fail")
# results2$Sex_Status <- ifelse(results2$Sex_P_adjust < 0.05, "pass", "fail")

# # Count those that pass
# sig <- which(results$Age_P_adjust < 0.05)
# sig <- results[sig,]
# dim(sig) # 800 - to 919 after study site and lag group adjusted for 

# # Count those that pass
# sig2 <- which(results2$Sex_P_adjust < 0.05)
# sig2 <- results2[sig2,]
# dim(sig2) # 805 - to 839 after study site and lag group adjusted for 


# Manually calculate p vals 
results$Age_Pcalc <- pchisq((results$Age_beta/results$Age_SE)^2, 1, lower.tail=F)
results2$Sex_Pcalc <- pchisq((results2$Sex_beta/results2$Sex_SE)^2, 1, lower.tail=F)

# # Plot p vals extracted vs calculated 
# plot(results$Age_P, results$Age_Pcalc) # these look good

# # Do p value adjustment 
# results$Age_P_adjust <- p.adjust(results$Age_Pcalc, method = "BH")
# results2$Sex_P_adjust <- p.adjust(results2$Sex_Pcalc, method = "BH")

# Threshold significance by P < 3.5x10-4

# Add identifier for those that pass/fail significance in each 
results$Age_Status <- ifelse(results$Age_Pcalc < 3.5e-4, "pass", "fail")
results2$Sex_Status <- ifelse(results2$Sex_Pcalc < 3.5e-4, "pass", "fail")

# Count those that pass
sig <- which(results$Age_Pcalc < 0.05)
sig <- results[sig,]
dim(sig) # 1213

# Count those that pass
sig2 <- which(results2$Sex_Pcalc < 0.05)
sig2 <- results2[sig2,]
dim(sig2) # 1210

# Order by new p vals calculated manually 
results <- results[order(results$Age_Pcalc),]
results2 <- results2[order(results2$Sex_Pcalc),]

# Save out results tables for further processing 
write.csv(results, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/Age_sex/age_results_file.csv", row.names = F)
write.csv(results2, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/Age_sex/sex_results_file.csv", row.names = F)



###########################################################################

### REPLICATION ASSESSMENT 

# Chosen studies to include 
# Sun et al - included (N 3,301)
# Ferkingstad et al - included (N 35,559)
# INTERVAL/LonGenity study - lehalier study: https://www.nature.com/articles/s41591-019-0673-2#Sec8 - chosen (N 4,263)

# Remaining studies either had N < 500 or looked at age/sex assocs with DNAm mediation)
# BLSA/GESTALT study - tanaka study: https://onlinelibrary.wiley.com/doi/full/10.1111/acel.12799 - (N 240)
# InCHIANTI study - https://elifesciences.org/articles/61073 (Age associated proteins mediation by DNA methylation - N 997)
# Ngo - https://www.ahajournals.org/doi/full/10.1161/CIRCULATIONAHA.116.021803 (N < 100)
# Menni - https://academic.oup.com/biomedgerontology/article/70/7/809/707747?login=true#supplementary-data (discovery N 206)

### First - look at how many age assocs were also sex significant 

screen

R

age <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/Age_sex/age_results_file.csv")
sex <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/Age_sex/sex_results_file.csv")

library(tidyverse)
library(readxl)

# Merge age and sex results together
agesex <- left_join(age, sex, by = "SeqId")

# Work out how many age assocs were also sex assocs 
age_sig <- filter(agesex, Age_Pcalc < 0.0003496503) # 587
sex_sig <- filter(agesex, Sex_Pcalc < 0.0003496503) # 
agesexsig <- filter(age_sig, Sex_Pcalc < 0.0003496503) # now 222 

### ASSESS REPLICATION OF PREVIOUS RESULTS 

## Ferkingstad et al 

# Accounting for multiple testing, the levels of 63% of the 4,907 proteins correlated positively with 
# age and 18% correlated negatively with age (Supplementary Table 1). Moreover, levels of 33% of the proteins
# were higher in men and 23% were higher in women.

# 5284 total in this study 
# 63 + 18 = 81% significant age
# 33 + 23 = 56% significant sex 

# 5284 / 100, then x 81 = 4280
# 5284 / 100, then x 56 = 2959

thr <- 0.05 / 4907

dat <- read_excel("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Replications/ferkingstad_age_sex.xlsx")
dat <- as.data.frame(dat)

dat$pvalage <- sub(',','.',dat$pvalage)
dat$pvalmale <- sub(',','.',dat$pvalmale)

dat$pvalage <- as.numeric(dat$pvalage)
dat$pvalmale <- as.numeric(dat$pvalmale)

datage <- dat[order(dat$pvalage),]
datsex <- dat[order(dat$pvalmale),]

# check thresholding is correct for  adjustment 
sub<- datage[which(datage$pvalage < thr),]
neg <- sub[which(sub$betaage < 0),] 
dim(neg)
pos <- sub[which(sub$betaage > 0),] 
dim(pos)
perc <- (dim(neg)[1] / 5284 ) * 100
perc
perc2 <- (dim(pos)[1] / 5284 ) * 100
perc2

sub <- datsex[which(datsex$pvalmale < thr),]
sub <- datsex[1:2959,]
neg <- sub[which(sub$betamale < 0),] 
dim(neg)
pos <- sub[which(sub$betamale > 0),] 
dim(pos)
perc <- (dim(neg)[1] / 5284 ) * 100
perc
perc2 <- (dim(pos)[1] / 5284 ) * 100
perc2

# Create subset of significant results to compare to 
datage<- datage[which(datage$pvalage < thr),] # 4306
datsex <- datsex[which(datsex$pvalmale < thr),] # 2921

# Create matchable IDs
datage$SeqId <- sub('_','-',datage$SeqId)
datsex$SeqId <- sub('_','-',datsex$SeqId)

write.csv(datage, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Replications/ferkingstad_age.csv", row.names = F)
write.csv(datsex, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Replications/ferkingstad_sex.csv", row.names = F)


## Sun et al 
sun <- read_excel("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Replications/sun_age_sex_BMI.xlsx")
sun <- as.data.frame(sun)
names(sun)[1] <- "SeqId"

# Create matchable IDs 
# anno <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Annotations_protein_file/annotation_formatted_for_paper.csv")
library(stringr)
x = strsplit(sun$SeqId, ".", fixed = T)
head(x)
ids = lapply(x, function(x){paste0(x[2], "-", x[3])})
ids = unlist(ids)
sun$SeqId <- ids

names(sun)[7:9] <-  c("agebeta", "ageSE", "ageP")
names(sun)[10:12] <- c("femalebeta", "femaleSE", "femaleP")

# convert p vals back from -log10(P)
sun$ageP <- 10^(-sun$ageP)
sun$femaleP <- 10^(-sun$femaleP)

sunage <- sun[order(sun$ageP),]
sunsex <- sun[order(sun$femaleP),]

write.csv(sunage, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Replications/sun_age.csv", row.names = F)
write.csv(sunsex, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Replications/sun_sex.csv", row.names = F)


## Lehallier et al 
# leh <- read_excel("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Replications/INTERVAL_lehalier_results.xlsx")
# leh <- as.data.frame(leh)

# dim(leh) # 2925 - matches paper!

# # # Sort out the naming situation for joining to our results - no longer needed thanks to linker used above 
# library(stringr)
# x = strsplit(leh$ID, ".", fixed = T)
# head(x)
# ids = lapply(x, function(x){paste0(x[2], "-", x[3])})
# ids = unlist(ids)
# leh$SeqId <- ids

# # Now merge in the seqIds and make sure the conversion has been completed for all of them
# anno <- read_excel("Y:/Protein_DNAm_Proxies/Manuscript_revision/Annotations_for_reference.xlsx")
# anno <- as.data.frame(anno)
# anno <- anno[c(1,18,13,4,2)]
# anno$check_col <- anno$SeqId
# anno <- anno[c(1,6,2:5)]
# leh <- left_join(leh, anno, by = "SeqId")
# write.csv(leh, "lehalier_renamed_joint_to_anno_using_extraction_method_chopped.csv", row.names = F) # check the conversions and manually edit any issues 

# # Read in the linker file to sort out the SeqIds 
# link <- read_excel("Y:/Danni/stradl_markers/pheWAS/03_lmekin_rank_inverse_280321/results_age_sex/replication_assessment/lehalier_linker_file.xlsx")
# link$SeqId <- as.character(link$SeqId)
# link <- as.data.frame(link)
# link$SeqId <- as.character(link$SeqId)
# link$SeqId <- gsub("\\.", "-", link$SeqId) # format the SeqIds to be consistent with ours 

# # Join the seqids into the lehalier human results file 
# leh <- left_join(leh, link, by = "ID")

# Read in the correct SeqIds file as edited to remove any issues with naming 
leh <- read_excel("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Replications/lehalier_renamed_joint_to_anno_using_extraction_method_chopped_edited_excel_file.xlsx")
leh <- as.data.frame(leh)
length(unique(leh$SeqId)) # 2925 - all good

# Correct the NA values for 4 issue cells read in 
t <- na.omit(leh)
diff <-  which(!leh$SeqId %in% t$SeqId)
issue <- leh[c(2924:2925),]

# Add in the correct values from the table in place of the incorrect NA values 
issue[1,2] <- 1.97626258336499E-323
issue[1,3] <- 2.89028402817129E-320

issue[2,2] <- 9.88131291682493E-324
issue[2,3] <- 2.89028402817129E-320

# Join the right 4 cells back to main dataset 
leh <- rbind(t, issue) # 2925
leh <- leh[c(14,2:7)]
names(leh) <- c("SeqId", "leh_Age_P", "leh_Age_P_adj", "leh_Age_beta", "leh_Sex_P", "leh_Sex_P_adj", "leh_Sex_beta")
write.csv(leh, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Replications/lehalier_renamed_and_complete.csv", row.names = F)


################

## Read age sex comparison files back into R for each study 

screen

R

library(tidyverse)
library(readxl)

sunage <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Replications/sun_age.csv")
sunsex <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Replications/sun_sex.csv")
ferkage <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Replications/ferkingstad_age.csv")
ferksex <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Replications/ferkingstad_sex.csv")
leh <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Replications/lehalier_renamed_and_complete.csv")

# get sets of significant results
lehage <- leh[which(leh$leh_Age_P_adj <= 0.05),] # 1379 - matches their paper exactly  
lehsex <- leh[which(leh$leh_Sex_P_adj <= 0.05),] # 1651 - matches their paper exactly 

sunage <- sunage[which(sunage$ageP <= 1e-5),] # 662 
sunsex <- sunsex[which(sunsex$femaleP <= 1e-5),] # 1249

dim(ferkage) # 4306 
dim(ferksex) # 2921

# sun - flip effect direction
# our betas are sex(M), lehallier is sex(M) and sun is sex (F)

sunsex[,10][sapply(sunsex[,10], is.numeric)] <- sunsex[,10][sapply(sunsex[,10], is.numeric)] * -1

## Read in our age sex results 

age <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/Age_sex/age_results_file.csv")
sex <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/Age_sex/sex_results_file.csv")

# Merge age and sex results together
agesex <- left_join(age, sex, by = "SeqId")

# Work out how many age assocs were also sex assocs 
age_sig <- filter(agesex, Age_Pcalc < 0.0003496503) # 587
sex_sig <- filter(agesex, Sex_Pcalc < 0.0003496503) # 
agesexsig <- filter(age_sig, Sex_Pcalc < 0.0003496503) # now 222 


## Now index whether there are proteins that replicated across studies 

#age 
age <- age[c(1,6,2,3,8)]
age$Sig <- ifelse(age$Age_Pcalc < 0.0003496503, "yes", "no") # 587 associations in our study 
age$Lehallier_rep <- "no"
age$Ferkingstad_rep <- "no"
age$Sun_rep <- "no"

# age vs lehallier 
age <- left_join(age, lehage, by = "SeqId")
age <- left_join(age, ferkage, by = "SeqId")
age <- left_join(age, sunage, by = "SeqId")

# sign(age$Age_beta[i])==sign(age$leh_Age_beta)[i]

for(i in 1:4235){
  prot <- as.character(age$SeqId[i])
  if(age$Sig[i] == "yes" & prot %in% lehage$SeqId == TRUE & sign(age$Age_beta[i])==sign(age$leh_Age_beta)[i]){
    age[i,7] <- "yes"
  } else {
    print("null")
  }

  if(age$Sig[i] == "yes" & prot %in% ferkage$SeqId == TRUE & sign(age$Age_beta[i])==sign(age$betaage)[i]){
    age[i,8] <- "yes"
  } else {
    print("null")
  }

  if(age$Sig[i] == "yes" & prot %in% sunage$SeqId == TRUE & sign(age$Age_beta[i])==sign(age$agebeta)[i]){
    age[i,9] <- "yes"
  } else {
    print("null")
  }
}


# Work out total number of comparable associations between them:
sig <- age[age$Sig == "yes",]
sig <- sig$SeqId

length(which(sig %in% ferkage$SeqId)) # 584
length(which(sig %in% sunage$SeqId)) # 175
length(which(sig %in% lehage$SeqId)) # 291

# Of 587 assocs, there were this many that replicated in the cohorts:
table(age$Ferkingstad_rep)
table(age$Sun_rep)
table(age$Lehallier_rep)


# before adjusting for direction of betas:
#   no  yes
# 3652  583

# > 
#   no  yes
# 4060  175

# > 
#   no  yes
# 3945  290

# after adjusting for direction of betas:

#   no  yes
# 3702  580

#   no  yes
# 4111  158

#   no  yes
# 3963  273


# so of 1050 comparable assocs, 1021 replicate 
# this is 97% of comparable assocs for age 

## sex
sex <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/Age_sex/sex_results_file.csv")
sex <- sex[c(1,6,2,3,8)]
sex$Sig <- ifelse(sex$Sex_Pcalc < 0.0003496503, "yes", "no") # 545 associations in our study 
sex$Lehallier_rep <- "no"
sex$Ferkingstad_rep <- "no"
sex$Sun_rep <- "no"

lehsex <- lehsex[c(1,5,6,7)]
ferksex <- ferksex[c(1,15,16)]
sunsex <- sunsex[c(1,10,11,12)]

sex <- left_join(sex, lehsex, by = "SeqId")
sex <- left_join(sex, ferksex, by = "SeqId")
sex <- left_join(sex, sunsex, by = "SeqId")


for(i in 1:4235){
  prot <- as.character(sex$SeqId[i])
  if(sex$Sig[i] == "yes" & prot %in% lehsex$SeqId == TRUE & sign(sex$Sex_beta[i])==sign(sex$leh_Sex_beta)[i]){
    sex[i,7] <- "yes"
  } else {
    print("null")
  }

  if(sex$Sig[i] == "yes" & prot %in% ferksex$SeqId == TRUE & sign(sex$Sex_beta[i])==sign(sex$betamale)[i]){
    sex[i,8] <- "yes"
  } else {
    print("null")
  }

  if(sex$Sig[i] == "yes" & prot %in% sunsex$SeqId == TRUE & sign(sex$Sex_beta[i])==sign(sex$femalebeta)[i]){
    sex[i,9] <- "yes"
  } else {
    print("null")
  }
}


# Work out total number of comparable associations between them:
sig <- sex[sex$Sig == "yes",]
sig <- sig$SeqId

length(which(sig %in% ferksex$SeqId)) # 514
length(which(sig %in% sunsex$SeqId)) # 241
length(which(sig %in% lehsex$SeqId)) # 268

# Of 545 assocs, there were this many that replicated in the cohorts:
table(sex$Ferkingstad_rep)
table(sex$Sun_rep)
table(sex$Lehallier_rep)

#   no  yes
# 3732  505

#   no  yes
# 4005  232

#   no  yes
# 3973  264

# Of 1023 possible comparisons, there were 1001 replications
# this is 98% 

write.csv(sex, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Replications/replication_sex.csv", row.names = F)
write.csv(age, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Replications/replication_age.csv", row.names = F)


###################

### JOIN TO SUPPL TABLE AND SAVE 

screen

R

age <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Replications/replication_age.csv")
sex <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Replications/replication_sex.csv")

age <- age[c(1:9)]
sex <- sex[c(1:9)]

sex <- sex[-2]

join <- left_join(age, sex, by = "SeqId")

# Rename and save 
names(join) <- c("SeqId", "Entrez Gene Name",
  "Age Beta", "Age SE", "Age P", "Age Significant",
  "Age Replicated Lehellier", "Age Replicated Ferkingstad", "Age Replicated Sun",
  "Sex Beta", "Sex SE", "Sex P", "Sex Significant",
  "Sex Replicated Lehellier", "Sex Replicated Ferkingstad", "Sex Replicated Sun")

write.csv(join, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Replications/replication_joint.csv", row.names = F)


