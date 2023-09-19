################################################################################################################################

### PheWAS analyses - Cognitive, APOE and Imaging

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

### COGNITIVE

#####################################################################################################################################

## Cognitive scores 
# phenotype <- c("digit_symbol")
# phenotype <- c("verbal_total")
# phenotype <- c( "gf")
# phenotype <- c("g")
# phenotype <- c("LM")
# phenotype <- c("mr_correct")
# phenotype <- c("vocabulary")

length <- 4235
phenotype <- c("digit_symbol")

for (j in 1:length(phenotype)){
  pheno <- phenotype[j]
  results <- data.frame(SeqId = 1:length, phenotype = 1:length, n = 1:length, beta = 1:length, SE = 1:length, P = 1:length)
  for(i in 1:4235){
    tryCatch({ 
      prot_name <- as.character(colnames(markers[i]))
  
      mod <- lmekin(scale(prot[,prot_name]) ~ scale(prot$st_age) + as.factor(prot$sex) + scale(prot[,pheno]) + as.factor(prot$combined) + as.factor(study_site) + as.factor(lag_group) + (1|prot$GS_id), data = prot, varlist = kin_model*2)

      print(i)
      print(prot_name)
      print(pheno)

      results[i,1] <- prot_name
      results[i,2] <- pheno
      results[i,3] <- mod$n
      results[i,4] <- extract_coxme_table(mod)[4,1]
      results[i,5] <- extract_coxme_table(mod)[4,2]
      results[i,6] <- extract_coxme_table(mod)[4,4]
    }, error = function(e) cat("skipped"))
  }
  write.csv(results, paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/Cognitive/", pheno, "_result.csv"))
}


length <- 4235
phenotype <- c("verbal_total")

for (j in 1:length(phenotype)){
  pheno <- phenotype[j]
  results <- data.frame(SeqId = 1:length, phenotype = 1:length, n = 1:length, beta = 1:length, SE = 1:length, P = 1:length)
  for(i in 1:4235){
    tryCatch({ 
      prot_name <- as.character(colnames(markers[i]))
  
      mod <- lmekin(scale(prot[,prot_name]) ~ scale(prot$st_age) + as.factor(prot$sex) + scale(prot[,pheno]) + + as.factor(prot$combined) + as.factor(study_site) + as.factor(lag_group) + (1|prot$GS_id), data = prot, varlist = kin_model*2)

      print(i)
      print(prot_name)
      print(pheno)

      results[i,1] <- prot_name
      results[i,2] <- pheno
      results[i,3] <- mod$n
      results[i,4] <- extract_coxme_table(mod)[4,1]
      results[i,5] <- extract_coxme_table(mod)[4,2]
      results[i,6] <- extract_coxme_table(mod)[4,4]
    }, error = function(e) cat("skipped"))
  }
  write.csv(results, paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/Cognitive/", pheno, "_result.csv"))
}


length <- 4235
phenotype <- c("gf")

for (j in 1:length(phenotype)){
  pheno <- phenotype[j]
  results <- data.frame(SeqId = 1:length, phenotype = 1:length, n = 1:length, beta = 1:length, SE = 1:length, P = 1:length)
  for(i in 1:4235){
    tryCatch({ 
      prot_name <- as.character(colnames(markers[i]))
  
      mod <- lmekin(scale(prot[,prot_name]) ~ scale(prot$st_age) + as.factor(prot$sex) + scale(prot[,pheno]) + + as.factor(prot$combined) + as.factor(study_site) + as.factor(lag_group) + (1|prot$GS_id), data = prot, varlist = kin_model*2)

      print(i)
      print(prot_name)
      print(pheno)

      results[i,1] <- prot_name
      results[i,2] <- pheno
      results[i,3] <- mod$n
      results[i,4] <- extract_coxme_table(mod)[4,1]
      results[i,5] <- extract_coxme_table(mod)[4,2]
      results[i,6] <- extract_coxme_table(mod)[4,4]
    }, error = function(e) cat("skipped"))
  }
  write.csv(results, paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/Cognitive/", pheno, "_result.csv"))
}


length <- 4235
phenotype <- c("g")

for (j in 1:length(phenotype)){
  pheno <- phenotype[j]
  results <- data.frame(SeqId = 1:length, phenotype = 1:length, n = 1:length, beta = 1:length, SE = 1:length, P = 1:length)
  for(i in 1:4235){
    tryCatch({ 
      prot_name <- as.character(colnames(markers[i]))
  
      mod <- lmekin(scale(prot[,prot_name]) ~ scale(prot$st_age) + as.factor(prot$sex) + scale(prot[,pheno]) + + as.factor(prot$combined) + as.factor(study_site) + as.factor(lag_group) + (1|prot$GS_id), data = prot, varlist = kin_model*2)

      print(i)
      print(prot_name)
      print(pheno)

      results[i,1] <- prot_name
      results[i,2] <- pheno
      results[i,3] <- mod$n
      results[i,4] <- extract_coxme_table(mod)[4,1]
      results[i,5] <- extract_coxme_table(mod)[4,2]
      results[i,6] <- extract_coxme_table(mod)[4,4]
    }, error = function(e) cat("skipped"))
  }
  write.csv(results, paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/Cognitive/", pheno, "_result.csv"))
}


length <- 4235
phenotype <- c("LM")

for (j in 1:length(phenotype)){
  pheno <- phenotype[j]
  results <- data.frame(SeqId = 1:length, phenotype = 1:length, n = 1:length, beta = 1:length, SE = 1:length, P = 1:length)
  for(i in 1:4235){
    tryCatch({ 
      prot_name <- as.character(colnames(markers[i]))
  
      mod <- lmekin(scale(prot[,prot_name]) ~ scale(prot$st_age) + as.factor(prot$sex) + scale(prot[,pheno]) + + as.factor(prot$combined) + as.factor(study_site) + as.factor(lag_group) + (1|prot$GS_id), data = prot, varlist = kin_model*2)

      print(i)
      print(prot_name)
      print(pheno)

      results[i,1] <- prot_name
      results[i,2] <- pheno
      results[i,3] <- mod$n
      results[i,4] <- extract_coxme_table(mod)[4,1]
      results[i,5] <- extract_coxme_table(mod)[4,2]
      results[i,6] <- extract_coxme_table(mod)[4,4]
    }, error = function(e) cat("skipped"))
  }
  write.csv(results, paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/Cognitive/", pheno, "_result.csv"))
}


length <- 4235
phenotype <- c("mr_correct")

for (j in 1:length(phenotype)){
  pheno <- phenotype[j]
  results <- data.frame(SeqId = 1:length, phenotype = 1:length, n = 1:length, beta = 1:length, SE = 1:length, P = 1:length)
  for(i in 1:4235){
    tryCatch({ 
      prot_name <- as.character(colnames(markers[i]))
  
      mod <- lmekin(scale(prot[,prot_name]) ~ scale(prot$st_age) + as.factor(prot$sex) + scale(prot[,pheno]) + + as.factor(prot$combined) + as.factor(study_site) + as.factor(lag_group) + (1|prot$GS_id), data = prot, varlist = kin_model*2)

      print(i)
      print(prot_name)
      print(pheno)

      results[i,1] <- prot_name
      results[i,2] <- pheno
      results[i,3] <- mod$n
      results[i,4] <- extract_coxme_table(mod)[4,1]
      results[i,5] <- extract_coxme_table(mod)[4,2]
      results[i,6] <- extract_coxme_table(mod)[4,4]
    }, error = function(e) cat("skipped"))
  }
  write.csv(results, paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/Cognitive/", pheno, "_result.csv"))
}


length <- 4235
phenotype <- c("vocabulary")

for (j in 1:length(phenotype)){
  pheno <- phenotype[j]
  results <- data.frame(SeqId = 1:length, phenotype = 1:length, n = 1:length, beta = 1:length, SE = 1:length, P = 1:length)
  for(i in 1:4235){
    tryCatch({ 
      prot_name <- as.character(colnames(markers[i]))
  
      mod <- lmekin(scale(prot[,prot_name]) ~ scale(prot$st_age) + as.factor(prot$sex) + scale(prot[,pheno]) + + as.factor(prot$combined) + as.factor(study_site) + as.factor(lag_group) + (1|prot$GS_id), data = prot, varlist = kin_model*2)

      print(i)
      print(prot_name)
      print(pheno)

      results[i,1] <- prot_name
      results[i,2] <- pheno
      results[i,3] <- mod$n
      results[i,4] <- extract_coxme_table(mod)[4,1]
      results[i,5] <- extract_coxme_table(mod)[4,2]
      results[i,6] <- extract_coxme_table(mod)[4,4]
    }, error = function(e) cat("skipped"))
  }
  write.csv(results, paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/Cognitive/", pheno, "_result.csv"))
}


##############################################################################################################

### RUN APOE MODELS 

##############################################################################################################

# Now we can run APOE status through the loops as numeric continuum 
length <- 4235
set <- c("set")
phenotype <- c("APOE")

for (j in 1:length(set)){
  pheno <- "APOE"
  results <- data.frame(SeqId = 1:length, phenotype = 1:length, n = 1:length, beta = 1:length, SE = 1:length, P = 1:length)
  for(i in 1:1000){
    tryCatch({ 
      prot_name <- as.character(colnames(markers[i]))

      mod <- lmekin(scale(prot[,prot_name]) ~ scale(prot$st_age) + as.factor(prot$sex) + scale(prot[,pheno]) + as.factor(prot$combined) + as.factor(study_site) + as.factor(lag_group) + (1|prot$GS_id), data = prot, varlist = kin_model*2)

      print(i)
      print(prot_name)
      print(pheno)

      results[i,1] <- prot_name
      results[i,2] <- paste0(pheno)
      results[i,3] <- mod$n
      results[i,4] <- extract_coxme_table(mod)[4,1]
      results[i,5] <- extract_coxme_table(mod)[4,2]
      results[i,6] <- extract_coxme_table(mod)[4,4]

    }, error = function(e) cat("skipped"))
  }
  write.csv(results, paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/APOE/batches/e2_b1.csv"), row.names = F)
}


# Now we can run APOE status through the loops as numeric continuum 
length <- 4235
set <- c("set")
phenotype <- c("APOE")

for (j in 1:length(set)){
  pheno <- "APOE"
  results <- data.frame(SeqId = 1:length, phenotype = 1:length, n = 1:length, beta = 1:length, SE = 1:length, P = 1:length)
  for(i in 1001:2000){
    tryCatch({ 
      prot_name <- as.character(colnames(markers[i]))

      mod <- lmekin(scale(prot[,prot_name]) ~ scale(prot$st_age) + as.factor(prot$sex) + scale(prot[,pheno]) + as.factor(prot$combined) + as.factor(study_site) + as.factor(lag_group) + (1|prot$GS_id), data = prot, varlist = kin_model*2)

      print(i)
      print(prot_name)
      print(pheno)

      results[i,1] <- prot_name
      results[i,2] <- paste0(pheno)
      results[i,3] <- mod$n
      results[i,4] <- extract_coxme_table(mod)[4,1]
      results[i,5] <- extract_coxme_table(mod)[4,2]
      results[i,6] <- extract_coxme_table(mod)[4,4]

    }, error = function(e) cat("skipped"))
  }
  write.csv(results, paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/APOE/batches/e2_b2.csv"), row.names = F)
}


# Now we can run APOE status through the loops as numeric continuum 
length <- 4235
set <- c("set")
phenotype <- c("APOE")

for (j in 1:length(set)){
  pheno <- "APOE"
  results <- data.frame(SeqId = 1:length, phenotype = 1:length, n = 1:length, beta = 1:length, SE = 1:length, P = 1:length)
  for(i in 2001:3000){
    tryCatch({ 
      prot_name <- as.character(colnames(markers[i]))

      mod <- lmekin(scale(prot[,prot_name]) ~ scale(prot$st_age) + as.factor(prot$sex) + scale(prot[,pheno]) + as.factor(prot$combined) + as.factor(study_site) + as.factor(lag_group) + (1|prot$GS_id), data = prot, varlist = kin_model*2)

      print(i)
      print(prot_name)
      print(pheno)

      results[i,1] <- prot_name
      results[i,2] <- paste0(pheno)
      results[i,3] <- mod$n
      results[i,4] <- extract_coxme_table(mod)[4,1]
      results[i,5] <- extract_coxme_table(mod)[4,2]
      results[i,6] <- extract_coxme_table(mod)[4,4]

    }, error = function(e) cat("skipped"))
  }
  write.csv(results, paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/APOE/batches/e2_b3.csv"), row.names = F)
}


# Now we can run APOE status through the loops as numeric continuum 
length <- 4235
set <- c("set")
phenotype <- c("APOE")

for (j in 1:length(set)){
  pheno <- "APOE"
  results <- data.frame(SeqId = 1:length, phenotype = 1:length, n = 1:length, beta = 1:length, SE = 1:length, P = 1:length)
  for(i in 3001:4000){
    tryCatch({ 
      prot_name <- as.character(colnames(markers[i]))

      mod <- lmekin(scale(prot[,prot_name]) ~ scale(prot$st_age) + as.factor(prot$sex) + scale(prot[,pheno]) + as.factor(prot$combined) + as.factor(study_site) + as.factor(lag_group) + (1|prot$GS_id), data = prot, varlist = kin_model*2)

      print(i)
      print(prot_name)
      print(pheno)

      results[i,1] <- prot_name
      results[i,2] <- paste0(pheno)
      results[i,3] <- mod$n
      results[i,4] <- extract_coxme_table(mod)[4,1]
      results[i,5] <- extract_coxme_table(mod)[4,2]
      results[i,6] <- extract_coxme_table(mod)[4,4]

    }, error = function(e) cat("skipped"))
  }
  write.csv(results, paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/APOE/batches/e2_b4.csv"), row.names = F)
}


# Now we can run APOE status through the loops as numeric continuum 
length <- 4235
set <- c("set")
phenotype <- c("APOE")

for (j in 1:length(set)){
  pheno <- "APOE"
  results <- data.frame(SeqId = 1:length, phenotype = 1:length, n = 1:length, beta = 1:length, SE = 1:length, P = 1:length)
  for(i in 4001:4235){
    tryCatch({ 
      prot_name <- as.character(colnames(markers[i]))

      mod <- lmekin(scale(prot[,prot_name]) ~ scale(prot$st_age) + as.factor(prot$sex) + scale(prot[,pheno]) + as.factor(prot$combined) + as.factor(study_site) + as.factor(lag_group) + (1|prot$GS_id), data = prot, varlist = kin_model*2)

      print(i)
      print(prot_name)
      print(pheno)

      results[i,1] <- prot_name
      results[i,2] <- paste0(pheno)
      results[i,3] <- mod$n
      results[i,4] <- extract_coxme_table(mod)[4,1]
      results[i,5] <- extract_coxme_table(mod)[4,2]
      results[i,6] <- extract_coxme_table(mod)[4,4]

    }, error = function(e) cat("skipped"))
  }
  write.csv(results, paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/APOE/batches/e2_b5.csv"), row.names = F)
}


##############################################################################################################

### RUN IMAGING

##############################################################################################################

length <- 4235
phenotype <- c("brain_accel")

for (j in 1:length(phenotype)){
  pheno <- phenotype[j]
  results <- data.frame(SeqId = 1:length, phenotype = 1:length, n = 1:length, beta = 1:length, SE = 1:length, P = 1:length)
  for(i in 1:4235){
    tryCatch({ 
      prot_name <- as.character(colnames(markers[i]))
  
      mod <- lmekin(scale(prot[,prot_name]) ~ scale(prot$st_age) + as.factor(prot$sex) + scale(prot[,pheno]) + as.factor(prot$combined) + as.factor(study_site) + as.factor(lag_group) + (1|prot$GS_id), na.action = na.exclude, data = prot, varlist = kin_model*2)

      print(i)
      print(prot_name)
      print(pheno)

      results[i,1] <- prot_name
      results[i,2] <- pheno
      results[i,3] <- mod$n
      results[i,4] <- extract_coxme_table(mod)[4,1]
      results[i,5] <- extract_coxme_table(mod)[4,2]
      results[i,6] <- extract_coxme_table(mod)[4,4]
    }, error = function(e) cat("skipped"))
  }
  write.csv(results, paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/Imaging/", pheno, "_result.csv"))
}


length <- 4235
phenotype <- c("gFA")

for (j in 1:length(phenotype)){
  pheno <- phenotype[j]
  results <- data.frame(SeqId = 1:length, phenotype = 1:length, n = 1:length, beta = 1:length, SE = 1:length, P = 1:length)
  for(i in 1:4235){
    tryCatch({ 
      prot_name <- as.character(colnames(markers[i]))
  
      mod <- lmekin(scale(prot[,prot_name]) ~ scale(prot$st_age) + as.factor(prot$sex) + scale(prot[,pheno]) + as.factor(prot$combined) + as.factor(study_site) + as.factor(lag_group) + (1|prot$GS_id), na.action = na.exclude, data = prot, varlist = kin_model*2)

      print(i)
      print(prot_name)
      print(pheno)

      results[i,1] <- prot_name
      results[i,2] <- pheno
      results[i,3] <- mod$n
      results[i,4] <- extract_coxme_table(mod)[4,1]
      results[i,5] <- extract_coxme_table(mod)[4,2]
      results[i,6] <- extract_coxme_table(mod)[4,4]
    }, error = function(e) cat("skipped"))
  }
  write.csv(results, paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/Imaging/", pheno, "_result.csv"))
}


length <- 4235
phenotype <- c("gMD")

for (j in 1:length(phenotype)){
  pheno <- phenotype[j]
  results <- data.frame(SeqId = 1:length, phenotype = 1:length, n = 1:length, beta = 1:length, SE = 1:length, P = 1:length)
  for(i in 1:4235){
    tryCatch({ 
      prot_name <- as.character(colnames(markers[i]))
  
      mod <- lmekin(scale(prot[,prot_name]) ~ scale(prot$st_age) + as.factor(prot$sex) + scale(prot[,pheno]) + as.factor(prot$combined) + as.factor(study_site) + as.factor(lag_group) + (1|prot$GS_id), na.action = na.exclude, data = prot, varlist = kin_model*2)

      print(i)
      print(prot_name)
      print(pheno)

      results[i,1] <- prot_name
      results[i,2] <- pheno
      results[i,3] <- mod$n
      results[i,4] <- extract_coxme_table(mod)[4,1]
      results[i,5] <- extract_coxme_table(mod)[4,2]
      results[i,6] <- extract_coxme_table(mod)[4,4]
    }, error = function(e) cat("skipped"))
  }
  write.csv(results, paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/Imaging/", pheno, "_result_updated_gMD.csv"))
}




length <- 4235
phenotype <- c("Fazekas_Score_Total")

for (j in 1:length(phenotype)){
  pheno <- phenotype[j]
  results <- data.frame(SeqId = 1:length, phenotype = 1:length, n = 1:length, beta = 1:length, SE = 1:length, P = 1:length)
  for(i in 1:4235){
    tryCatch({ 
      prot_name <- as.character(colnames(markers[i]))
  
      mod <- lmekin(scale(prot[,prot_name]) ~ scale(prot$st_age) + as.factor(prot$sex) + scale(prot[,pheno]) + as.factor(prot$combined) + as.factor(study_site) + as.factor(lag_group) + (1|prot$GS_id), na.action = na.exclude, data = prot, varlist = kin_model*2)

      print(i)
      print(prot_name)
      print(pheno)

      results[i,1] <- prot_name
      results[i,2] <- pheno
      results[i,3] <- mod$n
      results[i,4] <- extract_coxme_table(mod)[4,1]
      results[i,5] <- extract_coxme_table(mod)[4,2]
      results[i,6] <- extract_coxme_table(mod)[4,4]
    }, error = function(e) cat("skipped"))
  }
  write.csv(results, paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/Imaging/", pheno, "_result.csv"))
}


length <- 4235
phenotype <- c("WBV_No_Ventricles")

for (j in 1:length(phenotype)){
  pheno <- phenotype[j]
  results <- data.frame(SeqId = 1:length, phenotype = 1:length, n = 1:length, beta = 1:length, SE = 1:length, P = 1:length)
  for(i in 1:4235){
    tryCatch({ 
      prot_name <- as.character(colnames(markers[i]))
  
      mod <- lmekin(scale(prot[,prot_name]) ~ scale(prot$st_age) + as.factor(prot$sex) + scale(prot[,pheno]) + as.factor(prot$combined) + as.factor(prot$Site)*prot$Estimated_ICV + as.factor(lag_group) + as.factor(prot$Edited) + as.factor(prot$Batch) + (1|prot$GS_id), na.action = na.exclude, data = prot, varlist = kin_model*2)

      print(i)
      print(prot_name)
      print(pheno)

      results[i,1] <- prot_name
      results[i,2] <- pheno
      results[i,3] <- mod$n
      results[i,4] <- extract_coxme_table(mod)[4,1]
      results[i,5] <- extract_coxme_table(mod)[4,2]
      results[i,6] <- extract_coxme_table(mod)[4,4]
    }, error = function(e) cat("skipped"))
  }
  write.csv(results, paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/Imaging/", pheno, "_result.csv"))
}



length <- 4235
phenotype <- c("Global_GM_Volume")

for (j in 1:length(phenotype)){
  pheno <- phenotype[j]
  results <- data.frame(SeqId = 1:length, phenotype = 1:length, n = 1:length, beta = 1:length, SE = 1:length, P = 1:length)
  for(i in 1:4235){
    tryCatch({ 
      prot_name <- as.character(colnames(markers[i]))
  
      mod <- lmekin(scale(prot[,prot_name]) ~ scale(prot$st_age) + as.factor(prot$sex) + scale(prot[,pheno]) + as.factor(prot$combined) + as.factor(prot$Site)*prot$Estimated_ICV + as.factor(lag_group) + as.factor(prot$Edited) + as.factor(prot$Batch) + (1|prot$GS_id), na.action = na.exclude, data = prot, varlist = kin_model*2)

      print(i)
      print(prot_name)
      print(pheno)

      results[i,1] <- prot_name
      results[i,2] <- pheno
      results[i,3] <- mod$n
      results[i,4] <- extract_coxme_table(mod)[4,1]
      results[i,5] <- extract_coxme_table(mod)[4,2]
      results[i,6] <- extract_coxme_table(mod)[4,4]
    }, error = function(e) cat("skipped"))
  }
  write.csv(results, paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/Imaging/", pheno, "_result.csv"))
}


length <- 4235
phenotype <- c("WMH_Volume_Total")

for (j in 1:length(phenotype)){
  pheno <- phenotype[j]
  results <- data.frame(SeqId = 1:length, phenotype = 1:length, n = 1:length, beta = 1:length, SE = 1:length, P = 1:length)
  for(i in 1:4235){
    tryCatch({ 
      prot_name <- as.character(colnames(markers[i]))
  
      mod <- lmekin(scale(prot[,prot_name]) ~ scale(prot$st_age) + as.factor(prot$sex) + scale(prot[,pheno]) + as.factor(prot$combined) + as.factor(prot$Site)*prot$Estimated_ICV + as.factor(lag_group) + as.factor(prot$Reviewer) + (1|prot$GS_id), na.action = na.exclude, data = prot, varlist = kin_model*2)

      print(i)
      print(prot_name)
      print(pheno)

      results[i,1] <- prot_name
      results[i,2] <- pheno
      results[i,3] <- mod$n
      results[i,4] <- extract_coxme_table(mod)[4,1]
      results[i,5] <- extract_coxme_table(mod)[4,2]
      results[i,6] <- extract_coxme_table(mod)[4,4]
    }, error = function(e) cat("skipped"))
  }
  write.csv(results, paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/Imaging/", pheno, "_result_updated_WMHV_with_ICV_and_site_and_editor.csv"))
}

# 27842.pts-11.ccace-proc12       (Attached)
#         27738.pts-8.ccace-proc12        (Attached)
#         26339.pts-3.ccace-proc12        (Detached)
#         2720.pts-3.ccace-proc12 (Detached)
#         25132.pts-33.ccace-proc12       (Detached)
