##############################################################################################################

### PROCESS PheWAS models

##############################################################################################################

## This script joins batch PheWAS files together for each brain health trait (APOE, cognitive, imaging)
## It then loops through the joint files to annotate them with protein information and significance thresholding
## Any comparison lookups are conducted with previous studies of relevance
## Finally, results tables are joint and saved for supplementary files

##############################################################################################################

### APOE EXTRACTION

##############################################################################################################

cd /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS/00_Results_collation/joint_files_for_annotation/

screen

R

library(tidyverse)
library(readxl)

## APOE

# Will need to join up the batches of APOE runs for e2 and e4 
e21 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/APOE/batches/e2_b1.csv")
e22 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/APOE/batches/e2_b2.csv")
e23 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/APOE/batches/e2_b3.csv")
e24 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/APOE/batches/e2_b4.csv")
e25 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/APOE/batches/e2_b5.csv")

# Cut the info you need and bind together (should have put NA's in but too late now)
e21 <- e21[1:1000,]
e22 <- e22[1001:2000,]
e23 <- e23[2001:3000,]
e24 <- e24[3001:4000,]
e25 <- e25[4001:4235,]

e2 <- rbind(e21, e22)
e2 <- rbind(e2, e23)
e2 <- rbind(e2, e24)
e2 <- rbind(e2, e25)

# Save out combined files into main repowhere processing code can run the annotations and pass/fail FDR 
write.csv(e2, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/APOE/APOE_JOINT.csv")

# Run the loop to generate annotated results files above with only these 2 APOE files in the folder and youll have them added to results

location <- "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/APOE/"

location_output <- "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/APOE/APOE_processed/"
loop = list.files(paste0(location), ".")

pheno_count <- data.frame(Phenotype = 1:9, Associations = 1:9, Protein_Gene = 1:9)
pheno_count2 <- data.frame(Phenotype = 1:9, Associations = 1:9, Protein_Gene = 1:9)

# for(i in 1){

i <- 1
  ### DO FDR version
  j <- loop[i]
  title = gsub(".csv", "", j)
  title <- as.character(title)

  results <- read.csv(paste0(location, j))
  results <- results[-1]

  # manually calculate P values
  results$Pcalc <- pchisq((results$beta/results$SE)^2, 1, lower.tail=F)

  results <- results[order(results$Pcalc),] # order by P value 
  results$P_adjust <- p.adjust(results$Pcalc, method = "BH") # FDR adjustment 
  results$Status <- ifelse(results$P_adjust < 0.05, "pass", "fail") # Add identifier for those that pass/fail significance in each 
  results <- results[-6]

  anno <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Annotations_protein_file/annotation_formatted_for_paper.csv")

  joint <- left_join(results, anno, by = "SeqId") # join in annotations
  joint <- joint[c(1,10,2:8,9,11)] # order results 
  write.csv(joint, paste0(location_output, title, "_annotated_FDR.csv"), row.names = F)

  sig <- which(results$P_adjust < 0.05) # count those that pass 
  count <- length(sig)[1]
  pheno_count[i,1] <- title
  pheno_count[i,2] <- count

  sign <- joint[sig,] # filter to those with just significant P adjusted values 
  names(sign)[2] <- "Gene_name"
  proteins <- sign$Gene_name
  # proteins <- sign$SeqId
  str <- str_c(proteins, collapse = ", ")
  pheno_count[i,3] <- str
  print(dim(results))

  ############################################################################

  j <- loop[i]
  title = gsub(".csv", "", j)
  title <- as.character(title)

  results <- read.csv(paste0(location, j))
  results <- results[-1]

  # manually calculate P values
  results$Pcalc <- pchisq((results$beta/results$SE)^2, 1, lower.tail=F)

  results <- results[order(results$Pcalc),] # order by P value 
  results$Status <- ifelse(results$Pcalc < 0.0003496503, "pass", "fail") # Add identifier for those that pass/fail significance in each 
  results <- results[-6]

  anno <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Annotations_protein_file/annotation_formatted_for_paper.csv")

  joint <- left_join(results, anno, by = "SeqId") # join in annotations
  joint <- joint[c(1,9,2:8,10)] # order results 
  write.csv(joint, paste0(location_output, title, "_annotated_THR.csv"), row.names = F)

  sig <- which(results$Pcalc < 0.0003496503) # count those that pass 
  count <- length(sig)[1]
  pheno_count2[i,1] <- title
  pheno_count2[i,2] <- count

  sign <- joint[sig,] # filter to those with just significant P adjusted values 
  names(sign)[2] <- "Gene_name"
  proteins <- sign$Gene_name
  # proteins <- sign$SeqId
  str <- str_c(proteins, collapse = ", ")
  pheno_count2[i,3] <- str
  print(dim(results))
# }

write.csv(pheno_count, paste0(location_output, "APOE_FDR_gene_id.csv"), row.names = F)
write.csv(pheno_count2, paste0(location_output, "APOE_THR_gene_id.csv"), row.names = F)

### APOE format for suppl tables 
APOE <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/APOE/APOE_processed/APOE_JOINT_annotated_THR.csv")
APOE <- APOE[c(1,2,10,3:7)]
names(APOE) <- c("SeqId", "Gene of Protein", "UniProt Full Name", "Phenotype", "n", "Beta", "SE", "P")
write.csv(APOE, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/00_Results_collation/APOE_JOINT_formatted.csv", row.names = F)


##############################################################################################################

### COGNITIVE EXTRACTION

##############################################################################################################

cd /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS/03_Cognitive/batches/

screen

R

library(tidyverse)
library(readxl)

# Run the loop to generate annotated results files above with only these 2 APOE files in the folder and youll have them added to results

location <- "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/Cognitive/"

location_output <- "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/00_Results_collation/cognitive/"
loop = list.files(paste0(location), ".")

pheno_count <- data.frame(Phenotype = 1:9, Associations = 1:9, Protein_Gene = 1:9)
pheno_count2 <- data.frame(Phenotype = 1:9, Associations = 1:9, Protein_Gene = 1:9)

for(i in 1:length(loop)){
# i <- 2
  ### DO FDR version
  j <- loop[i]
  title = gsub(".csv", "", j)
  title <- as.character(title)

  results <- read.csv(paste0(location, j))
  results <- results[-1]

  # manually calculate P values
  results$Pcalc <- pchisq((results$beta/results$SE)^2, 1, lower.tail=F)

  results <- results[order(results$Pcalc),] # order by P value 
  results$P_adjust <- p.adjust(results$Pcalc, method = "BH") # FDR adjustment 
  results$Status <- ifelse(results$P_adjust < 0.05, "pass", "fail") # Add identifier for those that pass/fail significance in each 
  results <- results[-6]

  anno <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Annotations_protein_file/annotation_formatted_for_paper.csv")

  joint <- left_join(results, anno, by = "SeqId") # join in annotations
  joint <- joint[c(1,10,2:8,9,11)] # order results 
  write.csv(joint, paste0(location_output, title, "_annotated_FDR.csv"), row.names = F)

  sig <- which(results$P_adjust < 0.05) # count those that pass 
  count <- length(sig)[1]
  pheno_count[i,1] <- title
  pheno_count[i,2] <- count

  sign <- joint[sig,] # filter to those with just significant P adjusted values 
  names(sign)[2] <- "Gene_name"
  proteins <- sign$Gene_name
  # proteins <- sign$SeqId
  str <- str_c(proteins, collapse = ", ")
  pheno_count[i,3] <- str
  print(dim(results))
  print(dim(sign))

  ############################################################################

  ### 0.05 / 239 = 0.000209205 at 85% componenets 
  j <- loop[i]
  title = gsub(".csv", "", j)
  title <- as.character(title)

  results <- read.csv(paste0(location, j))
  results <- results[-1]

  # manually calculate P values
  results$Pcalc <- pchisq((results$beta/results$SE)^2, 1, lower.tail=F)

  results <- results[order(results$Pcalc),] # order by P value 
  results$Status <- ifelse(results$Pcalc < 0.0003496503, "pass", "fail") # Add identifier for those that pass/fail significance in each 
  results <- results[-6]

  anno <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Annotations_protein_file/annotation_formatted_for_paper.csv")

  joint <- left_join(results, anno, by = "SeqId") # join in annotations
  joint <- joint[c(1,9,2:8,10)] # order results 
  write.csv(joint, paste0(location_output, title, "_annotated_THR.csv"), row.names = F)

  sig <- which(results$Pcalc < 0.0003496503) # count those that pass 
  count <- length(sig)[1]
  pheno_count2[i,1] <- title
  pheno_count2[i,2] <- count

  sign <- joint[sig,] # filter to those with just significant P adjusted values 
  names(sign)[2] <- "Gene_name"
  proteins <- sign$Gene_name
  # proteins <- sign$SeqId
  str <- str_c(proteins, collapse = ", ")
  pheno_count2[i,3] <- str
  print(dim(results))
  print(dim(sign)[1])
}

write.csv(pheno_count, paste0(location_output, "COG_FDR_gene_id.csv"), row.names = F)
write.csv(pheno_count2, paste0(location_output, "COG_THR_gene_id.csv"), row.names = F)

### COG format for suppl tables 
digit <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/00_Results_collation/cognitive/digit_symbol_result_annotated_THR.csv")
digit <- digit[c(1,2,10,3:7)]
digit$phenotype <- "Processing Speed - Digit Symbol Score"
names(digit) <- c("SeqId", "Gene of Protein", "UniProt Full Name", "Phenotype 1", "n", "Beta", "SE", "P")

g <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/00_Results_collation/cognitive/g_result_annotated_THR.csv")
g <- g[c(1,3:7)]
g$phenotype <- "General Cognitive Ability - g Score"
names(g) <- c("SeqId", "Phenotype 2", "n", "Beta", "SE", "P")
g <- g[match(digit$SeqId, g$SeqId),]
g <- g[-1]

gf <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/00_Results_collation/cognitive/gf_result_annotated_THR.csv")
gf <- gf[c(1,3:7)]
gf$phenotype <- "General Fluid Cognitive Ability - gf Score"
names(gf) <- c("SeqId", "Phenotype 2", "n", "Beta", "SE", "P")
gf <- gf[match(digit$SeqId, gf$SeqId),]
gf <- gf[-1]

LM <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/00_Results_collation/cognitive/LM_result_annotated_THR.csv")
LM <- LM[c(1,3:7)]
LM$phenotype <- "Logical Memory - Logical Memory Score"
names(LM) <- c("SeqId", "Phenotype 2", "n", "Beta", "SE", "P")
LM <- LM[match(digit$SeqId, LM$SeqId),]
LM <- LM[-1]

NV <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/00_Results_collation/cognitive/mr_correct_result_annotated_THR.csv")
NV <- NV[c(1,3:7)]
NV$phenotype <- "Non-Verbal Reasoning - Matrix Reasoning Score"
names(NV) <- c("SeqId", "Phenotype 2", "n", "Beta", "SE", "P")
NV <- NV[match(digit$SeqId, NV$SeqId),]
NV <- NV[-1]

ver <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/00_Results_collation/cognitive/verbal_total_result_annotated_THR.csv")
ver <- ver[c(1,3:7)]
ver$phenotype <- "Verbal Reasoning - Verbal Total Score"
names(ver) <- c("SeqId", "Phenotype 2", "n", "Beta", "SE", "P")
ver <- ver[match(digit$SeqId, ver$SeqId),]
ver <- ver[-1]

voc <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/00_Results_collation/cognitive/vocabulary_result_annotated_THR.csv")
voc <- voc[c(1,3:7)]
voc$phenotype <- "Vocabulary - Vocabulary Score"
names(voc) <- c("SeqId", "Phenotype 2", "n", "Beta", "SE", "P")
voc <- voc[match(digit$SeqId, voc$SeqId),]
voc <- voc[-1]


# Join together MATCHED files 
join <- cbind(digit, g)
join <- cbind(join, gf)
join <- cbind(join, LM)
join <- cbind(join, NV)
join <- cbind(join, ver)
join <- cbind(join, voc)

# Write out suppl table file 
write.csv(join, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/00_Results_collation/COG_joint_suppl_table.csv", row.names = F)


##############################################################################################################

### PROCESS THE RESULTS FOR IMAGING MODELS

##############################################################################################################

# cd /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS/05_Imaging/batches/

# screen

# R

library(tidyverse)
library(readxl)

# Run the loop to generate annotated results files above with only these 2 APOE files in the folder and youll have them added to results

location <- "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/Imaging/"

location_output <- "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/00_Results_collation/imaging/"
loop = list.files(paste0(location), ".")

pheno_count <- data.frame(Phenotype = 1:9, Associations = 1:9, Protein_Gene = 1:9)
pheno_count2 <- data.frame(Phenotype = 1:9, Associations = 1:9, Protein_Gene = 1:9)

for(i in 1:length(loop)){
# i <- 2
  ### DO FDR version
  j <- loop[i]
  title = gsub(".csv", "", j)
  title <- as.character(title)

  results <- read.csv(paste0(location, j))
  results <- results[-1]

  # manually calculate P values
  results$Pcalc <- pchisq((results$beta/results$SE)^2, 1, lower.tail=F)

  results <- results[order(results$Pcalc),] # order by P value 
  results$P_adjust <- p.adjust(results$Pcalc, method = "BH") # FDR adjustment 
  results$Status <- ifelse(results$P_adjust < 0.05, "pass", "fail") # Add identifier for those that pass/fail significance in each 
  results <- results[-6]

  anno <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Annotations_protein_file/annotation_formatted_for_paper.csv")

  joint <- left_join(results, anno, by = "SeqId") # join in annotations
  joint <- joint[c(1,10,2:8,9,11)] # order results 
  write.csv(joint, paste0(location_output, title, "_annotated_FDR.csv"), row.names = F)

  sig <- which(results$P_adjust < 0.05) # count those that pass 
  count <- length(sig)[1]
  pheno_count[i,1] <- title
  pheno_count[i,2] <- count

  sign <- joint[sig,] # filter to those with just significant P adjusted values 
  names(sign)[2] <- "Gene_name"
  proteins <- sign$Gene_name
  # proteins <- sign$SeqId
  str <- str_c(proteins, collapse = ", ")
  pheno_count[i,3] <- str
  print(dim(results))
  print(dim(sign))

  ############################################################################

  ### 0.05 / 239 = 0.000209205 at 85% componenets 
  j <- loop[i]
  title = gsub(".csv", "", j)
  title <- as.character(title)

  results <- read.csv(paste0(location, j))
  results <- results[-1]

  # manually calculate P values
  results$Pcalc <- pchisq((results$beta/results$SE)^2, 1, lower.tail=F)

  results <- results[order(results$Pcalc),] # order by P value 
  results$Status <- ifelse(results$Pcalc < 0.0003496503, "pass", "fail") # Add identifier for those that pass/fail significance in each 
  results <- results[-6]

  anno <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Annotations_protein_file/annotation_formatted_for_paper.csv")

  joint <- left_join(results, anno, by = "SeqId") # join in annotations
  joint <- joint[c(1,9,2:8,10)] # order results 
  write.csv(joint, paste0(location_output, title, "_annotated_THR.csv"), row.names = F)

  sig <- which(results$Pcalc < 0.0003496503) # count those that pass 
  count <- length(sig)[1]
  pheno_count2[i,1] <- title
  pheno_count2[i,2] <- count

  sign <- joint[sig,] # filter to those with just significant P adjusted values 
  names(sign)[2] <- "Gene_name"
  proteins <- sign$Gene_name
  # proteins <- sign$SeqId
  str <- str_c(proteins, collapse = ", ")
  pheno_count2[i,3] <- str
  print(dim(results))
  print(dim(sign)[1])
}

write.csv(pheno_count, paste0(location_output, "IMG_FDR_gene_id.csv"), row.names = F)
write.csv(pheno_count2, paste0(location_output, "IMG_THR_gene_id.csv"), row.names = F)

######################################################################################

### COLLATE TO A SUMMARY TABLE FOR SUPPL RESULTS 

# Brain age accel 
# Global GM vol
# gFA
# gMD
# Cerebrum total vol 
# Fazekas score total 
# WBV with ventricles 

# We want 1-4235 proteins on the left, then each table with phenotype joined in on the right in layers in this order 

brain <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/00_Results_collation/imaging/brain_accel_result_annotated_THR.csv")
brain <- brain[c(1,2,10,3:7)]
brain$phenotype <- "Brain Age Acceleration"
names(brain) <- c("SeqId", "Gene of Protein", "UniProt Full Name", "Phenotype 1", "n", "Beta", "SE", "P")

GGM <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/00_Results_collation/imaging/Global_GM_Volume_result_annotated_THR.csv")
GGM <- GGM[c(1,3:7)]
GGM$phenotype <- "Global Grey Matter Volume"
names(GGM) <- c("SeqId", "Phenotype 2", "n", "Beta", "SE", "P")
GGM <- GGM[match(brain$SeqId, GGM$SeqId),]
GGM <- GGM[-1]

GFA <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/00_Results_collation/imaging/gFA_result_annotated_THR.csv")
GFA <- GFA[c(1,3:7)]
GFA$phenotype <- "General Fractional Anisotropy"
names(GFA) <- c("SeqId", "Phenotype 3", "n", "Beta", "SE", "P")
GFA <- GFA[match(brain$SeqId, GFA$SeqId),]
GFA <- GFA[-1]

GMD <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/00_Results_collation/imaging/gMD_result_updated_gMD_annotated_THR.csv")
GMD <- GMD[c(1,3:7)]
GMD$phenotype <- "General Mean Diffusivity"
names(GMD) <- c("SeqId", "Phenotype 4", "n", "Beta", "SE", "P")
GMD <- GMD[match(brain$SeqId, GMD$SeqId),]
GMD <- GMD[-1]

FAZ <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/00_Results_collation/imaging/Fazekas_Score_Total_result_annotated_THR.csv")
FAZ <- FAZ[c(1,3:7)]
FAZ$phenotype <- "Fazkeas White Matter Score"
names(FAZ) <- c("SeqId", "Phenotype 5", "n", "Beta", "SE", "P")
FAZ <- FAZ[match(brain$SeqId, FAZ$SeqId),]
FAZ <- FAZ[-1]

WM <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/00_Results_collation/imaging/WMH_Volume_Total_result_updated_WMHV_with_ICV_and_site_and_editor_annotated_THR.csv")
WM <- WM[c(1,3:7)]
WM$phenotype <- "White Matter Hyperintensity Volume"
names(WM) <- c("SeqId", "Phenotype 5", "n", "Beta", "SE", "P")
WM <- WM[match(brain$SeqId, WM$SeqId),]
WM <- WM[-1]

WBV <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/00_Results_collation/imaging/WBV_No_Ventricles_result_annotated_THR.csv")
WBV <- WBV[c(1,3:7)]
WBV$phenotype <- "Whole Brain Volume"
names(WBV) <- c("SeqId", "Phenotype 6", "n", "Beta", "SE", "P")
WBV <- WBV[match(brain$SeqId, WBV$SeqId),]
WBV <- WBV[-1]

# Join together MATCHED files 
join <- cbind(brain, GGM)
join <- cbind(join, GFA)
join <- cbind(join, GMD)
join <- cbind(join, FAZ)
join <- cbind(join, WM)
join <- cbind(join, WBV)

# Write out suppl table file 
write.csv(join, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/00_Results_collation/joint_imaging_suppl_table.csv", row.names = F)


######################################################################################

### COLLATE PHeWAS FILES 

screen 

R

library(tidyverse)
library(readxl)

### APOE 
apoe <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/APOE/APOE_processed/APOE_JOINT_annotated_THR.csv")
apoe$type <- "APOE"
list1 <- apoe[which(apoe$Status == "pass"),]
list1 <- list1[c("SeqId", "phenotype", "type", "Gene.Name.Name", "beta", "SE", "Pcalc", "UniProt.Full.Name")] # 16 proteins for APOE 

### IMAGING 
brain <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/00_Results_collation/imaging/brain_accel_result_annotated_THR.csv")
brain$type <- "imaging"
brain$phenotype <- "Brain Age Acceleration"
brain <- brain[which(brain$Status == "pass"),]
list1 <- rbind(list1, brain[c("SeqId", "phenotype", "type", "Gene.Name.Name", "beta", "SE", "Pcalc", "UniProt.Full.Name")])

GGM <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/00_Results_collation/imaging/Global_GM_Volume_result_annotated_THR.csv")
GGM$type <- "imaging"
GGM$phenotype <- "Global Grey Matter Volume"
GGM <- GGM[which(GGM$Status == "pass"),]
list1 <- rbind(list1, GGM[c("SeqId", "phenotype", "type", "Gene.Name.Name", "beta", "SE", "Pcalc", "UniProt.Full.Name")])

GFA <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/00_Results_collation/imaging/gFA_result_annotated_THR.csv")
GFA$type <- "imaging"
GFA$phenotype <- "General Fractional Anisotropy"
GFA <- GFA[which(GFA$Status == "pass"),]
list1 <- rbind(list1, GFA[c("SeqId", "phenotype", "type", "Gene.Name.Name", "beta", "SE", "Pcalc", "UniProt.Full.Name")])

GMD <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/00_Results_collation/imaging/gMD_result_updated_gMD_annotated_THR.csv")
GMD$type <- "imaging"
GMD$phenotype <- "General Mean Diffusivity"
GMD <- GMD[which(GMD$Status == "pass"),]
list1 <- rbind(list1, GMD[c("SeqId", "phenotype", "type", "Gene.Name.Name", "beta", "SE", "Pcalc", "UniProt.Full.Name")])

FAZ <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/00_Results_collation/imaging/Fazekas_Score_Total_result_annotated_THR.csv")
FAZ$type <- "imaging"
FAZ$phenotype <- "Fazekas White Matter Hyperintensity Score"
FAZ <- FAZ[which(FAZ$Status == "pass"),]
list1 <- rbind(list1, FAZ[c("SeqId", "phenotype", "type", "Gene.Name.Name", "beta", "SE", "Pcalc", "UniProt.Full.Name")])

WM <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/00_Results_collation/imaging/WMH_Volume_Total_result_updated_WMHV_with_ICV_and_site_and_editor_annotated_THR.csv")
WM$type <- "imaging"
WM <- WM[which(WM$Status == "pass"),]
list1 <- rbind(list1, WM[c("SeqId", "phenotype", "type", "Gene.Name.Name", "beta", "SE", "Pcalc", "UniProt.Full.Name")])

WBV <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/00_Results_collation/imaging/WBV_No_Ventricles_result_annotated_THR.csv")
WBV$type <- "imaging"
WBV$phenotype <- "Whole Brain Volume"
WBV <- WBV[which(WBV$Status == "pass"),]
list1 <- rbind(list1, WBV[c("SeqId", "phenotype", "type", "Gene.Name.Name", "beta", "SE", "Pcalc", "UniProt.Full.Name")])


### COG SCORES
g <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/00_Results_collation/cognitive/g_result_annotated_THR.csv")
g$type <- "cognitive"
g$phenotype <- "General Cognitive Ability"
g <- g[which(g$Status == "pass"),]
list1 <- rbind(list1, g[c("SeqId", "phenotype", "type", "Gene.Name.Name", "beta", "SE", "Pcalc", "UniProt.Full.Name")])

gf <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/00_Results_collation/cognitive/gf_result_annotated_THR.csv")
gf$type <- "cognitive"
gf$phenotype <- "General Fluid Cognitive Ability"
gf <- gf[which(gf$Status == "pass"),]
list1 <- rbind(list1, gf[c("SeqId", "phenotype", "type", "Gene.Name.Name", "beta", "SE", "Pcalc", "UniProt.Full.Name")])

digit_symbol <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/00_Results_collation/cognitive/digit_symbol_result_annotated_THR.csv")
digit_symbol$type <- "cognitive"
digit_symbol$phenotype <- "Processing Speed"
digit_symbol <- digit_symbol[which(digit_symbol$Status == "pass"),]
list1 <- rbind(list1, digit_symbol[c("SeqId", "phenotype", "type", "Gene.Name.Name", "beta", "SE", "Pcalc", "UniProt.Full.Name")])

LM <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/00_Results_collation/cognitive/LM_result_annotated_THR.csv")
LM$type <- "cognitive"
LM$phenotype <- "Logical Memory"
LM <- LM[which(LM$Status == "pass"),]
list1 <- rbind(list1, LM[c("SeqId", "phenotype", "type", "Gene.Name.Name", "beta", "SE", "Pcalc", "UniProt.Full.Name")])

mr_correct <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/00_Results_collation/cognitive/mr_correct_result_annotated_THR.csv")
mr_correct$type <- "cognitive"
mr_correct$phenotype <- "Non-Verbal Reasoning"
mr_correct <- mr_correct[which(mr_correct$Status == "pass"),]
list1 <- rbind(list1, mr_correct[c("SeqId", "phenotype", "type", "Gene.Name.Name", "beta", "SE", "Pcalc", "UniProt.Full.Name")])

verbal_total <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/00_Results_collation/cognitive/verbal_total_result_annotated_THR.csv")
verbal_total$type <- "cognitive"
verbal_total$phenotype <- "Verbal Reasoning"
verbal_total <- verbal_total[which(verbal_total$Status == "pass"),]
list1 <- rbind(list1, verbal_total[c("SeqId", "phenotype", "type", "Gene.Name.Name", "beta", "SE", "Pcalc", "UniProt.Full.Name")])

vocabulary <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/00_Results_collation/cognitive/vocabulary_result_annotated_THR.csv")
vocabulary$type <- "cognitive"
vocabulary$phenotype <- "Vocabulary"
vocabulary <- vocabulary[which(vocabulary$Status == "pass"),]
list1 <- rbind(list1, vocabulary[c("SeqId", "phenotype", "type", "Gene.Name.Name", "beta", "SE", "Pcalc", "UniProt.Full.Name")])

# Work out how many assocs for each phenotype and modality 
# apoe is already 11 as seen above 

dim(list1) # 405 total associations for phenos 
length(unique(list1[,4])) # 191 proteins unique to these associations

apoe <- list1 %>% filter(type == "APOE")
dim(apoe) # 14 APOE 
length(unique(apoe[,4])) # 14 unique proteins

cog <- list1 %>% filter(type == "cognitive")
dim(cog) # 296 cognitive 
length(unique(cog[,4])) # 142 unique proteins

im <- list1 %>% filter(type == "imaging")
dim(im) # 95 imaging
length(unique(im[,4])) # 60 unique proteins

# Save out list of unique proteins involved in neuro associations
list_save <- unique(list1[,4]) %>% as.data.frame()
write.csv(list_save, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/protein_list_save.csv")

#########################################

# Look at overlap between phenos 

which(apoe$SeqId %in% cog$SeqId) 

overlap <- apoe[which(apoe$SeqId %in% cog$SeqId),]

#       SeqId phenotype type Gene.Name.Name       beta         SE        Pcalc
# 4  17671-58      APOE APOE           ING4 -0.3486652 0.02930666 1.224827e-32
# 11  2797-56      APOE APOE           APOB  0.1735545 0.03084653 1.840241e-08
# 12  4337-49      APOE APOE            CRP -0.1218857 0.03101596 8.502616e-05
#                UniProt.Full.Name
# 4  Inhibitor of growth protein 4
# 11          Apolipoprotein B-100
# 12            C-reactive protein


overlap2 <- cog[which(cog$SeqId %in% apoe$SeqId),]

# > overlap2
#        SeqId                 phenotype      type Gene.Name.Name       beta
# 55   4337-49 General Cognitive Ability cognitive            CRP -0.1171013
# 209  2797-56          Processing Speed cognitive           APOB  0.1292759
# 252  4337-49      Non-Verbal Reasoning cognitive            CRP -0.1260488
# 273 17671-58      Non-Verbal Reasoning cognitive           ING4 -0.1146359
#             SE        Pcalc             UniProt.Full.Name
# 55  0.03106586 1.636059e-04            C-reactive protein
# 209 0.03479386 2.028186e-04          Apolipoprotein B-100
# 252 0.03099763 4.774474e-05            C-reactive protein
# 273 0.03131558 2.515596e-04 Inhibitor of growth protein 4


which(apoe$SeqId %in% im$SeqId) # 0 

which(cog$SeqId %in% im$SeqId) # 49

overlap3 <- cog[which(cog$SeqId %in% im$SeqId),]

overlap4 <- im[which(im$SeqId %in% cog$SeqId),]

unique(overlap4$SeqId)

unique(overlap4[,4])

# So 22 unique somamers and 22 proteins overlap between imaging and cognitive variables in FDR 

# Try to make a table to show the overlapping proteins 

overlap3 <- overlap3[-3]
names(overlap3)[2] <- "Cognitive"

overlap4 <- overlap4[-3]
names(overlap4)[2] <- "Imaging"

join <- merge(overlap3, overlap4, by = "SeqId", all = TRUE)

# Look at whether direction of effect is same or different 
join$Direction <- ifelse(join$beta.y < 0 & join$beta.x < 0 | join$beta.y > 0 & join$beta.x > 0, "Matched", "Opposite")

# Look at which is present in our EWAS with signal 
cpgs <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/cis_trans/FULL_cistrans_table_formatted_naming_second_threshold.csv")

join <- left_join(join, cpgs, by = "SeqId")
# sort out naming and order for table save out 
join2 <- join[c(1,3,13,2,8,10,11,12,4,5,6,14,15)]
join2 <- join2 %>% group_by(Gene.Name.Name.x)
join2 <- join2[order(-join2$Pcalc.y),]
join2 <- as.data.frame(join2)
names(join2) <- c("SeqId", "Gene of Protein", "UniProt Full Name", "Cognitive Phenotype", "Imaging Phenotype", "Imaging Beta", "Imaging SE", "Imaging P", "Cognitive Beta", "Cognitive SE", "Cognitive P", "Direction of Effect", "EWAS Signal")
write.csv(join2, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/common_proteins/common_proteins_V2.csv", row.names = F)

# Look at the full and see how many are in our EWAS dataset 
cpgs <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/cis_trans/FULL_cistrans_table_formatted_naming_second_threshold.csv")
cpgs <- cpgs[which(cpgs$SeqId %in% list1$SeqId),] 
dim(cpgs) # 35 pQTMs
length(unique(cpgs$Gene.of.Protein)) # 17 unique proteins
length(unique(cpgs$CpG)) # 31 unique CpGs 
length(unique(cpgs$Gene.of.CpG)) # 20 unqiue genes of CpGs 

trans <- cpgs[which(cpgs$Association.Type == "TRANS"),]
dim(trans) # 15

cis <- cpgs[which(cpgs$Association.Type == "CIS"),]
dim(cis) # 20

write.csv(cpgs, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/cis_trans/slice_neuro_final_35_pQTMs.csv", row.names = F)

# 31 unique CpGs 

# Write out a formatted table with the pQTMs that are selected for neurological markers overall 
chrom <- cpgs
chrom <- chrom %>% group_by(chrom$Gene.Name.Name)
chrom <- as.data.frame(chrom)
names(chrom) <- c("CpG", "Chromsome of CpG", "Gene of CpG", "CpG Position", "Orientation",
  "Beta", "SE", "P", "SeqId", "Gene of Protein", "UniProt", "UniProt Full Name", "Gene Start", "Gene End", "Chromosome of Gene", "Association Type")
# Add in the SNP info 
table <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Interpretation/interpretation_090621/PROTEIN_summary_with_pQTLs_linked.csv")
table <- table[c(1,4,9)]
join <- left_join(chrom, table, by = "SeqId")
write.csv(join, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/cis_trans/slice_neuro_final_35_pQTMs_pQTLs_added.csv", row.names = F)

# Isolate a list of CpG genes and protein genes to feed into enrichment analyses 
table <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/cis_trans/slice_neuro_final_35_pQTMs_pQTLs_added.csv")
list11 <- unique(table$Gene.of.CpG)
list22 <- unique(table$Gene.of.Protein)


# Now go to mQTL/eQTL script to assess lookup 


###############################################################################

### Additional step: positive and negative direction of effect table 

# Get lists of pos/neg proteins in associations - list1 has the 405 total associations combined to use 

phenos <- unique(list1$phenotype)

#  [1] "APOE"
#  [2] "Brain Age Acceleration"
#  [3] "Global Grey Matter Volume"
#  [4] "General Fractional Anisotropy"
#  [5] "General Mean Diffusivity"
#  [6] "Fazekas White Matter Hyperintensity Score"
#  [7] "WMH_Volume_Total"
#  [8] "Whole Brain Volume"
#  [9] "General Cognitive Ability"
# [10] "General Fluid Cognitive Ability"
# [11] "Processing Speed"
# [12] "Logical Memory"
# [13] "Non-Verbal Reasoning"
# [14] "Verbal Reasoning"
# [15] "Vocabulary"

table <- data.frame(Phenotype = 1:15, Positive = 1:15, Negative = 1:15)

for (i in 1:length(phenos)){
  subset_pheno <- as.character(phenos[i])
  subset <- list1[which(list1$phenotype %in% subset_pheno),]
  neg <- subset[which(subset$beta < 0),]
  pos <- subset[which(subset$beta > 0),]
  print(dim(subset)[1] == (dim(neg)[1] + dim(pos)[1])) # check that the neg and pos associations add to full amount 
  neg_list <- neg[,4]
  pos_list <- pos[,4]
  neg_list <- str_c(neg_list, collapse = ", ")
  pos_list <- str_c(pos_list, collapse = ", ")
  table[i,1] <- subset_pheno
  table[i,2] <- pos_list
  table[i,3] <- neg_list
}

write.csv(table, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/pos_neg_summary_PheWAS.csv", row.names = F)


# Write out a record of the full 405 assocs for suppl tables
write.csv(list1, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/list1_405_assocs_pheWAS.csv", row.names = F)
