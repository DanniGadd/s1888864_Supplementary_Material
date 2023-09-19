############################################################################################

### MOA 774 with genetic ORM

############################################################################################

# Run 31/01/21

# EWAS with MOA on 774 individuals, using a genetic ORM to adjust for relatedness
# Adjustement for basic covariates:
# age
# sex
# wave
# batch 
# Depression status (enriched in STRADL)
# The pQTLs have been regressed from proteins
# WBC estimates have been added as covariates to DNAm
# EXTRA: BMI and smoking added as covariates to DNAm

# This model will produce the fully adjusted EWAS results.

############################################################################################################

### Methylation file - 774 (778, excluding 4 individuals that did not have depression status identified)

############################################################################################################

cd /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/01_Basic

screen

R

## Read in methylation files and transpose so CpGs are columns and IDs are rows
library(data.table) 
w2 = readRDS("/Cluster_Filespace/Marioni_Group/Rob/Somalogic/wave2-STRADL-mvals.rds")
w3 = readRDS("/Cluster_Filespace/Marioni_Group/Rob/Somalogic/wave3-STRADL-mvals.rds")
w2 = w2[which(row.names(w2) %in% row.names(w3)),]
w3 = w3[which(row.names(w3) %in% row.names(w2)),]

## Combine w2 and w3 
meth = cbind(w2,w3)

## Add in GS ids instead of DNAm_ids 
stradl_ids <- read.table("/Cluster_Filespace/Marioni_Group/Rob/Somalogic/STRADL_DNAm_target_REM_17April2020.txt",header =T)
id.778 <- read.csv("/Cluster_Filespace/Marioni_Group/Rob/Somalogic/IDs_778.csv")

## Subset to the people in the ID file 
stradl_ids <- stradl_ids[which(stradl_ids$GS_id %in% id.778$x),]
ids = stradl_ids$DNAm_id
meth <- meth[,match(ids, colnames(meth))]
table(stradl_ids$DNAm_id == colnames(meth))
colnames(meth) <- stradl_ids$GS_id

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
phenotypes = merge(phenotypes, dep, by = "Stradl_id")

# Subset to 774 individuals with complete depression information
sample <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/774_IDs_meth_order.csv")
phenotypes <- phenotypes[which(phenotypes$Stradl_id %in% sample$Stradl_id),]

# Add WBC info into phenotypes
w2cells <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/w2_cells_joint.csv")
w3cells <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/w3_cells_joint.csv")
merged <- rbind(w2cells, w3cells)
names(merged)[1] <- "GS_id"

phenotypes = merge(phenotypes, merged, by = "GS_id")

# Join in epismoker info for those in STRADL to linker file, then to cov file 
epi <- readRDS("/Cluster_Filespace/Marioni_Group/Danni/stradl_epismoker.rds")
names(epi)[1] <- "DNAm_id"
target <- read.table("/Cluster_Filespace/Marioni_Group/Rob/Somalogic/STRADL_DNAm_target_REM_17April2020.txt",header =T)
target <- target[c(2,3)]
names(target)[2] <- "DNAm_id"
library(tidyverse)
epi <- left_join(epi, target, by = "DNAm_id")

phenotypes = merge(phenotypes, epi, by = "GS_id")

# Add BMI
BMI <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS/01_Phenotype_collation/prot_file_270121.csv")
BMI <- BMI[c("GS_id", "BMI")]

phenotypes = merge(phenotypes, BMI, by = "GS_id")

# Subset meth file to 774 individuals with complete information
meth <- meth[,which(colnames(meth) %in% phenotypes$GS_id)]

# Match ID order between phenotype and meth files 
ids= colnames(meth)
phenotypes = phenotypes[match(ids, as.character(phenotypes$GS_id)),]
table(colnames(meth) == phenotypes$GS_id)

# Transpose methylation data
meth = t(meth)

# Resdualise meth matrix 
for(i in 1:(ncol(meth))){
  print(i)
  meth[,i] <- resid(lm(meth[,i] ~ phenotypes$age_stradl + as.factor(phenotypes$sex) + as.factor(phenotypes$Batch) + as.factor(phenotypes$wave) + as.factor(phenotypes$combined) + phenotypes$Bcell + phenotypes$CD4T + phenotypes$CD8T + phenotypes$Gran + phenotypes$NK + phenotypes$Mono + phenotypes$Eos + phenotypes$smokingScore + phenotypes$BMI, na.action = na.exclude))
}

meth <- as.data.frame(meth)

## Combined IDS of individuals (from the relatedness subtyping script for GRM)
osca_dat <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/774_IDs_meth_order.csv")

# match meth file to the order in the FID IID 774 file 
ids = osca_dat$GS_id
meth <- meth[match(ids, row.names(meth)),]
table(osca_dat$GS_id == row.names(meth))

## Combine IDs for FID and IID and rearrange columns so it goes FID, IID, cpg1, cpg2... 
tmp <- meth
tmp$IID <- osca_dat$IID
tmp$FID <- osca_dat$FID
tmp <- tmp[,c(ncol(tmp), ncol(tmp)-1, 1:(ncol(tmp)-2))]
tmp$FID <- as.character(tmp$FID)
tmp$IID <- as.character(tmp$IID)
identical(tmp$IID, as.character(osca_dat$IID)) # true 

# This file has everything regressed for MOA analyses
fwrite(tmp, file="/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Methylation_preps/methylation_774_meth_order_FULL_310121.txt", row.names=F, quote = F, sep=' ')


############################################################################################################

### Phenotype file - 778

############################################################################################################

cd /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/01_Basic

screen

R

# Read in order of meth files we are matching to - combined IDS of individuals (from the relatedness subtyping script for GRM)
order <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/774_IDs_meth_order.csv")

# Read in protein file which has GS id and stradl id, then all 4,235 proteins that have been preprocessed
pheno2 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Protein_preps/Phenotype_file_774_pQTLs_regressed.csv", check.names = F)
names(pheno2)[1] <- "ID"
pheno2 <- pheno2[-2]

# Check to make sure the id in the proteins file is the same as the id in the ID reference file 
identical(pheno2$ID, order$GS_id) # TRUE

# Assign FID and IID from the ID reference file 
pheno2$IID <- order$IID
pheno2$FID <- order$FID

# Order so proteins are first 
pheno2 <- pheno2[c(2:4236,1,4237,4238)]

identical(pheno2$IID, order$IID) # TRUE

# Write out phenotype files for a heritability run in the 4235
location <- "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/04_FULL/batches/"

## Write out protein files so FID, IID, phenotype 
for(i in 1:350){ 
name <- names(pheno2)[i]
write.table(pheno2[,c("FID","IID",name)], file = paste0(location, "batch1/",name,".phen"), col.names = F, row.names=F, quote = F, sep=' ')
}

## Write out protein files so FID, IID, phenotype 
for(i in 351:700){ 
name <- names(pheno2)[i]
write.table(pheno2[,c("FID","IID",name)], file = paste0(location, "batch2/",name,".phen"), col.names = F, row.names=F, quote = F, sep=' ')
}

## Write out protein files so FID, IID, phenotype 
for(i in 701:1050){ 
name <- names(pheno2)[i]
write.table(pheno2[,c("FID","IID",name)], file = paste0(location, "batch3/",name,".phen"), col.names = F, row.names=F, quote = F, sep=' ')
}

## Write out protein files so FID, IID, phenotype 
for(i in 1051:1400){ 
name <- names(pheno2)[i]
write.table(pheno2[,c("FID","IID",name)], file = paste0(location, "batch4/",name,".phen"), col.names = F, row.names=F, quote = F, sep=' ')
}

## Write out protein files so FID, IID, phenotype 
for(i in 1401:1750){ 
name <- names(pheno2)[i]
write.table(pheno2[,c("FID","IID",name)], file = paste0(location, "batch5/",name,".phen"), col.names = F, row.names=F, quote = F, sep=' ')
}

## Write out protein files so FID, IID, phenotype 
for(i in 1751:2100){ 
name <- names(pheno2)[i]
write.table(pheno2[,c("FID","IID",name)], file = paste0(location, "batch6/",name,".phen"), col.names = F, row.names=F, quote = F, sep=' ')
}

## Write out protein files so FID, IID, phenotype 
for(i in 2101:2450){ 
name <- names(pheno2)[i]
write.table(pheno2[,c("FID","IID",name)], file = paste0(location, "batch7/",name,".phen"), col.names = F, row.names=F, quote = F,sep=' ')
}

## Write out protein files so FID, IID, phenotype 
for(i in 2451:2800){ 
name <- names(pheno2)[i]
write.table(pheno2[,c("FID","IID",name)], file = paste0(location, "batch8/",name,".phen"), col.names = F, row.names=F, quote = F, sep=' ')
}

## Write out protein files so FID, IID, phenotype 
for(i in 2801:3150){ 
name <- names(pheno2)[i]
write.table(pheno2[,c("FID","IID",name)], file = paste0(location, "batch9/",name,".phen"), col.names = F, row.names=F, quote = F,sep=' ')
}

## Write out protein files so FID, IID, phenotype 
for(i in 3151:3500){ 
name <- names(pheno2)[i]
write.table(pheno2[,c("FID","IID",name)], file = paste0(location, "batch10/",name,".phen"), col.names = F, row.names=F, quote = F, sep=' ')
}

## Write out protein files so FID, IID, phenotype 
for(i in 3501:3850){ 
name <- names(pheno2)[i]
write.table(pheno2[,c("FID","IID",name)], file = paste0(location, "batch11/",name,".phen"), col.names = F, row.names=F, quote = F,sep=' ')
}

## Write out protein files so FID, IID, phenotype 
for(i in 3851:4235){ 
name <- names(pheno2)[i]
write.table(pheno2[,c("FID","IID",name)], file = paste0(location, "batch12/",name,".phen"), col.names = F, row.names=F, quote = F,sep=' ')
}

# The phenotype (protein) files are prepped per protein in the 12 batch folders to run in parallel



############################################################################################################

### MAKE BINARY FILES AND ORM - 774

############################################################################################################

### This doesnt change as the meth file is the same (only phenotype has had eGFR not regressed from it so meth files can be used again)
cd /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/01_Basic

# Make Binary Methylation File
osca_Linux --efile /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Methylation_preps/methylation_774_meth_order_FULL_310121.txt --methylation-m --make-bod --out /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Methylation_preps/methylation_774_meth_order_FULL_310121

# Regenerate the opi file with correct annotations 
anno = readRDS("/Cluster_Filespace/Marioni_Group/Daniel/EPIC_AnnotationObject_df.rds")
opi1 = read.table("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Methylation_preps/methylation_774_meth_order_FULL_310121.opi", header=F, stringsAsFactors=F)
opi <- anno[opi1$V2,c("chr","Name", "pos","UCSC_RefGene_Name", "strand")]
opi$chr <- gsub("chr", "", opi$chr)
opi$chr <- gsub("X", "23", opi$chr)
opi$chr <- gsub("Y", "24", opi$chr)
opi$chr <- as.numeric(opi$chr)
opi$UCSC_RefGene_Name  <- as.factor(opi$UCSC_RefGene_Name )
opi$strand <- as.factor(opi$strand)
opi[which(opi$UCSC_RefGene_Name ==""), "UCSC_RefGene_Name"] <- NA
write.table(opi, file="/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Methylation_preps/methylation_774_meth_order_FULL_310121.opi",  # make sure to keep the same filename as before
                col.names=F, 
                row.names=F, 
                quote=F, sep='\t')

# check that the annotations are present after the update
head /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Methylation_preps/methylation_774_meth_order_FULL_310121.opi


##########################################################################################

# BATCH EWAS 

cd /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/04_FULL/batches/batch1/
for i in *.phen 

do

  A=$( echo $i | cut -d"/" -f3)
  B=$( echo $A | cut -d. -f1)

osca_Linux --moa --orm /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/GRM/ormtest774 --befile /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Methylation_preps/methylation_774_meth_order_FULL_310121 --pheno $i --out ${B}_MOA_774_with_GRM --methylation-m


done 


cd /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/04_FULL/batches/batch2/
for i in *.phen 

do

  A=$( echo $i | cut -d"/" -f3)
  B=$( echo $A | cut -d. -f1)

osca_Linux --moa --orm /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/GRM/ormtest774 --befile /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Methylation_preps/methylation_774_meth_order_FULL_310121 --pheno $i --out ${B}_MOA_774_with_GRM --methylation-m


done 


cd /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/04_FULL/batches/batch3/
for i in *.phen 

do

  A=$( echo $i | cut -d"/" -f3)
  B=$( echo $A | cut -d. -f1)

osca_Linux --moa --orm /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/GRM/ormtest774 --befile /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Methylation_preps/methylation_774_meth_order_FULL_310121 --pheno $i --out ${B}_MOA_774_with_GRM --methylation-m


done 


cd /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/04_FULL/batches/batch4/
for i in *.phen 

do

  A=$( echo $i | cut -d"/" -f3)
  B=$( echo $A | cut -d. -f1)

osca_Linux --moa --orm /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/GRM/ormtest774 --befile /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Methylation_preps/methylation_774_meth_order_FULL_310121 --pheno $i --out ${B}_MOA_774_with_GRM --methylation-m


done 


cd /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/04_FULL/batches/batch5/
for i in *.phen 

do

  A=$( echo $i | cut -d"/" -f3)
  B=$( echo $A | cut -d. -f1)

osca_Linux --moa --orm /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/GRM/ormtest774 --befile /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Methylation_preps/methylation_774_meth_order_FULL_310121 --pheno $i --out ${B}_MOA_774_with_GRM --methylation-m


done 


cd /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/04_FULL/batches/batch6/
for i in *.phen 

do

  A=$( echo $i | cut -d"/" -f3)
  B=$( echo $A | cut -d. -f1)

osca_Linux --moa --orm /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/GRM/ormtest774 --befile /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Methylation_preps/methylation_774_meth_order_FULL_310121 --pheno $i --out ${B}_MOA_774_with_GRM --methylation-m


done 


cd /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/04_FULL/batches/batch7/
for i in *.phen 

do

  A=$( echo $i | cut -d"/" -f3)
  B=$( echo $A | cut -d. -f1)

osca_Linux --moa --orm /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/GRM/ormtest774 --befile /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Methylation_preps/methylation_774_meth_order_FULL_310121 --pheno $i --out ${B}_MOA_774_with_GRM --methylation-m


done 


cd /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/04_FULL/batches/batch8/
for i in *.phen 

do

  A=$( echo $i | cut -d"/" -f3)
  B=$( echo $A | cut -d. -f1)

osca_Linux --moa --orm /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/GRM/ormtest774 --befile /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Methylation_preps/methylation_774_meth_order_FULL_310121 --pheno $i --out ${B}_MOA_774_with_GRM --methylation-m


done 


cd /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/04_FULL/batches/batch9/
for i in *.phen 

do

  A=$( echo $i | cut -d"/" -f3)
  B=$( echo $A | cut -d. -f1)

osca_Linux --moa --orm /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/GRM/ormtest774 --befile /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Methylation_preps/methylation_774_meth_order_FULL_310121 --pheno $i --out ${B}_MOA_774_with_GRM --methylation-m


done 


cd /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/04_FULL/batches/batch10/
for i in *.phen 

do

  A=$( echo $i | cut -d"/" -f3)
  B=$( echo $A | cut -d. -f1)

osca_Linux --moa --orm /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/GRM/ormtest774 --befile /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Methylation_preps/methylation_774_meth_order_FULL_310121 --pheno $i --out ${B}_MOA_774_with_GRM --methylation-m


done 


cd /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/04_FULL/batches/batch11/
for i in *.phen 

do

  A=$( echo $i | cut -d"/" -f3)
  B=$( echo $A | cut -d. -f1)

osca_Linux --moa --orm /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/GRM/ormtest774 --befile /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Methylation_preps/methylation_774_meth_order_FULL_310121 --pheno $i --out ${B}_MOA_774_with_GRM --methylation-m


done 


cd /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/04_FULL/batches/batch12/
for i in *.phen 

do

  A=$( echo $i | cut -d"/" -f3)
  B=$( echo $A | cut -d. -f1)

osca_Linux --moa --orm /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/GRM/ormtest774 --befile /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Methylation_preps/methylation_774_meth_order_FULL_310121 --pheno $i --out ${B}_MOA_774_with_GRM --methylation-m


done 



####################################################################################################

### PROCESSING THE EWAS RESULTS FOR BATCHES ABOVE - MOA 778 with genetic ORM - AND covs - AND pQTLs

####################################################################################################

# Processing the 170521 results 

# Split this up for speed into 12 screens 

# Go through each protein and subset to the significant sites 

screen

R

library(data.table)

list <- c("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/04_FULL/batches/batch1/")

sig_hits <- "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/04_FULL/batches/sig_hits/"

batch <- c("b1")

for(i in 1:length(list)){

  path <- list[i]
  batch <- batch[i]
  setwd(paste0(path))

  L <- list.files(".", ".mlma")

  # results tables for each protein (all vs sig)
  for(j in 1:length(L)){
    file <- fread(L[j], header = T)
    file <- as.data.frame(file)
    name <- as.character(L[j])
    name <- gsub("_.*", "", name)
    file$SeqId <- name
    # write.csv(file, paste0(path, name, "_results_table.csv"), row.names = F)
    sub <- file[-which(file$p > 0.000000036),]
    # write.csv(sub, paste0(path, name, "_results_sig.csv"), row.names = F)
    write.csv(sub, paste0(sig_hits, batch, "_", name, "_results_sig.csv"), row.names = F)
  }
}


screen

R

library(data.table)

list <- c("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/04_FULL/batches/batch2/")

sig_hits <- "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/04_FULL/batches/sig_hits/"

batch <- c("b2")

for(i in 1:length(list)){

  path <- list[i]
  batch <- batch[i]
  setwd(paste0(path))

  L <- list.files(".", ".mlma")

  # results tables for each protein (all vs sig)
  for(j in 1:length(L)){
    file <- fread(L[j], header = T)
    file <- as.data.frame(file)
    name <- as.character(L[j])
    name <- gsub("_.*", "", name)
    file$SeqId <- name
    # write.csv(file, paste0(path, name, "_results_table.csv"), row.names = F)
    sub <- file[-which(file$p > 0.000000036),]
    # write.csv(sub, paste0(path, name, "_results_sig.csv"), row.names = F)
    write.csv(sub, paste0(sig_hits, batch, "_", name, "_results_sig.csv"), row.names = F)
  }
}


screen

R

library(data.table)

list <- c("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/04_FULL/batches/batch3/")

sig_hits <- "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/04_FULL/batches/sig_hits/"

batch <- c("b3")

for(i in 1:length(list)){

  path <- list[i]
  batch <- batch[i]
  setwd(paste0(path))

  L <- list.files(".", ".mlma")

  # results tables for each protein (all vs sig)
  for(j in 1:length(L)){
    file <- fread(L[j], header = T)
    file <- as.data.frame(file)
    name <- as.character(L[j])
    name <- gsub("_.*", "", name)
    file$SeqId <- name
    # write.csv(file, paste0(path, name, "_results_table.csv"), row.names = F)
    sub <- file[-which(file$p > 0.000000036),]
    # write.csv(sub, paste0(path, name, "_results_sig.csv"), row.names = F)
    write.csv(sub, paste0(sig_hits, batch, "_", name, "_results_sig.csv"), row.names = F)
  }
}



screen

R

library(data.table)

list <- c("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/04_FULL/batches/batch4/")

sig_hits <- "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/04_FULL/batches/sig_hits/"

batch <- c("b4")

for(i in 1:length(list)){

  path <- list[i]
  batch <- batch[i]
  setwd(paste0(path))

  L <- list.files(".", ".mlma")

  # results tables for each protein (all vs sig)
  for(j in 1:length(L)){
    file <- fread(L[j], header = T)
    file <- as.data.frame(file)
    name <- as.character(L[j])
    name <- gsub("_.*", "", name)
    file$SeqId <- name
    # write.csv(file, paste0(path, name, "_results_table.csv"), row.names = F)
    sub <- file[-which(file$p > 0.000000036),]
    # write.csv(sub, paste0(path, name, "_results_sig.csv"), row.names = F)
    write.csv(sub, paste0(sig_hits, batch, "_", name, "_results_sig.csv"), row.names = F)
  }
}


screen

R

library(data.table)

list <- c("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/04_FULL/batches/batch5/")

sig_hits <- "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/04_FULL/batches/sig_hits/"

batch <- c("b5")

for(i in 1:length(list)){

  path <- list[i]
  batch <- batch[i]
  setwd(paste0(path))

  L <- list.files(".", ".mlma")

  # results tables for each protein (all vs sig)
  for(j in 1:length(L)){
    file <- fread(L[j], header = T)
    file <- as.data.frame(file)
    name <- as.character(L[j])
    name <- gsub("_.*", "", name)
    file$SeqId <- name
    # write.csv(file, paste0(path, name, "_results_table.csv"), row.names = F)
    sub <- file[-which(file$p > 0.000000036),]
    # write.csv(sub, paste0(path, name, "_results_sig.csv"), row.names = F)
    write.csv(sub, paste0(sig_hits, batch, "_", name, "_results_sig.csv"), row.names = F)
  }
}


screen

R

library(data.table)

list <- c("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/04_FULL/batches/batch6/")

sig_hits <- "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/04_FULL/batches/sig_hits/"

batch <- c("b6")

for(i in 1:length(list)){

  path <- list[i]
  batch <- batch[i]
  setwd(paste0(path))

  L <- list.files(".", ".mlma")

  # results tables for each protein (all vs sig)
  for(j in 1:length(L)){
    file <- fread(L[j], header = T)
    file <- as.data.frame(file)
    name <- as.character(L[j])
    name <- gsub("_.*", "", name)
    file$SeqId <- name
    # write.csv(file, paste0(path, name, "_results_table.csv"), row.names = F)
    sub <- file[-which(file$p > 0.000000036),]
    # write.csv(sub, paste0(path, name, "_results_sig.csv"), row.names = F)
    write.csv(sub, paste0(sig_hits, batch, "_", name, "_results_sig.csv"), row.names = F)
  }
}


screen

R

library(data.table)

list <- c("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/04_FULL/batches/batch7/")

sig_hits <- "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/04_FULL/batches/sig_hits/"

batch <- c("b7")

for(i in 1:length(list)){

  path <- list[i]
  batch <- batch[i]
  setwd(paste0(path))

  L <- list.files(".", ".mlma")

  # results tables for each protein (all vs sig)
  for(j in 1:length(L)){
    file <- fread(L[j], header = T)
    file <- as.data.frame(file)
    name <- as.character(L[j])
    name <- gsub("_.*", "", name)
    file$SeqId <- name
    # write.csv(file, paste0(path, name, "_results_table.csv"), row.names = F)
    sub <- file[-which(file$p > 0.000000036),]
    # write.csv(sub, paste0(path, name, "_results_sig.csv"), row.names = F)
    write.csv(sub, paste0(sig_hits, batch, "_", name, "_results_sig.csv"), row.names = F)
  }
}


screen

R

library(data.table)

list <- c("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/04_FULL/batches/batch8/")

sig_hits <- "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/04_FULL/batches/sig_hits/"

batch <- c("b8")

for(i in 1:length(list)){

  path <- list[i]
  batch <- batch[i]
  setwd(paste0(path))

  L <- list.files(".", ".mlma")

  # results tables for each protein (all vs sig)
  for(j in 1:length(L)){
    file <- fread(L[j], header = T)
    file <- as.data.frame(file)
    name <- as.character(L[j])
    name <- gsub("_.*", "", name)
    file$SeqId <- name
    # write.csv(file, paste0(path, name, "_results_table.csv"), row.names = F)
    sub <- file[-which(file$p > 0.000000036),]
    # write.csv(sub, paste0(path, name, "_results_sig.csv"), row.names = F)
    write.csv(sub, paste0(sig_hits, batch, "_", name, "_results_sig.csv"), row.names = F)
  }
}



screen

R

library(data.table)

list <- c("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/04_FULL/batches/batch9/")

sig_hits <- "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/04_FULL/batches/sig_hits/"

batch <- c("b9")

for(i in 1:length(list)){

  path <- list[i]
  batch <- batch[i]
  setwd(paste0(path))

  L <- list.files(".", ".mlma")

  # results tables for each protein (all vs sig)
  for(j in 1:length(L)){
    file <- fread(L[j], header = T)
    file <- as.data.frame(file)
    name <- as.character(L[j])
    name <- gsub("_.*", "", name)
    file$SeqId <- name
    # write.csv(file, paste0(path, name, "_results_table.csv"), row.names = F)
    sub <- file[-which(file$p > 0.000000036),]
    # write.csv(sub, paste0(path, name, "_results_sig.csv"), row.names = F)
    write.csv(sub, paste0(sig_hits, batch, "_", name, "_results_sig.csv"), row.names = F)
  }
}



screen

R

library(data.table)

list <- c("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/04_FULL/batches/batch10/")

sig_hits <- "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/04_FULL/batches/sig_hits/"

batch <- c("b10")

for(i in 1:length(list)){

  path <- list[i]
  batch <- batch[i]
  setwd(paste0(path))

  L <- list.files(".", ".mlma")

  # results tables for each protein (all vs sig)
  for(j in 1:length(L)){
    file <- fread(L[j], header = T)
    file <- as.data.frame(file)
    name <- as.character(L[j])
    name <- gsub("_.*", "", name)
    file$SeqId <- name
    # write.csv(file, paste0(path, name, "_results_table.csv"), row.names = F)
    sub <- file[-which(file$p > 0.000000036),]
    # write.csv(sub, paste0(path, name, "_results_sig.csv"), row.names = F)
    write.csv(sub, paste0(sig_hits, batch, "_", name, "_results_sig.csv"), row.names = F)
  }
}


screen

R

library(data.table)

list <- c("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/04_FULL/batches/batch11/")

sig_hits <- "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/04_FULL/batches/sig_hits/"

batch <- c("b11")

for(i in 1:length(list)){

  path <- list[i]
  batch <- batch[i]
  setwd(paste0(path))

  L <- list.files(".", ".mlma")

  # results tables for each protein (all vs sig)
  for(j in 1:length(L)){
    file <- fread(L[j], header = T)
    file <- as.data.frame(file)
    name <- as.character(L[j])
    name <- gsub("_.*", "", name)
    file$SeqId <- name
    # write.csv(file, paste0(path, name, "_results_table.csv"), row.names = F)
    sub <- file[-which(file$p > 0.000000036),]
    # write.csv(sub, paste0(path, name, "_results_sig.csv"), row.names = F)
    write.csv(sub, paste0(sig_hits, batch, "_", name, "_results_sig.csv"), row.names = F)
  }
}


screen

R

library(data.table)

list <- c("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/04_FULL/batches/batch12/")

sig_hits <- "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/04_FULL/batches/sig_hits/"

batch <- c("b12")

for(i in 1:length(list)){

  path <- list[i]
  batch <- batch[i]
  setwd(paste0(path))

  L <- list.files(".", ".mlma")

  # results tables for each protein (all vs sig)
  for(j in 1:length(L)){
    file <- fread(L[j], header = T)
    file <- as.data.frame(file)
    name <- as.character(L[j])
    name <- gsub("_.*", "", name)
    file$SeqId <- name
    # write.csv(file, paste0(path, name, "_results_table.csv"), row.names = F)
    sub <- file[-which(file$p > 0.000000036),]
    # write.csv(sub, paste0(path, name, "_results_sig.csv"), row.names = F)
    write.csv(sub, paste0(sig_hits, batch, "_", name, "_results_sig.csv"), row.names = F)
  }
}

##########################################################################

# Collate top hits for proteins 

screen

R

path <- "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/04_FULL/batches/sig_hits/"

setwd(path)

L <- list.files(".", ".csv")
L # 4231 converged 

files <- lapply(L, read.csv)
names <- as.character(L)
batch <- gsub("_.*", "", names)
marker <- gsub("_results.*", "", names)
marker <- gsub(".*_", "", marker)
names(files) <- marker
osca <- do.call(rbind, files)
osca <- osca[c(9,1:8)]

# annotations file add in
library(tidyverse) 
library(readxl)
anno <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Annotations_protein_file/annotation_formatted_for_paper.csv")
comb <- left_join(osca, anno, by = "SeqId")

not <- anno[which(!anno$SeqId %in% marker),]


# non convergence models 
#        SeqId UniProt Gene.Name.Name
# 223  15509-2  P54802          NAGLU
# 899  15584-9  P36980          CFHR2
# 1702 4407-10  P26927           MST1
# 2670  6402-8  Q9UKJ1          PILRA
#                                     UniProt.Full.Name
# 223                     Alpha-N-acetylglucosaminidase
# 899             Complement factor H-related protein 2
# 1702            Hepatocyte growth factor-like protein
# 2670 Paired immunoglobulin-like type 2 receptor alpha


# Write out file at ewas threshold of significance 
write.csv(comb, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Results_collation/FULL_EWAS_results_3.6_threshold_020221.csv", row.names = F)

# Filter to CpG-protein correction value of adjusted p 0.05/143/772619 = 4.525521e-10
# 0.0000000004525521

comb2_filt2 <- comb[which(comb$p < 0.0000000004525521),]

dim(comb2_filt2) # 2928  - at cpg/protein adjusted significance 
comb2 <- comb2_filt2[order(comb2_filt2$p),]
write.csv(comb2, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Results_collation/FULL_EWAS_results_second_threshold_020221.csv", row.names = F)


# Filter to CpG-protein correction value of adjusted p 0.05/4235/772619 = 1.528098e-11
# > 0.05 / 4235 / 772619
# [1] 1.528098e-11

# 0.00000000001528098

comb2_filt2 <- comb[which(comb$p < 0.00000000001528098),]

dim(comb2_filt2) # 2080  - at cpg/protein adjusted significance 
comb2 <- comb2_filt2[order(comb2_filt2$p),]
write.csv(comb2, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Results_collation/FULL_EWAS_results_020221.csv", row.names = F)

########################################################################################################

### Subset to nested model structure across models and format tables for suppl reporting 

########################################################################################################

# FIRST FORMAT TABLES 

cd /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Results_collation/Formatted_tables/

screen 

R

library(readxl)
library(tidyverse)


# Load corrected basic model 
base <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Results_collation/Basic_EWAS_results_second_threshold_020221.csv")
base <- base[c(3,2,5,4,6:9,1,10,11,12)]
names(base) <- c("CpG", "Chromsome of CpG", "Gene of CpG", "CpG position", "Orientation",
  "Beta", "SE", "P", "SeqId", "UniProt", "Gene of protein", "UniProt Full Name") 
base <- base[order(base$P),]
write.csv(base, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Results_collation/Formatted_tables/Basic_EWAS_results_020221_formatted_second_threshold.csv", row.names = F)

# > dim(base)
# [1] 238245     12

# Load WBC SNP corrected model 
WBC <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Results_collation/WBC_EWAS_results_second_threshold_020221.csv")
WBC <- WBC[c(3,2,5,4,6:9,1,10,11,12)]
names(WBC) <- c("CpG", "Chromsome of CpG", "Gene of CpG", "CpG position", "Orientation",
  "Beta", "SE", "P", "SeqId", "UniProt", "Gene of protein", "UniProt Full Name")
WBC <- WBC[order(WBC$P),]
write.csv(WBC, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Results_collation/Formatted_tables/WBC_EWAS_results_020221_formatted_second_threshold.csv", row.names = F)

# > dim(WBC)
# [1] 3213   12


# Load full corrected model 
full <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Results_collation/FULL_EWAS_results_second_threshold_020221.csv")
full <- full[c(3,2,5,4,6:9,1,10,11,12)]
names(full) <-  c("CpG", "Chromsome of CpG", "Gene of CpG", "CpG position", "Orientation",
  "Beta", "SE", "P", "SeqId", "UniProt", "Gene of protein", "UniProt Full Name") # 2797 remain significant 
full <- full[order(full$P),]
write.csv(full, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Results_collation/Formatted_tables/FULL_EWAS_results_020221_formatted_second_threshold.csv", row.names = F)

# > dim(full)
# [1] 2928   12

# NEXT LOOK AT WHICH ARE CONSISTENT ACROSS MODELS 

# From base to middle model 
base <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Results_collation/Basic_EWAS_results_second_threshold_020221.csv")
WBC <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Results_collation/WBC_EWAS_results_second_threshold_020221.csv")
keep = paste(base$SeqId, base$Probe, sep = "_") # 3213
WBC$retain <- paste(WBC$SeqId, WBC$Probe, sep = "_")
WBC <- WBC[which(WBC$retain %in% keep),] # 2718 
WBC$retain <- NULL
write.csv(WBC, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Results_collation/Formatted_tables/retained_basic_to_middle_2718_of_3213_unformatted_second_threshold.csv", row.names = F)

# From middle model (filtered) to the final full model 
WBC <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Results_collation/WBC_EWAS_results_second_threshold_020221.csv")
full <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Results_collation/FULL_EWAS_results_second_threshold_020221.csv")
keep = paste(WBC$SeqId, WBC$Probe, sep = "_")
full$retain <- paste(full$SeqId, full$Probe, sep = "_")
full <- full[which(full$retain %in% keep),] 
full$retain <- NULL
full <- full[c(3,2,5,4,6:9,1,10,11,12)]
names(full) <- c("CpG", "Chromsome of CpG", "Gene of CpG", "CpG position", "Orientation", # 2847
  "Beta", "SE", "P", "SeqId", "Gene of protein", "UniProt", "UniProt Full Name") 
write.csv(full, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Results_collation/Formatted_tables/retained_middle_to_full_covs_2847_of_2928.csv", row.names = F)

# Look at how many we lose from the final model overall 
pqtl <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Results_collation/Formatted_tables/retained_basic_to_middle_2718_of_3213_unformatted_second_threshold.csv")
filt <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Results_collation/Formatted_tables/retained_middle_to_full_covs_2847_of_2928.csv")

pqtl$retain = paste(pqtl$SeqId, pqtl$Probe, sep = "_")
filt$retain <- paste(filt$SeqId, filt$CpG, sep = "_")

diff <- pqtl[which(!pqtl$retain %in% filt$retain),]
names(diff)[10] <- "gene"
unique(diff$gene)

diff <- pqtl[which(pqtl$retain %in% filt$retain),] # 2847


# How many proteins were involved in significant full model results

full <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Results_collation/Formatted_tables/FULL_EWAS_results_020221_formatted_second_threshold.csv")
length(unique(full[,11])) # 191 unique genes
length(unique(full[,9])) # 195 unique somamers 

# Subset so that we can identify which had multiple somamers
sub <- full[c(11,9)]
sub$match <- paste0(sub[,1], "_", sub[,2])
unique(sub$match) # 146
sub <- unique(sub)
sub <- sub[order(sub[,1]), ]

##############################################################################

### NOW LOOK AT WHETHER WE REPLICATE ZAGHLOOL ET AL 98 pQTM RESULTS VS FULL 

##############################################################################

# Read in zaghlool et al results file with 98 pQTMs and key stats 
library(readxl)
library(tidyverse)
zag <- read_excel("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_240521_results/formatted_tables/Zaghlool_41467_2019_13831_MOESM5_ESM_EDITED.xlsx")
zag <- as.data.frame(zag)
zag$proteinID <- sub("\\_.", "", zag$proteinID) # get seqIds for matching 
list <- zag$proteinID
zag$sig <- paste(zag$proteinID, zag$cpg, sep = "_")

# Work out how many comparable proteins there are accross the 2 studies that could in theory be replicated (i.e. which of the 98 pQTMs did we test for in our EWAS)
stradl <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/prep/Phenotype_file_778_eGFR_inlcuded_pQTLs_regressed.csv", check.names = F)
length(which(zag$proteinID %in% colnames(stradl))) # 94 matched 

zag_not <- zag[which(!zag$proteinID %in% colnames(stradl)),]

# 4 which were not comparable across the 2 studies (involving 3 proteins)
#    order   label proteinID        cpg            p       beta         se   N
# 29    29   NLRC5   5099-14 cg07839457 1.288322e-17 -0.2742358 0.03146097 935
# 70    70   NLRC5   5099-14 cg08159663 9.221454e-13 -0.2458608 0.03386822 800
# 83    83   PRTN3   3514-49 cg27535410 1.205885e-11 -0.2187312 0.03186051 938
# 94    94 SIGLEC5    5125-6 cg09488502 4.893690e-11  0.2202263 0.03307355 869
#                   sig
# 29 5099-14_cg07839457
# 70 5099-14_cg08159663
# 83 3514-49_cg27535410
# 94  5125-6_cg09488502

# make sure only 94 comparable included 
zag <- zag[which(zag$proteinID %in% colnames(stradl)),]

### make sure all cpgs in the zaghlool et al list are present n our meth dataset used 

# first read the meth data in as above 
## Read in methylation file used for BayesR and transpose so CpGs are columns and IDs are rows
library(data.table) 
w2 = readRDS("/Cluster_Filespace/Marioni_Group/Rob/Somalogic/wave2-STRADL-mvals.rds")
w3 = readRDS("/Cluster_Filespace/Marioni_Group/Rob/Somalogic/wave3-STRADL-mvals.rds")
w2 = w2[which(row.names(w2) %in% row.names(w3)),]
w3 = w3[which(row.names(w3) %in% row.names(w2)),]

## Combine w2 and w3 
meth = cbind(w2,w3)

## Add in GS ids instead of DNAm_ids 
stradl_ids <- read.table("/Cluster_Filespace/Marioni_Group/Rob/Somalogic/STRADL_DNAm_target_REM_17April2020.txt",header =T)
id.778 <- read.csv("/Cluster_Filespace/Marioni_Group/Rob/Somalogic/IDs_778.csv")

## Subset to the people in the ID file 
stradl_ids <- stradl_ids[which(stradl_ids$GS_id %in% id.778$x),]
ids = stradl_ids$DNAm_id
meth <- meth[,match(ids, colnames(meth))]
table(stradl_ids$DNAm_id == colnames(meth))
colnames(meth) <- stradl_ids$GS_id

# next look for common cpgs across both studies for the 98 pQTMs
meth <- meth 
dim(meth) #  772619    778
length(which(zag$cpg %in% rownames(meth))) # only 81 CpGs are common to both studies 

# lets look at which are not included 

not_cpg <- zag[which(!zag$cpg %in% rownames(meth)),]


# > not_cpg
#    order label proteinID        cpg            p       beta         se   N
# 7      7 ICAM5   5124-69 cg10604476 6.092957e-25  0.3556380 0.03334739 812
# 12    12 PAPPA   4148-49 cg07562835 6.053805e-24 -0.3245133 0.03126409 920
# 19    19 PAPPA   4148-49 cg22805603 2.312486e-20  0.2952428 0.03119564 938
# 26    26 PAPPA   4148-49 cg18463686 3.845588e-18  0.2802216 0.03160582 923
# 38    38 PAPPA   4148-49 cg07708453 3.103772e-16  0.2621269 0.03150946 938
# 45    45 PAPPA   4148-49 cg03149567 2.478359e-15 -0.2553335 0.03170958 929
# 49    49 PAPPA   4148-49 cg06756385 1.092719e-14 -0.2484275 0.03162756 938
# 50    50 NLRC5   3311-27 cg08159663 1.441578e-14 -0.2898182 0.03682426 664
# 69    69 NLRC5   3485-28 cg00218406 3.968618e-13 -0.2340531 0.03179146 934
# 77    77 NLRC5   3292-75 cg08159663 3.382490e-12 -0.2394780 0.03387316 800
# 78    78 PAPPA   4148-49 cg11456013 3.891687e-12 -0.2261974 0.03215273 918
# 85    85    C4   4481-34 cg13028630 1.375174e-11 -0.2457877 0.03577153 724
# 96    96 PAPPA   4148-49 cg26922697 7.399367e-11 -0.2103716 0.03193007 937
#                   sig
# 7  5124-69_cg10604476
# 12 4148-49_cg07562835
# 19 4148-49_cg22805603
# 26 4148-49_cg18463686
# 38 4148-49_cg07708453
# 45 4148-49_cg03149567
# 49 4148-49_cg06756385
# 50 3311-27_cg08159663
# 69 3485-28_cg00218406
# 77 3292-75_cg08159663
# 78 4148-49_cg11456013
# 85 4481-34_cg13028630
# 96 4148-49_cg26922697

# We need to therefore subset to just those associations which are comaprable (i.e proteins and CpGs)
# 94 pQTMs currently in zag file - now remove CpGs that arent comparable 
zag <- zag[which(zag$cpg %in% rownames(meth)),]
dim(zag) # 81  9
# so a possible 81 pQTMs are comparable across our studies 


# Read in results files 
res <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Results_collation/Formatted_tables/FULL_EWAS_results_020221_formatted_second_threshold.csv")
res <- res[which(res$SeqId %in% list),]
res$index <- paste(res$SeqId, res$CpG, sep = "_")

length(which(zag$sig %in% res$index)) # 26 present in our set (full)

# save the subset of 26 out for colour coding results 
save <- zag[which(zag$sig %in% res$index),]
write.csv(save, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Results_collation/Replication/22_rep_stringent_second_threshold.csv", row.names = F)

# What about at EWAS levels of significance 
res <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Results_collation/FULL_EWAS_results_3.6_threshold_020221.csv")
res <- res[which(res$SeqId %in% list),]
res$index <- paste(res$SeqId, res$Probe, sep = "_")

length(which(zag$sig %in% res$index)) # 42 replicate at lower trheshold of significance
rep <- zag[which(zag$sig %in% res$index),]

# Make a replication table 
table <- zag[c(2:7,9)]
names(res)[13] <- "sig"

table <- left_join(zag, res, by= "sig")

# save a copy of this table 
table <- table[c(3,2,4,6,7,5,16,17,18,21)]
colnames(table) <- c("SeqId", "Gene of Protein", "CpG", "Zaghlool Beta", "Zaghlool SE", "Zaghlool P", "Beta", "SE", "P", "UniProt Full Name")
write.csv(table, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Results_collation/Replication/42_rep_3_6_threshold.csv", row.names = F)

### Now we will need to find the effect stats for the associations which do not pass significance for the proteins in the table remaining 

# ill need to fish out the EWAS results for these proteins now 

# 3216-2 - 1 association - this protein is in batch 7 
p1 <- fread("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/04_FULL/batches/batch7/3216-2_MOA_774_with_GRM.mlma")
p1$SeqId <- "3216-2"

# 3485-28 - 2 associations - batch 8
p2 <- fread("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/04_FULL/batches/batch8/3485-28_MOA_774_with_GRM.mlma")
p2$SeqId <- "3485-28"

# 4990-87 - 1 association - batch 8 
p3 <- fread("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/04_FULL/batches/batch8/4990-87_MOA_774_with_GRM.mlma")
p3$SeqId <- "4990-87"

# 4148-49 - rest (32) are this protein - batch 8 
p4 <- fread("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/04_FULL/batches/batch8/4148-49_MOA_774_with_GRM.mlma")
p4$SeqId <- "4148-49"


# so of the 36, there were 32 for PAPPA and 4148-49 that didnt replicate 

# Im going to join the tables, name protein_cpg column and see if we can match in the stats for our replication table 

p <- rbind(p1,p2)
p <- rbind(p,p3)
p <- rbind(p,p4)

p2 <- as.data.frame(p)

# get joiner column
table$assoc <- paste(table$SeqId, table$CpG, sep = "_")
p2$assoc <- paste(p2$SeqId, p2$Probe, sep = "_")

# name appropriately 
p2 <- p2[,c(10,6,7,8)]
names(p2) <- c("assoc", "b_non", "SE_non", "p_non")

# Can we find the associations of interest 
names(table)[10] <- "Uni"
sub <- table[which(table$Uni %in% NA),] # 39 associations to fill in 

length(which(sub$assoc %in% p2$assoc)) # 39 present 

# try to merge in

merge <- left_join(table, p2, by = "assoc")

# write out this file and examine to make sure we have results for every instance of the 81 comparable 
write.csv(merge, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Results_collation/Replication/lookup_of_nominal_EWAS_replication_second_threshold.csv", row.names = F)

# this will be our replication suppl table - now do the suppl figure to illustrate it too 

###########################################################################################

### FIGURE COMPARISON

# to do this, i'll read back in the file that i edited to collate the missing vs non missing data added in 
library(readxl)
merge <- read_excel("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Results_collation/Replication/rep_table.xlsx")
merge <- as.data.frame(merge)

### Plot the replication of associations

names(merge)[4] <- "zb"
names(merge)[5] <- "zb_se"
names(merge)[7] <- "ob"
names(merge)[8] <- "ob_se"

library(ggplot2)

# z scores
merge$z_zb <- (merge$zb) / (merge$zb_se)
merge$z_ob <- (merge$ob) / (merge$ob_se)

# Colour annotations
names(merge)[10] <- "rep"
merge$Rep_both <- ifelse(merge$rep == "P < 4.5x10-10", 1, 0)
merge$Rep_mid <- ifelse(merge$rep == "P < 3.6x10-8", 1, 0)
merge$Rep_nom <- ifelse(merge$rep == "Nominal P < 0.05", 1, 0)

# Plot the z-transformed effect sizes where comparisons were available across the studies 
pdf("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Results_collation/Replication/replication_plot_second_threshold.pdf", width = 5, height = 5)
ggplot(merge, aes(z_zb, z_ob), scale="globalminmax") +
       geom_vline(xintercept = 0, linetype = 1) +
       geom_hline(yintercept = 0, linetype = 1) +
       geom_point() +
       theme_minimal() +
      xlab("z-score effect Zaghlool et al") +
      ylab("z-score effect") +
       theme(plot.title = element_text(size = 14)) +
   geom_abline(slope=1, intercept=0, linetype = 2) +
  geom_point(data=merge[merge$Rep_both =="1",], shape=21, fill="skyblue2", size=2)+
  geom_point(data=merge[merge$Rep_mid =="1",], shape=21, fill="palegreen3", size=2)+
  geom_point(data=merge[merge$Rep_nom =="1",], shape=21, fill="tan1", size=2) 
dev.off()


# + xlim(-2.5,2.5) + ylim(-2.5,2.5)

##############################################################################

### NOW LOOK AT LAMBDA VALUE CALCULATION FOR THE FULL MODELS   

##############################################################################

## OSCA Lambdas 

screen

R

library(data.table)

cpgs <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Results_collation/Formatted_tables/FULL_EWAS_results_020221_formatted_second_threshold.csv")
length <- dim(cpgs)[1]
list = unique(cpgs$SeqId) # 146 proteins to get lambdas for 
names(cpgs)[10] <- "Somamer"

output <- matrix(nrow= length, ncol = 3)
output <- as.data.frame(output)
names(output) <- c("Protein", "Gene of Protein", "Lambda")

# Run for each batch to get lambdas and join at the end 

for(i in 1:length(list)){ 
  tryCatch({ 
    ewas = fread(paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/04_FULL/batches/batch1/", list[[i]], "_MOA_774_with_GRM.mlma"))
    ewas = as.data.frame(ewas)

    gene = cpgs[cpgs$SeqId %in% list[[i]],"Somamer"]
    gene = unique(gene)

    ## Calculate Lambda 
    # For p-values, calculate chi-squared statistic
    chisq <- qchisq(1-ewas$p,1)

    ## Calculate lambda gc (λgc)
    lambda <- median(chisq)/qchisq(0.5,1)

    output[i, 1] <- list[[i]]
    output[i, 2] <- gene
    output[i, 3] <- lambda 
  }, error = function(e) cat("skipped"))
} 

write.csv(output, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Results_collation/lambdas/batches/batch1_inflation.csv", row.names = F)


for(i in 1:length(list)){ 
  tryCatch({ 
    ewas = fread(paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/04_FULL/batches/batch2/", list[[i]], "_MOA_774_with_GRM.mlma"))
    ewas = as.data.frame(ewas)

    gene = cpgs[cpgs$SeqId %in% list[[i]],"Somamer"]
    gene = unique(gene)

    ## Calculate Lambda 
    # For p-values, calculate chi-squared statistic
    chisq <- qchisq(1-ewas$p,1)

    ## Calculate lambda gc (λgc)
    lambda <- median(chisq)/qchisq(0.5,1)

    output[i, 1] <- list[[i]]
    output[i, 2] <- gene
    output[i, 3] <- lambda 
  }, error = function(e) cat("skipped"))
} 

write.csv(output, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Results_collation/lambdas/batches/batch2_inflation.csv", row.names = F)


for(i in 1:length(list)){ 
  tryCatch({ 
    ewas = fread(paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/04_FULL/batches/batch3/", list[[i]], "_MOA_774_with_GRM.mlma"))
    ewas = as.data.frame(ewas)

    gene = cpgs[cpgs$SeqId %in% list[[i]],"Somamer"]
    gene = unique(gene)

    ## Calculate Lambda 
    # For p-values, calculate chi-squared statistic
    chisq <- qchisq(1-ewas$p,1)

    ## Calculate lambda gc (λgc)
    lambda <- median(chisq)/qchisq(0.5,1)

    output[i, 1] <- list[[i]]
    output[i, 2] <- gene
    output[i, 3] <- lambda 
  }, error = function(e) cat("skipped"))
} 

write.csv(output, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Results_collation/lambdas/batches/batch3_inflation.csv", row.names = F)


for(i in 1:length(list)){ 
  tryCatch({ 
    ewas = fread(paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/04_FULL/batches/batch4/", list[[i]], "_MOA_774_with_GRM.mlma"))
    ewas = as.data.frame(ewas)

    gene = cpgs[cpgs$SeqId %in% list[[i]],"Somamer"]
    gene = unique(gene)

    ## Calculate Lambda 
    # For p-values, calculate chi-squared statistic
    chisq <- qchisq(1-ewas$p,1)

    ## Calculate lambda gc (λgc)
    lambda <- median(chisq)/qchisq(0.5,1)

    output[i, 1] <- list[[i]]
    output[i, 2] <- gene
    output[i, 3] <- lambda 
  }, error = function(e) cat("skipped"))
} 

write.csv(output, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Results_collation/lambdas/batches/batch4_inflation.csv", row.names = F)


for(i in 1:length(list)){ 
  tryCatch({ 
    ewas = fread(paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/04_FULL/batches/batch5/", list[[i]], "_MOA_774_with_GRM.mlma"))
    ewas = as.data.frame(ewas)

    gene = cpgs[cpgs$SeqId %in% list[[i]],"Somamer"]
    gene = unique(gene)

    ## Calculate Lambda 
    # For p-values, calculate chi-squared statistic
    chisq <- qchisq(1-ewas$p,1)

    ## Calculate lambda gc (λgc)
    lambda <- median(chisq)/qchisq(0.5,1)

    output[i, 1] <- list[[i]]
    output[i, 2] <- gene
    output[i, 3] <- lambda 
  }, error = function(e) cat("skipped"))
} 

write.csv(output, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Results_collation/lambdas/batches/batch5_inflation.csv", row.names = F)


for(i in 1:length(list)){ 
  tryCatch({ 
    ewas = fread(paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/04_FULL/batches/batch6/", list[[i]], "_MOA_774_with_GRM.mlma"))
    ewas = as.data.frame(ewas)

    gene = cpgs[cpgs$SeqId %in% list[[i]],"Somamer"]
    gene = unique(gene)

    ## Calculate Lambda 
    # For p-values, calculate chi-squared statistic
    chisq <- qchisq(1-ewas$p,1)

    ## Calculate lambda gc (λgc)
    lambda <- median(chisq)/qchisq(0.5,1)

    output[i, 1] <- list[[i]]
    output[i, 2] <- gene
    output[i, 3] <- lambda 
  }, error = function(e) cat("skipped"))
} 

write.csv(output, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Results_collation/lambdas/batches/batch6_inflation.csv", row.names = F)


for(i in 1:length(list)){ 
  tryCatch({ 
    ewas = fread(paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/04_FULL/batches/batch7/", list[[i]], "_MOA_774_with_GRM.mlma"))
    ewas = as.data.frame(ewas)

    gene = cpgs[cpgs$SeqId %in% list[[i]],"Somamer"]
    gene = unique(gene)

    ## Calculate Lambda 
    # For p-values, calculate chi-squared statistic
    chisq <- qchisq(1-ewas$p,1)

    ## Calculate lambda gc (λgc)
    lambda <- median(chisq)/qchisq(0.5,1)

    output[i, 1] <- list[[i]]
    output[i, 2] <- gene
    output[i, 3] <- lambda 
  }, error = function(e) cat("skipped"))
} 

write.csv(output, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Results_collation/lambdas/batches/batch7_inflation.csv", row.names = F)


for(i in 1:length(list)){ 
  tryCatch({ 
    ewas = fread(paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/04_FULL/batches/batch8/", list[[i]], "_MOA_774_with_GRM.mlma"))
    ewas = as.data.frame(ewas)

    gene = cpgs[cpgs$SeqId %in% list[[i]],"Somamer"]
    gene = unique(gene)

    ## Calculate Lambda 
    # For p-values, calculate chi-squared statistic
    chisq <- qchisq(1-ewas$p,1)

    ## Calculate lambda gc (λgc)
    lambda <- median(chisq)/qchisq(0.5,1)

    output[i, 1] <- list[[i]]
    output[i, 2] <- gene
    output[i, 3] <- lambda 
  }, error = function(e) cat("skipped"))
} 

write.csv(output, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Results_collation/lambdas/batches/batch8_inflation.csv", row.names = F)


for(i in 1:length(list)){ 
  tryCatch({ 
    ewas = fread(paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/04_FULL/batches/batch9/", list[[i]], "_MOA_774_with_GRM.mlma"))
    ewas = as.data.frame(ewas)

    gene = cpgs[cpgs$SeqId %in% list[[i]],"Somamer"]
    gene = unique(gene)

    ## Calculate Lambda 
    # For p-values, calculate chi-squared statistic
    chisq <- qchisq(1-ewas$p,1)

    ## Calculate lambda gc (λgc)
    lambda <- median(chisq)/qchisq(0.5,1)

    output[i, 1] <- list[[i]]
    output[i, 2] <- gene
    output[i, 3] <- lambda 
  }, error = function(e) cat("skipped"))
} 

write.csv(output, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Results_collation/lambdas/batches/batch9_inflation.csv", row.names = F)


for(i in 1:length(list)){ 
  tryCatch({ 
    ewas = fread(paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/04_FULL/batches/batch10/", list[[i]], "_MOA_774_with_GRM.mlma"))
    ewas = as.data.frame(ewas)

    gene = cpgs[cpgs$SeqId %in% list[[i]],"Somamer"]
    gene = unique(gene)

    ## Calculate Lambda 
    # For p-values, calculate chi-squared statistic
    chisq <- qchisq(1-ewas$p,1)

    ## Calculate lambda gc (λgc)
    lambda <- median(chisq)/qchisq(0.5,1)

    output[i, 1] <- list[[i]]
    output[i, 2] <- gene
    output[i, 3] <- lambda 
  }, error = function(e) cat("skipped"))
} 

write.csv(output, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Results_collation/lambdas/batches/batch10_inflation.csv", row.names = F)


for(i in 1:length(list)){ 
  tryCatch({ 
    ewas = fread(paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/04_FULL/batches/batch11/", list[[i]], "_MOA_774_with_GRM.mlma"))
    ewas = as.data.frame(ewas)

    gene = cpgs[cpgs$SeqId %in% list[[i]],"Somamer"]
    gene = unique(gene)

    ## Calculate Lambda 
    # For p-values, calculate chi-squared statistic
    chisq <- qchisq(1-ewas$p,1)

    ## Calculate lambda gc (λgc)
    lambda <- median(chisq)/qchisq(0.5,1)

    output[i, 1] <- list[[i]]
    output[i, 2] <- gene
    output[i, 3] <- lambda 
  }, error = function(e) cat("skipped"))
} 

write.csv(output, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Results_collation/lambdas/batches/batch11_inflation.csv", row.names = F)


for(i in 1:length(list)){ 
  tryCatch({ 
    ewas = fread(paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/04_FULL/batches/batch12/", list[[i]], "_MOA_774_with_GRM.mlma"))
    ewas = as.data.frame(ewas)

    gene = cpgs[cpgs$SeqId %in% list[[i]],"Somamer"]
    gene = unique(gene)

    ## Calculate Lambda 
    # For p-values, calculate chi-squared statistic
    chisq <- qchisq(1-ewas$p,1)

    ## Calculate lambda gc (λgc)
    lambda <- median(chisq)/qchisq(0.5,1)

    output[i, 1] <- list[[i]]
    output[i, 2] <- gene
    output[i, 3] <- lambda 
  }, error = function(e) cat("skipped"))
} 

write.csv(output, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Results_collation/lambdas/batches/batch12_inflation.csv", row.names = F)




## Collate the lambda values across batches to get the full set 
setwd("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Results_collation/lambdas/batches/")
L <- list.files("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Results_collation/lambdas/batches/", ".")

files <- lapply(L, read.csv)
files <- na.omit(files)
osca <- do.call(rbind, files)
osca <- na.omit(osca)
names(osca) <- c("Protein SeqId", "Gene of Protein", "Lambda")

# write.csv(osca, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_180821_results/lambdas_no_eGFR/joint/joint.csv", row.names = F)

osca2 <- osca[which(osca[,1] %in% list),]

unique(osca2[,1])

library(tidyverse)
osca3 <- osca2 %>% group_by("Protein SeqId")

# It seems we are duplicating rows that are identical  -remove duplicates

table(duplicated(osca))

unique <- unique(osca)
dim(unique) # 146 proteins - as expected 

unique <- unique[order(unique$Lambda),]

# Add gene names in from annotation file and reorder for saving 
anno <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Annotations_protein_file/annotation_formatted_for_paper.csv")
names(anno)[1] <- "Protein"
names(unique)[1] <- "Protein"
unique <- left_join(unique, anno, by = "Protein")
unique <- unique[c(1,3,4,5)]

write.csv(unique, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Results_collation/lambdas/joint_second_threshold.csv", row.names = F)


##########################################################################################################################

