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
# EXTRA: WBC estimates have been added as covariates 

# This model will produce the version of EWAS that represents the WBC adjustment.


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
  meth[,i] <- resid(lm(meth[,i] ~ phenotypes$age_stradl + as.factor(phenotypes$sex) + as.factor(phenotypes$Batch) + as.factor(phenotypes$wave) + as.factor(phenotypes$combined) + phenotypes$Bcell + phenotypes$CD4T + phenotypes$CD8T + phenotypes$Gran + phenotypes$NK + phenotypes$Mono + phenotypes$Eos, na.action = na.exclude))
}

meth<-as.data.frame(meth)

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
fwrite(tmp, file="/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Methylation_preps/methylation_774_meth_order_WBC_310121.txt", row.names=F, quote = F, sep=' ')


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
location <- "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/03_WBCs/batches/"

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
osca_Linux --efile /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Methylation_preps/methylation_774_meth_order_WBC_310121.txt --methylation-m --make-bod --out /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Methylation_preps/methylation_774_meth_order_WBC_310121

# Regenerate the opi file with correct annotations 
anno = readRDS("/Cluster_Filespace/Marioni_Group/Daniel/EPIC_AnnotationObject_df.rds")
opi1 = read.table("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Methylation_preps/methylation_774_meth_order_WBC_310121.opi", header=F, stringsAsFactors=F)
opi <- anno[opi1$V2,c("chr","Name", "pos","UCSC_RefGene_Name", "strand")]
opi$chr <- gsub("chr", "", opi$chr)
opi$chr <- gsub("X", "23", opi$chr)
opi$chr <- gsub("Y", "24", opi$chr)
opi$chr <- as.numeric(opi$chr)
opi$UCSC_RefGene_Name  <- as.factor(opi$UCSC_RefGene_Name )
opi$strand <- as.factor(opi$strand)
opi[which(opi$UCSC_RefGene_Name ==""), "UCSC_RefGene_Name"] <- NA
write.table(opi, file="/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Methylation_preps/methylation_774_meth_order_WBC_310121.opi",  # make sure to keep the same filename as before
                col.names=F, 
                row.names=F, 
                quote=F, sep='\t')

# check that the annotations are present after the update
head /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Methylation_preps/methylation_774_meth_order_WBC_310121.opi


##########################################################################################

# BATCH EWAS 

cd /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/03_WBCs/batches/batch1/
for i in *.phen 

do

  A=$( echo $i | cut -d"/" -f3)
  B=$( echo $A | cut -d. -f1)

osca_Linux --moa --orm /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/GRM/ormtest774 --befile /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Methylation_preps/methylation_774_meth_order_WBC_310121 --pheno $i --out ${B}_MOA_774_with_GRM --methylation-m


done 


cd /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/03_WBCs/batches/batch2/
for i in *.phen 

do

  A=$( echo $i | cut -d"/" -f3)
  B=$( echo $A | cut -d. -f1)

osca_Linux --moa --orm /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/GRM/ormtest774 --befile /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Methylation_preps/methylation_774_meth_order_WBC_310121 --pheno $i --out ${B}_MOA_774_with_GRM --methylation-m


done 


cd /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/03_WBCs/batches/batch3/
for i in *.phen 

do

  A=$( echo $i | cut -d"/" -f3)
  B=$( echo $A | cut -d. -f1)

osca_Linux --moa --orm /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/GRM/ormtest774 --befile /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Methylation_preps/methylation_774_meth_order_WBC_310121 --pheno $i --out ${B}_MOA_774_with_GRM --methylation-m


done 


cd /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/03_WBCs/batches/batch4/
for i in *.phen 

do

  A=$( echo $i | cut -d"/" -f3)
  B=$( echo $A | cut -d. -f1)

osca_Linux --moa --orm /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/GRM/ormtest774 --befile /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Methylation_preps/methylation_774_meth_order_WBC_310121 --pheno $i --out ${B}_MOA_774_with_GRM --methylation-m


done 


cd /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/03_WBCs/batches/batch5/
for i in *.phen 

do

  A=$( echo $i | cut -d"/" -f3)
  B=$( echo $A | cut -d. -f1)

osca_Linux --moa --orm /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/GRM/ormtest774 --befile /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Methylation_preps/methylation_774_meth_order_WBC_310121 --pheno $i --out ${B}_MOA_774_with_GRM --methylation-m


done 


cd /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/03_WBCs/batches/batch6/
for i in *.phen 

do

  A=$( echo $i | cut -d"/" -f3)
  B=$( echo $A | cut -d. -f1)

osca_Linux --moa --orm /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/GRM/ormtest774 --befile /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Methylation_preps/methylation_774_meth_order_WBC_310121 --pheno $i --out ${B}_MOA_774_with_GRM --methylation-m


done 


cd /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/03_WBCs/batches/batch7/
for i in *.phen 

do

  A=$( echo $i | cut -d"/" -f3)
  B=$( echo $A | cut -d. -f1)

osca_Linux --moa --orm /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/GRM/ormtest774 --befile /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Methylation_preps/methylation_774_meth_order_WBC_310121 --pheno $i --out ${B}_MOA_774_with_GRM --methylation-m


done 


cd /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/03_WBCs/batches/batch8/
for i in *.phen 

do

  A=$( echo $i | cut -d"/" -f3)
  B=$( echo $A | cut -d. -f1)

osca_Linux --moa --orm /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/GRM/ormtest774 --befile /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Methylation_preps/methylation_774_meth_order_WBC_310121 --pheno $i --out ${B}_MOA_774_with_GRM --methylation-m


done 


cd /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/03_WBCs/batches/batch9/
for i in *.phen 

do

  A=$( echo $i | cut -d"/" -f3)
  B=$( echo $A | cut -d. -f1)

osca_Linux --moa --orm /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/GRM/ormtest774 --befile /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Methylation_preps/methylation_774_meth_order_WBC_310121 --pheno $i --out ${B}_MOA_774_with_GRM --methylation-m


done 


cd /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/03_WBCs/batches/batch10/
for i in *.phen 

do

  A=$( echo $i | cut -d"/" -f3)
  B=$( echo $A | cut -d. -f1)

osca_Linux --moa --orm /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/GRM/ormtest774 --befile /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Methylation_preps/methylation_774_meth_order_WBC_310121 --pheno $i --out ${B}_MOA_774_with_GRM --methylation-m


done 


cd /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/03_WBCs/batches/batch11/
for i in *.phen 

do

  A=$( echo $i | cut -d"/" -f3)
  B=$( echo $A | cut -d. -f1)

osca_Linux --moa --orm /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/GRM/ormtest774 --befile /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Methylation_preps/methylation_774_meth_order_WBC_310121 --pheno $i --out ${B}_MOA_774_with_GRM --methylation-m


done 


cd /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/03_WBCs/batches/batch12/
for i in *.phen 

do

  A=$( echo $i | cut -d"/" -f3)
  B=$( echo $A | cut -d. -f1)

osca_Linux --moa --orm /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/GRM/ormtest774 --befile /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Methylation_preps/methylation_774_meth_order_WBC_310121 --pheno $i --out ${B}_MOA_774_with_GRM --methylation-m


done 



####################################################################################################

### PROCESSING THE EWAS RESULTS FOR BATCHES ABOVE - MOA 778 with genetic ORM - WBC - AND pQTLs - but no smoking/BMI

####################################################################################################

# Split this up for speed into 12 screens 

# Go through each protein and subset to the significant sites 


screen

R


library(data.table)

list <- c("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/03_WBCs/batches/batch1/")

sig_hits <- "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/03_WBCs/batches/sig_hits/"

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

list <- c("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/03_WBCs/batches/batch2/")

sig_hits <- "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/03_WBCs/batches/sig_hits/"


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

list <- c("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/03_WBCs/batches/batch3/")

sig_hits <- "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/03_WBCs/batches/sig_hits/"

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

list <- c("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/03_WBCs/batches/batch4/")

sig_hits <- "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/03_WBCs/batches/sig_hits/"


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

list <- c("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/03_WBCs/batches/batch5/")

sig_hits <- "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/03_WBCs/batches/sig_hits/"


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

list <- c("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/03_WBCs/batches/batch6/")

sig_hits <- "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/03_WBCs/batches/sig_hits/"


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

list <- c("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/03_WBCs/batches/batch7/")

sig_hits <- "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/03_WBCs/batches/sig_hits/"


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

list <- c("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/03_WBCs/batches/batch8/")

sig_hits <- "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/03_WBCs/batches/sig_hits/"


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

list <- c("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/03_WBCs/batches/batch9/")

sig_hits <- "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/03_WBCs/batches/sig_hits/"


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

list <- c("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/03_WBCs/batches/batch10/")

sig_hits <- "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/03_WBCs/batches/sig_hits/"


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

list <- c("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/03_WBCs/batches/batch11/")

sig_hits <- "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/03_WBCs/batches/sig_hits/"


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

list <- c("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/03_WBCs/batches/batch12/")

sig_hits <- "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/03_WBCs/batches/sig_hits/"


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

path <- "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/03_WBCs/batches/sig_hits/"

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

# Filter to CpG-protein correction value of adjusted p 0.05/143/772619 = 4.525521e-10
# 0.0000000004525521

comb2_filt2 <- comb[which(comb$p < 0.0000000004525521),]

dim(comb2_filt2) # 3213  - at cpg/protein adjusted significance 
comb2 <- comb2_filt2[order(comb2_filt2$p),]
write.csv(comb2, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Results_collation/WBC_EWAS_results_second_threshold_020221.csv", row.names = F)


# Filter to CpG-protein correction value of adjusted p 0.05/4235/772619 = 1.528098e-11
# > 0.05 / 4235 / 772619
# [1] 1.528098e-11

# 0.00000000001528098

comb2_filt2 <- comb[which(comb$p < 0.00000000001528098),]

dim(comb2_filt2) # 2230  - at cpg/protein adjusted significance 
comb2 <- comb2_filt2[order(comb2_filt2$p),]
write.csv(comb2, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Results_collation/WBC_EWAS_results_020221.csv", row.names = F)

##############################################################

### WORK OUT WHICH FULL PROTEINS DIDNT CONVERGE 

# Read in protein file which has GS id and stradl id, then all 4,235 proteins that have been preprocessed
pheno1 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/prep/proteins_4235_778_eGFR_included.csv", check.names = F)
pheno2 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/prep/Phenotype_file_778_eGFR_removed.csv", check.names = F)
names(pheno2)[1] <- "ID"
pheno2 <- pheno2[-2]
colnames(pheno2) <- colnames(pheno1)


path <- "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/03_WBCs/batches/sig_hits/"

setwd(path)

L <- list.files(".", ".csv")
L

files <- lapply(L, read.csv)
names <- as.character(L)
batch <- gsub("_.*", "", names)
marker <- gsub("_results.*", "", names)
marker <- gsub(".*_", "", marker)
names(files) <- marker
osca <- do.call(rbind, files)
osca <- osca[c(9,1:8)]

# load pheno2 above and do this to work out difference 
list <- names(files)
list2 <- pheno2[-which(colnames(pheno2) %in% list)]

# Plots - of proteins look okay 
proteins <- colnames(list2)[1:8]
proteins





