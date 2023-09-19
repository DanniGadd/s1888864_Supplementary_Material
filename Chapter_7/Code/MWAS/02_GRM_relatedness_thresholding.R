############################################################################################

## Identify related individuals based on 5% upper limits of GRM in STRADL 

############################################################################################

# This script generates lists with FID and IID and GS ID for the unrelated and related samples used

############################################################################################

## CHECK TO SEE WHO IS UNRELATED IN THE GENETIC SAMPLE

# Daniel sent this location for the GRM to use:
# /Cluster_Filespace/Marioni_Group/Daniel/Somalogic_GWAS/grm_1000_full

# Prune the GRM for relatedness by a cutoff of 0.05
gcta64 --grm /Cluster_Filespace/Marioni_Group/Daniel/Somalogic_GWAS/grm_1000_full --grm-cutoff 0.05 --make-grm --out /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/trimming_GRM_180321/grm_pruned_0.05_related

# 181 inidviduals removed from the GRM and 884 kept from the 1000 we had 

# Extract the GRM subject id of all the singletons by a cutoff of 0.05
gcta64 --grm /Cluster_Filespace/Marioni_Group/Daniel/Somalogic_GWAS/grm_1000_full --grm-singleton 0.05  --out /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/trimming_GRM_180321/grm_singletons

# Save IDs for the 884 unrelated in STRADL 
test <- read.delim("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/trimming_GRM_180321/grm_singletons.singleton.txt", header = F) # this is the IDs of all my 884 single people 

############################################################################################

## NOW GET A LIST OF THE UNRELATED INDIVIDUALS IN THE MWAS STUDY SAMPLE

# read in protein example dataset with GS id listed for the 778 individuals 
# pheno2 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/prep/proteins_the_whole_lot_total.csv", check.names = F)
pheno2 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/774_ids_MWAS.csv") # from 00_cell_proportion_estimates script in MWAS repo

# Add FAM and ID columns 
id <- read.csv("/Cluster_Filespace/Marioni_Group/Rob/Somalogic/IDs_778.csv") # 778 people in EWAS
fam <- read.table("/Cluster_Filespace/Marioni_Group/Rob/Somalogic/GWAS_Somalogic_data.fam") # 1064 total LBC family IDs included, with other info 

## Make sure order of IDs in Family ID and regular ID file match 
# So the second column of the fam file is the participant ID, which we will make sure matches here 
fam1 = fam[fam$V2 %in% id$x,] # subset the fam file to the participant IDs which are in the ID file 
matcher = id$x # create order based on the id file 
fam1 = fam1[match(matcher,fam1$V2),] # match order of participant IDs between the 2 files (first column is fam ID)

# Name properly 
osca_dat <- fam1[c(1,2)]
names(osca_dat) <- c("FID", "IID")

# Remove those that are excluded due to not having depression status 
osca_dat <- osca_dat[-which(!osca_dat$IID %in% pheno2$GS_id),]

# Check to make sure the id in the proteins file is the same
identical(pheno2$GS_id, osca_dat$IID) # TRUE 

pheno2$IID <- osca_dat$IID
pheno2$FID <- osca_dat$FID

identical(pheno2$GS_id, pheno2$IID) # TRUE 

# Subset phenotype file to only those that are unrelated (658) 
# read in the unrelated IIDs generated above 
test <- read.delim("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/trimming_GRM_180321/grm_singletons.singleton.txt", header = F) # this is the IDs of all my 884 single people 
ov <- which(pheno2$IID %in% test$V2)
pheno3 <- pheno2[ov,]

# > dim(pheno2)
# [1]  654

# 774 - 654 = 120 related individuals in the MWAS

# write out the IDs for the 654 people 
pheno3$FID<-sub("^","GS_", pheno3$FID)
pheno3$IID<-sub("^","GS_", pheno3$IID)
write.csv(pheno3, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/654_unrelated_IDs.csv", row.names = F)

# write out the IDs for the 774 people 
pheno2$FID<-sub("^","GS_", pheno2$FID)
pheno2$IID<-sub("^","GS_", pheno2$IID)
write.csv(pheno2, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/774_IDs_meth_order.csv", row.names = F)


