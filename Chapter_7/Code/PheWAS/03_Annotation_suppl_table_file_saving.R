
### RECORD OF CHANGES MADE TO ANNOTATIONS FILE FOR STRADL 

# Danni - 270121
# Renaming gene names for 4 septin family proteins that were in date format in the original file 

screen

R

# Read annotation file into R
anno <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Annotations_protein_file/SOMAscan_Assay_v4_Annotations_version3.3.2_Septin_Name_Correction_270121.csv")
anno <- anno[c(1,18,13,4,2)]

# Read in list of 4235 SOMAmers 
pheno2 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/prep/proteins_4235_778_eGFR_included.csv", check.names = F)
names <- colnames(pheno2)
anno <- anno[which(anno$SeqId %in% names),]
anno <- anno[c(1,4,2,3)]

# > length(unique(anno[,3]))
# [1] 4058 - unique proteins 

write.csv(anno, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Annotations_protein_file/annotation_formatted_for_paper.csv", row.names = F)



