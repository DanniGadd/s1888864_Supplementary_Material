################################################################################################################################
################################################################################################################################
############################## APOE status generation in GS ####################################################################
################################################################################################################################
################################################################################################################################

# Aim - generate APOE and GrimAge variables in STRADL waves x2 
setwd("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Aging_sex_profiling/")

library(tidyverse)
library(readxl)

### generate apoe alleles and merge with cgs ###
a = read.table("/Cluster_Filespace/Marioni_Group/GS/GS_dataset/APOE_snps/rs7412.txt", header=T)
b = read.table("/Cluster_Filespace/Marioni_Group/GS/GS_dataset/APOE_snps/rs429358.txt", header=T)
names(a)[5] <- "rs7412"
names(b)[5] <- "rs429358"
apoe = merge(a[,c(1,5)],b[,c(1,5)],by=c("StudyID"))
apoe$apoe <-
ifelse((apoe$rs7412=="T" & apoe$rs429358=="T" ), "e2e2",
ifelse((apoe$rs7412=="C" & apoe$rs429358=="T" ), "e3e3",
ifelse((apoe$rs7412=="C" & apoe$rs429358=="C" ), "e4e4",
ifelse((apoe$rs7412=="Both" & apoe$rs429358=="T" ), "e2e3",
ifelse((apoe$rs7412=="Both" & apoe$rs429358=="Both" ), "e2e4",
ifelse((apoe$rs7412=="C" & apoe$rs429358=="Both" ), "e3e4", NA))))))

# Save out the apoe status file in correct format 
write.csv(apoe, "/Cluster_Filespace/Marioni_Group/GS/GS_dataset/APOE_snps/apoe_GS.csv", row.names = F)



