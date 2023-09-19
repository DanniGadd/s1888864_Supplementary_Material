############################################################################################

## GRM GENERATION IN STRADL 

############################################################################################

# The following steps have been performed in this script to generate a GRM in STRADL:
# 1 - Recode SNP chip data to 0,1,2 for the 774 individuals (IDs without GS suffix)
# 2 - Generate myprofile.txt from recoded dta
# 3 - Generate orm from -efile (myprofile.txt) using --make-orm
# 4 - Edit orm.id file in R, appending GS_ suffixes
# 5 - Run moa with DNAm data, fitting orm
# One difference between mine and yours is I just used the chip data - did you use the imputed data for your grm?
# I'm not so sure I've even done it right. OSCA documentation only really talks about gene expression/DNAm data
# I just substituted it with a matrix of 0,1,2 genotype counts
# They reference GWAS data for things like eQTL mapping

setwd("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/GRM/")

# Previously 778 samples - now 774 samples
# Old file: o778 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_080421/test_files_comparison_090421/record_of_IDS_for_methylation_778_meth_order.csv")
o778 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/774_IDs_meth_order.csv")
o778 <- o778[c(3,5,4,1,2)]

keep_778  = o778  
keep_778[,2] = gsub("GS_", "", keep_778[,2])
keep_778[,3] = gsub("GS_", "", o778[,3])
write.table(keep_778[,2:3], file="keep_774.txt", sep=' ', quote=F, row.names=F, col.names=F)

# 778 SNP prep to orm

# Move to cluster and enter the following from terminal 

cd /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/GRM/

# plink19 --recodeA --bfile /GWAS_Source/GS_GWAS/GS20K_QC\'d/GS20K_QCpsychprotocol_SPH_04112015 --out osca_test_774 --keep keep_774.txt

# Move back to R and enter the following 

library(data.table)
snps = fread("osca_test_774.raw", header=T, stringsAsFactors=F)
snps = as.data.frame(snps)
snps[,3:6] = NULL

write.table(snps, file="myprofile_774.txt", row.names=F, sep=' ')

# Move back to cluster and enter the following in terminal 

cd /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/GRM/
osca_Linux --efile myprofile_774.txt --keep keep_774.txt --make-orm --out ormtest774

# Run quick test to ensure GRM format works with MOA
# i=/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_080421/relatedness_testing_120421/test_1_778/2780-35eGFR_included_778.phen
# osca_Linux --moa --befile /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_080421/test_files_comparison_090421/methylation_778_meth_order --pheno $i --qcovar /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_080421/test_files_comparison_090421/WBC_covariates_778_main_model.cov --fast-linear --out test_MOA_778_output --methylation-m
# osca_Linux --moa --efile myprofile.txt --pheno $i --qcovar /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_080421/test_files_comparison_090421/WBC_covariates_778_main_model.cov --fast-linear --out test1

# In R (set IDs to GS_XXX to match DNAm IDs)
orm = read.table("ormtest774.orm.id")
orm[,1] = paste0("GS_", orm[,1])
orm[,2] = paste0("GS_", orm[,2])

write.table(orm, file="ormtest774.orm.id", sep='\t', col.names=F, row.names=F, quote=F)

