############################################################################################
############################################################################################
############################### pQTLs for EWAS extraction - GS #############################
############################################################################################
############################################################################################

############################################################################################################

### pQTL identification - load in and format relevant SNPs for extraction 

############################################################################################################

# Read in the sun pQTL data 
library(data.table)
library(tidyverse)
library(readxl)
list <- read_excel("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_pQTL_extraction_GS/pQTL_list.xlsx")
list <- as.data.frame(list)

# Write out so we can join to suppl files from EWAS 
list2 <- list[-1,]
write.csv(list2, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Interpretation/interpretation_090621/sun_pQTLs.csv", row.names = F)

# write.table(snplist, file="snplist.txt", quote=F, row.names=F, col.names=F) - snplist is a df containing your SNP IDs
# read in the .bim file and your pqtl summary data and paste your chr:pos:allele data together
# check ordering of chr:pos:A1:A2 against your .bim file and switch it to chr:pos:A2:A1 if necessary
# pQTL annotations and pull out chr/pos/allele1/allele2

## Take a look at imputed data format 
head /Cluster_Filespace/Marioni_Group/Daniel/GWAS_AgeAccel/Meta_Analysis_Daniel/cojo_files/GS_HRC/GS20K_chr10_HRC.r1-1_nomono_I4_cpra.bim
# 10      10_60494_A_G    0       60494   G       A
# chr   chr_pos_allele2_allele1   pos    allele1  allele2 
# So for this SNP ill try to extract 10_60494_A_G and 10_60494_G_A to cover all bases 

## Format the SNP list for extraction such that it matches the imputed SNP format 
# First do the A1A2 format 
order <- list[c(7,8,11,12)]
table(is.na(order)) # only 4 entries from the first row to be removed 
order <- na.omit(order)
order <- order %>% unite("new", 1:4, remove = FALSE)

# Next do the reverse format for the alleles swapped around
order2 <- list[c(7,8,12,11)]
table(is.na(order2)) # only 4 entries from the first row to be removed 
order2 <- na.omit(order2)
order2 <- order2 %>% unite("new", 1:4, remove = FALSE)
write.table(order2, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_pQTL_extraction_GS/chr_extraction_imputed/order2_list.txt", quote = F, row.names = F, col.names = F)


# Now get both sets of lists and join together to save out for extraction list
order3 <- rbind(order, order2)
snplist <- order3[1]
write.table(snplist, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_pQTL_extraction_GS/chr_extraction_imputed/snplist.txt", quote = F, row.names = F, col.names = F)

###################################################################################################################

### IMPUTED DATA EXTRACTION

###################################################################################################################

### I've done each chr separately as I couldnt get the loop working and wanted to check each one at first 

cd /Cluster_Filespace/Marioni_Group/Daniel/GWAS_AgeAccel/Meta_Analysis_Daniel/cojo_files/GS_HRC/
chrs <- c("GS20K_chr1_HRC.r1-1_nomono_I4_cpra")
plink19 --recodeA --bfile GS20K_chr1_HRC.r1-1_nomono_I4_cpra --extract /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_pQTL_extraction_GS/chr_extraction_imputed/snplist.txt --out /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_pQTL_extraction_GS/chr_extraction_imputed/run_all_chr/chr_1_
chrs <- c("GS20K_chr2_HRC.r1-1_nomono_I4_cpra")
plink19 --recodeA --bfile GS20K_chr2_HRC.r1-1_nomono_I4_cpra --extract /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_pQTL_extraction_GS/chr_extraction_imputed/snplist.txt --out /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_pQTL_extraction_GS/chr_extraction_imputed/run_all_chr/chr_2_
chrs <- c("GS20K_chr3_HRC.r1-1_nomono_I4_cpra")
plink19 --recodeA --bfile GS20K_chr3_HRC.r1-1_nomono_I4_cpra --extract /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_pQTL_extraction_GS/chr_extraction_imputed/snplist.txt --out /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_pQTL_extraction_GS/chr_extraction_imputed/run_all_chr/chr_3_
chrs <- c("GS20K_chr4_HRC.r1-1_nomono_I4_cpra")
plink19 --recodeA --bfile GS20K_chr4_HRC.r1-1_nomono_I4_cpra --extract /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_pQTL_extraction_GS/chr_extraction_imputed/snplist.txt --out /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_pQTL_extraction_GS/chr_extraction_imputed/run_all_chr/chr_4_
chrs <- c("GS20K_chr5_HRC.r1-1_nomono_I4_cpra")
plink19 --recodeA --bfile GS20K_chr5_HRC.r1-1_nomono_I4_cpra --extract /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_pQTL_extraction_GS/chr_extraction_imputed/snplist.txt --out /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_pQTL_extraction_GS/chr_extraction_imputed/run_all_chr/chr_5_
chrs <- c("GS20K_chr6_HRC.r1-1_nomono_I4_cpra")
plink19 --recodeA --bfile GS20K_chr6_HRC.r1-1_nomono_I4_cpra --extract /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_pQTL_extraction_GS/chr_extraction_imputed/snplist.txt --out /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_pQTL_extraction_GS/chr_extraction_imputed/run_all_chr/chr_6_
chrs <- c("GS20K_chr7_HRC.r1-1_nomono_I4_cpra")
plink19 --recodeA --bfile GS20K_chr7_HRC.r1-1_nomono_I4_cpra --extract /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_pQTL_extraction_GS/chr_extraction_imputed/snplist.txt --out /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_pQTL_extraction_GS/chr_extraction_imputed/run_all_chr/chr_7_
chrs <- c("GS20K_chr8_HRC.r1-1_nomono_I4_cpra")
plink19 --recodeA --bfile GS20K_chr8_HRC.r1-1_nomono_I4_cpra --extract /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_pQTL_extraction_GS/chr_extraction_imputed/snplist.txt --out /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_pQTL_extraction_GS/chr_extraction_imputed/run_all_chr/chr_8_
chrs <- c("GS20K_chr9_HRC.r1-1_nomono_I4_cpra")
plink19 --recodeA --bfile GS20K_chr9_HRC.r1-1_nomono_I4_cpra --extract /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_pQTL_extraction_GS/chr_extraction_imputed/snplist.txt --out /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_pQTL_extraction_GS/chr_extraction_imputed/run_all_chr/chr_9_
chrs <- c("GS20K_chr10_HRC.r1-1_nomono_I4_cpra")
plink19 --recodeA --bfile GS20K_chr10_HRC.r1-1_nomono_I4_cpra --extract /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_pQTL_extraction_GS/chr_extraction_imputed/snplist.txt --out /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_pQTL_extraction_GS/chr_extraction_imputed/run_all_chr/chr_10_
chrs <- c("GS20K_chr11_HRC.r1-1_nomono_I4_cpra")
plink19 --recodeA --bfile GS20K_chr11_HRC.r1-1_nomono_I4_cpra --extract /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_pQTL_extraction_GS/chr_extraction_imputed/snplist.txt --out /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_pQTL_extraction_GS/chr_extraction_imputed/run_all_chr/chr_11_
chrs <- c("GS20K_chr12_HRC.r1-1_nomono_I4_cpra")
plink19 --recodeA --bfile GS20K_chr12_HRC.r1-1_nomono_I4_cpra --extract /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_pQTL_extraction_GS/chr_extraction_imputed/snplist.txt --out /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_pQTL_extraction_GS/chr_extraction_imputed/run_all_chr/chr_12_
chrs <- c("GS20K_chr13_HRC.r1-1_nomono_I4_cpra")
plink19 --recodeA --bfile GS20K_chr13_HRC.r1-1_nomono_I4_cpra --extract /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_pQTL_extraction_GS/chr_extraction_imputed/snplist.txt --out /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_pQTL_extraction_GS/chr_extraction_imputed/run_all_chr/chr_13_
chrs <- c("GS20K_chr14_HRC.r1-1_nomono_I4_cpra")
plink19 --recodeA --bfile GS20K_chr14_HRC.r1-1_nomono_I4_cpra --extract /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_pQTL_extraction_GS/chr_extraction_imputed/snplist.txt --out /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_pQTL_extraction_GS/chr_extraction_imputed/run_all_chr/chr_14_
chrs <- c("GS20K_chr15_HRC.r1-1_nomono_I4_cpra")
plink19 --recodeA --bfile GS20K_chr15_HRC.r1-1_nomono_I4_cpra --extract /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_pQTL_extraction_GS/chr_extraction_imputed/snplist.txt --out /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_pQTL_extraction_GS/chr_extraction_imputed/run_all_chr/chr_15_
chrs <- c("GS20K_chr16_HRC.r1-1_nomono_I4_cpra")
plink19 --recodeA --bfile GS20K_chr16_HRC.r1-1_nomono_I4_cpra --extract /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_pQTL_extraction_GS/chr_extraction_imputed/snplist.txt --out /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_pQTL_extraction_GS/chr_extraction_imputed/run_all_chr/chr_16_
chrs <- c("GS20K_chr17_HRC.r1-1_nomono_I4_cpra")
plink19 --recodeA --bfile GS20K_chr17_HRC.r1-1_nomono_I4_cpra --extract /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_pQTL_extraction_GS/chr_extraction_imputed/snplist.txt --out /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_pQTL_extraction_GS/chr_extraction_imputed/run_all_chr/chr_17_
chrs <- c("GS20K_chr18_HRC.r1-1_nomono_I4_cpra")
plink19 --recodeA --bfile GS20K_chr18_HRC.r1-1_nomono_I4_cpra --extract /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_pQTL_extraction_GS/chr_extraction_imputed/snplist.txt --out /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_pQTL_extraction_GS/chr_extraction_imputed/run_all_chr/chr_18_
chrs <- c("GS20K_chr19_HRC.r1-1_nomono_I4_cpra")
plink19 --recodeA --bfile GS20K_chr19_HRC.r1-1_nomono_I4_cpra --extract /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_pQTL_extraction_GS/chr_extraction_imputed/snplist.txt --out /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_pQTL_extraction_GS/chr_extraction_imputed/run_all_chr/chr_19_
chrs <- c("GS20K_chr20_HRC.r1-1_nomono_I4_cpra")
plink19 --recodeA --bfile GS20K_chr20_HRC.r1-1_nomono_I4_cpra --extract /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_pQTL_extraction_GS/chr_extraction_imputed/snplist.txt --out /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_pQTL_extraction_GS/chr_extraction_imputed/run_all_chr/chr_20_
chrs <- c("GS20K_chr21_HRC.r1-1_nomono_I4_cpra")
plink19 --recodeA --bfile GS20K_chr21_HRC.r1-1_nomono_I4_cpra --extract /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_pQTL_extraction_GS/chr_extraction_imputed/snplist.txt --out /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_pQTL_extraction_GS/chr_extraction_imputed/run_all_chr/chr_21_
chrs <- c("GS20K_chr22_HRC.r1-1_nomono_I4_cpra")
plink19 --recodeA --bfile GS20K_chr22_HRC.r1-1_nomono_I4_cpra --extract /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_pQTL_extraction_GS/chr_extraction_imputed/snplist.txt --out /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_pQTL_extraction_GS/chr_extraction_imputed/run_all_chr/chr_22_

## Get SNPs extracted for each chromosome and converted to .csv format 
setwd("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_pQTL_extraction_GS/chr_extraction_imputed/")
files <- list.files("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_pQTL_extraction_GS/chr_extraction_imputed/run_all_chr/", ".raw")

list_join <- list()

for (i in 1:length(files)){
	file <- as.character(files[i])
	SNP <- read.table(paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_pQTL_extraction_GS/chr_extraction_imputed/run_all_chr/", file), header = TRUE, stringsAsFactors = FALSE, check.names = F)
	head(SNP)[1:10]
	chr <- sub("chr_", "", file)
	chr <- sub("_.raw", "", chr)
	rownames(SNP) <- SNP$IID
	SNP <- SNP[-c(1:6)]
	colnames(SNP) <- gsub("_.$", "", colnames(SNP))
	write.csv(SNP, paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_pQTL_extraction_GS/chr_extraction_imputed/run_all_chr_converted_to_csv/", chr, ".csv"))
	list_join[[i]] <- SNP
	print(chr)
	print(dim(SNP))
}

setwd("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_pQTL_extraction_GS/chr_extraction_imputed/run_all_chr_converted_to_csv/")
dim(snplist) # 3960 in original SNP extraction list - doubled SNPs in list due to both allele combinations

# do a quick test to check row names are common across SNP files 
chr1 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_pQTL_extraction_GS/chr_extraction_imputed/run_all_chr_converted_to_csv/1.csv", check.names = F)
chr2 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_pQTL_extraction_GS/chr_extraction_imputed/run_all_chr_converted_to_csv/2.csv", check.names = F)
identical(row.names(chr1), row.names(chr2)) # TRUE
j1 <- cbind(chr1, chr2)

# bind the SNPs for all chr together and look at how many we have extracted
test <- do.call(cbind, list_join) 
dim(test) # 837

# Work out which SNPs are in the original extraction list but arent in the extracted SNP file 
which(order$new %in% colnames(test)) 
length(which(order2$new %in% colnames(test))) # 1690
length(which(snplist$new %in% colnames(test))) # 1690

# how many unique SNPs were in the list for extraction "order2" which was the correct format 
length(unique(order2$new)) # 1021

# so i potentially have 837 of a possible 1021 unique SNPs 
ov <- which(order2$new %in% colnames(test))
order3 <- order2[-ov,]
order4 <- order3[1]
order4 <- order4 %>% unique() # 184 SNPs not in the extraction (837 + 184 = 1021, our total unique SNPs in the list)
write.table(order4, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_pQTL_extraction_GS/chr_extraction_imputed/snplist_not_extracted.txt", quote = F, row.names = F, col.names = F)

# So of the 1021 SNPs in the sun et al list that are unique, there are 837 in the extraction and 184 (mostly indels) that arent extracted 
# Need to check whether its allele formatting thats hindering these, or whether it looks like there are no pQTLs present in GS for them 
# Match by chr_pos to see if i can pick up the alleles of interest here 

# Lets have a look at the ones for chr 22 
chr <- read.delim("/Cluster_Filespace/Marioni_Group/Daniel/GWAS_AgeAccel/Meta_Analysis_Daniel/cojo_files/GS_HRC/GS20K_chr22_HRC.r1-1_nomono_I4_cpra.bim", header = F, check.names = F)
which(chr$V4 %in% "17586715") # 0
which(chr$V4 %in% "33159919") # 0

# Chr 21
chr <- read.delim("/Cluster_Filespace/Marioni_Group/Daniel/GWAS_AgeAccel/Meta_Analysis_Daniel/cojo_files/GS_HRC/GS20K_chr21_HRC.r1-1_nomono_I4_cpra.bim", header = F, check.names = F)
which(chr$V4 %in% "22839703") # 0
which(chr$V4 %in% "47366802") # 0

chr <- read.delim("/Cluster_Filespace/Marioni_Group/Daniel/GWAS_AgeAccel/Meta_Analysis_Daniel/cojo_files/GS_HRC/GS20K_chr21_HRC.r1-1_nomono_I4_cpra.bim", header = F, check.names = F)
which(chr$V4 %in% "045946888") # 

# Chr 20
chr <- read.delim("/Cluster_Filespace/Marioni_Group/Daniel/GWAS_AgeAccel/Meta_Analysis_Daniel/cojo_files/GS_HRC/GS20K_chr20_HRC.r1-1_nomono_I4_cpra.bim", header = F, check.names = F)
which(chr$V4 %in% "2780762") # 0
which(chr$V4 %in% "3676788") # 0
which(chr$V4 %in% "24911330") # 0
which(chr$V4 %in% "29573649") # 0
which(chr$V4 %in% "32819223") # 0
which(chr$V4 %in% "62347189") # 0


# Chr 19
chr <- read.delim("/Cluster_Filespace/Marioni_Group/Daniel/GWAS_AgeAccel/Meta_Analysis_Daniel/cojo_files/GS_HRC/GS20K_chr20_HRC.r1-1_nomono_I4_cpra.bim", header = F, check.names = F)
which(chr$V4 %in% "836043") # 0
which(chr$V4 %in% "7781677") # 0
which(chr$V4 %in% "10190336") # 0
which(chr$V4 %in% "41933790") # 0
which(chr$V4 %in% "43679088") # 0
which(chr$V4 %in% "43696022") # 0
which(chr$V4 %in% "43711272") # 0
which(chr$V4 %in% "45413233") # 0
which(chr$V4 %in% "45424351") # 0
which(chr$V4 %in% "51222294") # 0
which(chr$V4 %in% "52004074") # 0
which(chr$V4 %in% "52130148") # 0
which(chr$V4 %in% "54337999") # 0
which(chr$V4 %in% "54338067") # 0
which(chr$V4 %in% "55249060") # 0
which(chr$V4 %in% "55374472") # 0


# It looks like we have a case in which the positional info doesnt match up for the missin SNPs either - so it is likely that GS does not have these present!
# Continue with the set of 837 SNPs for regressions
# Write this file out for use:

write.csv(test, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_pQTL_extraction_GS/chr_extraction_imputed/run_all_chr_joint_to_single_file/pQTLs_170521.csv", row.names = T)

