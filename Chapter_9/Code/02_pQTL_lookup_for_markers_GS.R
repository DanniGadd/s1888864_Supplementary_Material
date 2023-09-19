####################################################################################

### pQTL extraction for each proteins via GWAS results 

####################################################################################

# In this script the GWAS summary statistics are loaded for carolines GWAS runs for each protein
# We then save out COJO-formatted files and perform COJO on the linear regression statistics 
# Sentinel SNPs for GDF15 and NtproBNP are then extracted from the imputed GS dataset 
# These will be used to regress out the genetic effects of the protein when generating EpiScores 

cd /Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/GWAS_caroline/

screen

R

library(readxl)
library(tidyverse)
library(data.table)

# Load in the GWAS sum stats available for each protein 
GDF <- fread("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/GWAS_caroline/GS20K_gdf15_rnk_My8MHH5F8V_imp.stats.gz")
BNP <- fread("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/GWAS_caroline/GS20K_nt.probnp_rnk_1ePbnCoIl6_imp.stats.gz")

#see most significant SNP in each case 
GDF[which.min(GDF$P_BOLT_LMM),]
BNP[which.min(BNP$P_BOLT_LMM),]

# > GDF[which.min(GDF$P_BOLT_LMM),]
#                SNP CHR       BP  GENPOS ALLELE1 ALLELE0   A1FREQ     INFO
# 1: 19_18500722_G_A  19 18500722 0.45168       G       A 0.739515 0.936485
#    CHISQ_LINREG P_LINREG      BETA        SE CHISQ_BOLT_LMM_INF P_BOLT_LMM_INF
# 1:       1111.1 1.3e-243 -0.338444 0.0105398            1031.11       3.1e-226
#    CHISQ_BOLT_LMM P_BOLT_LMM
# 1:        1048.24   5.9e-230
# > BNP[which.min(BNP$P_BOLT_LMM),]
#               SNP CHR       BP   GENPOS ALLELE1 ALLELE0   A1FREQ    INFO
# 1: 1_11919271_A_G   1 11919271 0.239924       A       G 0.577163 0.99473
#    CHISQ_LINREG P_LINREG      BETA         SE CHISQ_BOLT_LMM_INF P_BOLT_LMM_INF
# 1:      417.552  8.3e-93 -0.198019 0.00978029            409.931        3.8e-91
#    CHISQ_BOLT_LMM P_BOLT_LMM
# 1:        408.752    6.9e-91

# Check to see which P caroline filtered main results from for each protein 
GDF1 <- GDF[GDF$P_BOLT_LMM < 5e-8,]
table(GDF1$CHR) # 314
gdfhits <- read.delim("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/GWAS_caroline/GS20K_gdf15_rnk_My8MHH5F8V_tophits.tsv")
bnphits <- read.delim("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/GWAS_caroline/GS20K_nt.probnp_rnk_1ePbnCoIl6_tophits.tsv")

# > dim(bnphits)
# [1] 911   9
# > dim(gdfhits)
# [1] 314   9

# Caroline filtered her results by the LMM P value, however there is no beta/se for these summary statistics so ill look to filter by linear regression now
GDF2 <- GDF[GDF$P_LINREG < 5e-8,]
BNP2 <- BNP[BNP$P_LINREG < 5e-8,]

# > dim(GDF2)
# [1] 330  18
# > dim(BNP2)
# [1] 943  18

which(GDF2$BP %in% "18388612")

# As these are not too dissimilar, should we therefore use the linear results to generate sentinel SNP estimates?

# Check the chromosomes required for each protein if so
table(gdfhits$CHR) # chr 19 - 314 SNPs
table(bnphits$CHR) # chr 1,4,8,12,18 (SNP counts below)
#   1   4   8  12  18
# 511 120  42 237   1

table(GDF2$CHR) # chr 19 - 330 SNPs
table(BNP2$CHR) # chr 1,4,8,12,18 (SNP counts below)
#   1   4   8  12  18
# 517 125  59 241   1


# For now, ill also put the maximum N for each protein as N - confirm with Caroline 
GDF$N <- 18414
BNP$N <- 17863

# Subset and rename correct columns - using linear regression P as it has accompanying beta/se values 
GDF_cojo <- GDF[,c("SNP", "ALLELE1", "ALLELE0", "A1FREQ", "BETA", "SE", "P_LINREG", "N")]
BNP_cojo <- BNP[,c("SNP", "ALLELE1", "ALLELE0", "A1FREQ", "BETA", "SE", "P_LINREG", "N")]

names(GDF_cojo) <- c("SNP", "A1", "A2", "freq", "b", "se", "p", "N")
names(BNP_cojo) <- c("SNP", "A1", "A2", "freq", "b", "se", "p", "N")

write.table(GDF_cojo, "/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/GWAS_caroline/GDF_COJO_table.txt", row.names = F, quote = F)
write.table(BNP_cojo, "/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/GWAS_caroline/BNP_COJO_table.txt", row.names = F, quote = F)

####################################################################################

### COJO - RUN FOR EACH CHROMOSOME OF INTEREST USING BIM/BED/FAM DATA FROM GS 

## Take a look at imputed data format 
head /Cluster_Filespace/Marioni_Group/Daniel/GWAS_AgeAccel/Meta_Analysis_Daniel/cojo_files/GS_HRC/GS20K_chr10_HRC.r1-1_nomono_I4_cpra.bim

## GDF15 

# COJO step (chr 19)
gcta64  --bfile /Cluster_Filespace/Marioni_Group/Daniel/GWAS_AgeAccel/Meta_Analysis_Daniel/cojo_files/GS_HRC/GS20K_chr19_HRC.r1-1_nomono_I4_cpra  --chr 19 --cojo-file /Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/GWAS_caroline/GDF_COJO_table.txt --cojo-slct --cojo-actual-geno  --out /Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/GWAS_caroline/GDF_cojo_test_file

# Take a look at results files 
a <- read.delim("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/GWAS_caroline/GDF_cojo_test_file.jma.cojo")

# SNPs selected for GDF15 chr 19 
#   Chr             SNP       bp refA     freq         b        se            p
# 1  19 19_18500722_G_A 18500722    G 0.739515 -0.338444 0.0105398 3.09146e-226
# 2  19 19_18501034_C_T 18501034    C 0.834774 -0.236603 0.0121434  1.49843e-84
# 3  19 19_18501035_G_A 18501035    G 0.772361  0.092834 0.0110012  3.21343e-17
#       n freq_geno        bJ     bJ_se           pJ      LD_r
# 1 18414  0.754543 -0.465062 0.0110644  0.00000e+00 -0.247411
# 2 18414  0.841079 -0.409408 0.0127872 6.33253e-225 -0.235338
# 3 18414  0.779403 -0.138766 0.0114419  7.51309e-34  0.000000

## BNP - chr 1,4,8,12,18

# COJO step (chr 1)
gcta64  --bfile /Cluster_Filespace/Marioni_Group/Daniel/GWAS_AgeAccel/Meta_Analysis_Daniel/cojo_files/GS_HRC/GS20K_chr1_HRC.r1-1_nomono_I4_cpra  --chr 1 --cojo-file /Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/GWAS_caroline/BNP_COJO_table.txt --cojo-slct --cojo-actual-geno  --out /Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/GWAS_caroline/BNPchr1_cojo_test_file

t <- read.delim("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/GWAS_caroline/BNPchr1_cojo_test_file.jma.cojo")

#   Chr            SNP       bp refA     freq          b         se           p
# 1   1 1_11905981_A_G 11905981    A 0.846378  0.0148017 0.01321960 2.62851e-01
# 2   1 1_11908146_C_G 11908146    C 0.987356 -0.4447220 0.04712550 3.83709e-21
# 3   1 1_11919271_A_G 11919271    A 0.577163 -0.1980190 0.00978029 3.79532e-91
#       n freq_geno        bJ     bJ_se           pJ       LD_r
# 1 17863  0.847070  0.129488 0.0143263  1.58828e-19 -0.0435945
# 2 17863  0.989642 -0.289431 0.0474113  1.03000e-09  0.1233390
# 3 17863  0.584839 -0.227221 0.0106300 2.26181e-101  0.0000000

# COJO step (chr 4)
gcta64  --bfile /Cluster_Filespace/Marioni_Group/Daniel/GWAS_AgeAccel/Meta_Analysis_Daniel/cojo_files/GS_HRC/GS20K_chr4_HRC.r1-1_nomono_I4_cpra  --chr 4 --cojo-file /Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/GWAS_caroline/BNP_COJO_table.txt --cojo-slct --cojo-actual-geno  --out /Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/GWAS_caroline/BNPchr4_cojo_test_file

d <- read.delim("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/GWAS_caroline/BNPchr4_cojo_test_file.jma.cojo")

#   Chr             SNP        bp refA     freq          b        se           p
# 1   4   4_9209010_C_T   9209010    C 0.545908 -0.0553814 0.0146905 1.63327e-04
# 2   4  4_73458197_T_A  73458197    T 0.507910  0.0575759 0.0133329 1.57212e-05
# 3   4 4_103188709_C_T 103188709    C 0.918300 -0.1973130 0.0175578 2.65592e-29
#       n freq_geno         bJ      bJ_se          pJ        LD_r
# 1 17863  0.535668 -0.0560828 0.00959350 5.03787e-09 -0.00165795
# 2 17863  0.511781  0.0569261 0.00955536 2.56111e-09 -0.00524643
# 3 17863  0.921051 -0.1975330 0.01777450 1.08121e-28  0.00000000

# COJO step (chr 8)
gcta64  --bfile /Cluster_Filespace/Marioni_Group/Daniel/GWAS_AgeAccel/Meta_Analysis_Daniel/cojo_files/GS_HRC/GS20K_chr8_HRC.r1-1_nomono_I4_cpra  --chr 8 --cojo-file /Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/GWAS_caroline/BNP_COJO_table.txt --cojo-slct --cojo-actual-geno  --out /Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/GWAS_caroline/BNPchr8_cojo_test_file

j <- read.delim("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/GWAS_caroline/BNPchr8_cojo_test_file.jma.cojo")


#   Chr            SNP       bp refA     freq         b        se           p
# 1   8 8_22266530_C_T 22266530    C 0.698742 0.0608844 0.0104014 4.81389e-09
#       n freq_geno        bJ     bJ_se          pJ LD_r
# 1 17863  0.697958 0.0608844 0.0104693 6.04377e-09    0


# COJO step (chr 12)
gcta64  --bfile /Cluster_Filespace/Marioni_Group/Daniel/GWAS_AgeAccel/Meta_Analysis_Daniel/cojo_files/GS_HRC/GS20K_chr12_HRC.r1-1_nomono_I4_cpra  --chr 12 --cojo-file /Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/GWAS_caroline/BNP_COJO_table.txt --cojo-slct --cojo-actual-geno  --out /Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/GWAS_caroline/BNPchr12_cojo_test_file

w <- read.delim("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/GWAS_caroline/BNPchr12_cojo_test_file.jma.cojo")

#   Chr             SNP       bp refA     freq         b        se           p
# 1  12 12_89829728_G_T 89829728    G 0.882861 0.0714746 0.0151624 2.42976e-06
# 2  12 12_89839072_C_T 89839072    C 0.769347 0.1389930 0.0113512 1.79152e-34
#       n freq_geno       bJ     bJ_se          pJ      LD_r
# 1 17863  0.885808 0.111931 0.0153493 3.04933e-13 -0.196184
# 2 17863  0.766723 0.155558 0.0115788 3.78371e-41  0.000000


# COJO step (chr 18)
gcta64  --bfile /Cluster_Filespace/Marioni_Group/Daniel/GWAS_AgeAccel/Meta_Analysis_Daniel/cojo_files/GS_HRC/GS20K_chr18_HRC.r1-1_nomono_I4_cpra  --chr 18 --cojo-file /Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/GWAS_caroline/BNP_COJO_table.txt --cojo-slct --cojo-actual-geno  --out /Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/GWAS_caroline/BNPchr18_cojo_test_file

h <- read.delim("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/GWAS_caroline/BNPchr18_cojo_test_file.jma.cojo")


#   Chr             SNP       bp refA     freq         b        se           p
# 1  18 18_73681727_G_A 73681727    G 0.969302 -0.169896 0.0297451 1.11833e-08
#       n freq_geno        bJ    bJ_se         pJ LD_r
# 1 17863  0.972818 -0.169896 0.029733 1.1031e-08    0

# Make a joint table with all hits that are sentinel 
a$Protein <- "GDF15"
t$Protein <- "NtproBNP"
d$Protein <- "NtproBNP"
j$Protein <- "NtproBNP"
w$Protein <- "NtproBNP"
h$Protein <- "NtproBNP"

joint <- rbind(a,t)
joint <- rbind(joint,d)
joint <- rbind(joint,j)
joint <- rbind(joint,w)
joint <- rbind(joint,h)

# Save out table as index for regressions in episcore generation
write.csv(joint, "/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/GWAS_caroline/pQTL_extraction/pQTL_extraction_index_table.csv", row.names = F)

####################################################################################

### AFTER pQTLs IDENTIFIED, PULL THEM OUT FROM GS GENETIC DATASETS TO ADJUST FOR IN SCORE GENERATION

# Create a list of the pQTLs implicated from cojo for each protein and save out in specific format for input 
GDFlist <- data.frame(X = c("19_18500722_G_A", "19_18501034_C_T", "19_18501035_G_A"))
write.table(GDFlist, "/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/GWAS_caroline/GDF15_snplist.txt", quote = F, row.names = F, col.names = F)

BNPlist <- data.frame(X = c("1_11905981_A_G", "1_11908146_C_G", "1_11919271_A_G", "4_9209010_C_T", "4_73458197_T_A", "4_103188709_C_T",
"8_22266530_C_T", "12_89829728_G_T", "12_89839072_C_T", "18_73681727_G_A"))
write.table(BNPlist, "/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/GWAS_caroline/BNP_snplist.txt", quote = F, row.names = F, col.names = F)

# Extract available SNPs from the GS imputed data 

# GDF15 - chr 19 
cd /Cluster_Filespace/Marioni_Group/Daniel/GWAS_AgeAccel/Meta_Analysis_Daniel/cojo_files/GS_HRC/
plink19 --recodeA --bfile GS20K_chr19_HRC.r1-1_nomono_I4_cpra --extract /Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/GWAS_caroline/GDF15_snplist.txt --out /Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/GWAS_caroline/pQTL_extraction/GDF15_chr_19_

# BNP - chr 1,4,8,12,18
cd /Cluster_Filespace/Marioni_Group/Daniel/GWAS_AgeAccel/Meta_Analysis_Daniel/cojo_files/GS_HRC/
plink19 --recodeA --bfile GS20K_chr1_HRC.r1-1_nomono_I4_cpra --extract /Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/GWAS_caroline/BNP_snplist.txt --out /Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/GWAS_caroline/pQTL_extraction/BNP_chr_1_

cd /Cluster_Filespace/Marioni_Group/Daniel/GWAS_AgeAccel/Meta_Analysis_Daniel/cojo_files/GS_HRC/
plink19 --recodeA --bfile GS20K_chr4_HRC.r1-1_nomono_I4_cpra --extract /Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/GWAS_caroline/BNP_snplist.txt --out /Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/GWAS_caroline/pQTL_extraction/BNP_chr_4_

cd /Cluster_Filespace/Marioni_Group/Daniel/GWAS_AgeAccel/Meta_Analysis_Daniel/cojo_files/GS_HRC/
plink19 --recodeA --bfile GS20K_chr8_HRC.r1-1_nomono_I4_cpra --extract /Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/GWAS_caroline/BNP_snplist.txt --out /Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/GWAS_caroline/pQTL_extraction/BNP_chr_8_

cd /Cluster_Filespace/Marioni_Group/Daniel/GWAS_AgeAccel/Meta_Analysis_Daniel/cojo_files/GS_HRC/
plink19 --recodeA --bfile GS20K_chr12_HRC.r1-1_nomono_I4_cpra --extract /Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/GWAS_caroline/BNP_snplist.txt --out /Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/GWAS_caroline/pQTL_extraction/BNP_chr_12_

cd /Cluster_Filespace/Marioni_Group/Daniel/GWAS_AgeAccel/Meta_Analysis_Daniel/cojo_files/GS_HRC/
plink19 --recodeA --bfile GS20K_chr18_HRC.r1-1_nomono_I4_cpra --extract /Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/GWAS_caroline/BNP_snplist.txt --out /Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/GWAS_caroline/pQTL_extraction/BNP_chr_18_


## Get SNPs extracted for each chromosome and converted to .csv format 
setwd("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/GWAS_caroline/pQTL_extraction/")

## GDF extraction

files <- list.files("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/GWAS_caroline/pQTL_extraction/", ".raw")

list_join <- list()

i <- 6
  file <- as.character(files[i])
  SNP <- read.table(paste0("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/GWAS_caroline/pQTL_extraction/", file), header = TRUE, stringsAsFactors = FALSE, check.names = F)
  head(SNP)
  chr <- sub("BNP_chr_", "", file)
  chr <- sub("_.raw", "", chr)
  rownames(SNP) <- SNP$IID
  SNP <- SNP[-c(1:6)]
  colnames(SNP) <- gsub("_.$", "", colnames(SNP))
  write.csv(SNP, paste0("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/GWAS_caroline/pQTL_extraction/GDF_formatted/", chr, ".csv"))

## BNP extraction

for (i in 1:5){
  file <- as.character(files[i])
  SNP <- read.table(paste0("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/GWAS_caroline/pQTL_extraction/", file), header = TRUE, stringsAsFactors = FALSE, check.names = F)
  head(SNP)
  chr <- sub("BNP_chr_", "", file)
  chr <- sub("_.raw", "", chr)
  rownames(SNP) <- SNP$IID
  SNP <- SNP[-c(1:6)]
  colnames(SNP) <- gsub("_.$", "", colnames(SNP))
  write.csv(SNP, paste0("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/GWAS_caroline/pQTL_extraction/BNP_formatted/", chr, ".csv"))
  list_join[[i]] <- SNP
  print(chr)
  print(dim(SNP))
}

# [1] "1"
# [1] 20032     3
# [1] "12"
# [1] 20032     2
# [1] "18"
# [1] 20032     1
# [1] "4"
# [1] 20032     3
# [1] "8"
# [1] 20032     1 # all 10 extracted 

# bind the SNPs for all chr together and look at how many we have extracted
test <- do.call(cbind, list_join) 
dim(test) # 20032    10
write.csv(test, "/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/GWAS_caroline/pQTL_extraction/BNP_formatted/all_BNP.csv")

# Check that the files with SNPs read in okay for use 
BNP_ex <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/GWAS_caroline/pQTL_extraction/BNP_formatted/all_BNP.csv", check.names = F)
GDF_ex <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/GWAS_caroline/pQTL_extraction/GDF_formatted/GDF15_chr_19.csv", check.names = F)


####################################################################################

