####################################################################################

### Calculate the epismoker measure for LBC1936 - original script 

####################################################################################

source("http://bioconductor.org/biocLite.R")
install.packages("devtools") # if you don't have the package, run install.packages("devtools")
library(devtools)
install_github("sailalithabollepalli/EpiSmokEr")
install.packages("IlluminaHumanMethylation450kmanifest")

library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(minfi)
library(htmlTable)
library(rmarkdown)

suppressPackageStartupMessages({
library(EpiSmokEr)  
library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(minfi)
library(htmlTable)
library(rmarkdown)
})

# Read in target file with matched IDs for methylation data 
load("/Cluster_Filespace/Marioni_Group/LBC/LBC_methylation/Beta_3525_norm_bgcorrect_0.001BetaThreshold_probefilter.RObject")
result2 <- epismoker(dataset=dat, method = "SSc")

# Save out data scores for smoking
write.rds("/Cluster_Filespace/Marioni_Group/Danni/lbc_epismoker.rds")


####################################################################################

### SORTING OUT THE INSTALLATION ISSUE 

####################################################################################

# # to normally install we would do:
# source("http://bioconductor.org/biocLite.R")
# install.packages("devtools") # if you don't have the package, run install.packages("devtools")
# library(devtools)
# install_github("sailalithabollepalli/EpiSmokEr")
# library(EpiSmokEr)

# Bioc manager update we now do:
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.12")
BiocManager::install(c("devtools"))

# Now we can install this package
BiocManager::install("IlluminaHumanMethylation450kmanifest")

# load all the other packages?
library(devtools)
library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(minfi)
library(htmlTable)
library(rmarkdown)

# Now install EpiSmokEr
install_github("sailalithabollepalli/EpiSmokEr")
library(EpiSmokEr)

# Read in target file with matched IDs for methylation data 
load("/Cluster_Filespace/Marioni_Group/LBC/LBC_methylation/Beta_3525_norm_bgcorrect_0.001BetaThreshold_probefilter.RObject")
result2 <- epismoker(dataset=dat, method = "SSc")

# Save out data scores for smoking
# write.rds("/Cluster_Filespace/Marioni_Group/Danni/lbc_epismoker.rds")

####################################################################################

### APPLY TO STRADL - the updated installed version

####################################################################################

# Taken from the EWAS models script in stradl_markers on 22_02_21 to get DNAm data 

## Read in STRADL methylation data and sample sheets 
W2 <- readRDS("/Cluster_Filespace/Marioni_Group/Rob/Somalogic/wave2-STRADL-mvals.rds")
W3 <- readRDS("/Cluster_Filespace/Marioni_Group/Rob/Somalogic/wave3-STRADL-mvals.rds")
sheet_w2 <- read.csv("/Cluster_Filespace/Marioni_Group/Rob/Somalogic/wave2_batch.csv")
sheet_w3 <- read.csv("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/wave3-final/samplesheet.final.csv")

# Read in target file 
target <- read.delim("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/Validation/STRADL_DNAm_target_REM_17April2020.txt")
IDs <- target[c(1,3)]

# Subset methylation so it contains only the 778 in the EWAS
ov <- which(colnames(W2) %in% IDs$DNAm_id)
wave2 <- W2[,ov]

ov2 <- which(colnames(W3) %in% IDs$DNAm_id)
wave3 <- W3[,ov2]

# > dim(wave3)
# [1] 773860    299
# > dim(wave2)
# [1] 793706    479


## Combine STRADL methylation data 
overlap <- which(rownames(wave3) %in% rownames(wave2))
length(overlap) # 772619

w3 <- wave3[overlap,]

overlap <- which(rownames(wave2) %in% rownames(w3))
length(overlap) # 772619

w2 <- wave2[overlap,]

# > dim(w2)
# [1] 772619    479
# > dim(w3)
# [1] 772619    299

identical(rownames(w2), rownames(w3)) # TRUE

join <- cbind(w2, w3)
dim(join) 


# Now apply EpiSmoker to calculate scores 
result3 <- epismoker(dataset=join, method = "SSc")

# Save out data scores for smoking
saveRDS(result3, "/Cluster_Filespace/Marioni_Group/Danni/stradl_epismoker.rds")

