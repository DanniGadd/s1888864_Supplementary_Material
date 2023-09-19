############################################################################################

### File transfer - preps

############################################################################################

# Save a test file in the desired format for transfer test run

screen

R

library(data.table)
library(tidyverse) 
library(readxl)

# annotations file
anno <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Annotations_protein_file/annotation_formatted_for_paper.csv")

# Set location
path <- "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/04_FULL/batches/batch1/"
setwd(paste0(path))
L <- list.files(".", ".mlma")

# Select one test file
j <- 1

# prep the file with seqid and gene name and save out as csv for ewas catalog
file <- fread("10000-28_MOA_774_with_GRM.mlma", header = T)
    file <- as.data.frame(file)
    name <- as.character(L[j])
    name <- gsub("_.*", "", name)
    file$SeqId <- name
    anno_sub <- anno[which(anno$SeqId %in% name),]
    protein_gene <- as.character(anno_sub[,3])
    file$Protein_gene <- protein_gene
    file <- file[order(file$p),]
write.csv(file, paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_File_transfer/Files/", "MWAS_", name, "_", protein_gene, ".csv"), row.names = F)

# Location for file transfer is therefore: 
# /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_File_transfer/Files/


###########################################################################################

# Save out full files for all proteins ready to transfer

# Batch 1 

screen

R

library(data.table)
library(tidyverse) 
library(readxl)

# annotations file
anno <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Annotations_protein_file/annotation_formatted_for_paper.csv")

# Set location
path <- "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/04_FULL/batches/batch1/"
setwd(paste0(path))
L <- list.files(".", ".mlma")

# Select one test file
for (j in 1:length(L)){
file <- fread(L[j], header = T)
    file <- as.data.frame(file)
    name <- as.character(L[j])
    name <- gsub("_.*", "", name)
    file$SeqId <- name
    anno_sub <- anno[which(anno$SeqId %in% name),]
    protein_gene <- as.character(anno_sub[,3])
    file$Protein_gene <- protein_gene
    file <- file[order(file$p),]
write.csv(file, paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_File_transfer/Fully_adjusted_model_files/", "MWAS_", name, "_", protein_gene, ".csv"), row.names = F)
}


# Batch 2

screen

R

library(data.table)
library(tidyverse) 
library(readxl)

# annotations file
anno <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Annotations_protein_file/annotation_formatted_for_paper.csv")

# Set location
path <- "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/04_FULL/batches/batch2/"
setwd(paste0(path))
L <- list.files(".", ".mlma")

# Select one test file
for (j in 1:length(L)){
file <- fread(L[j], header = T)
    file <- as.data.frame(file)
    name <- as.character(L[j])
    name <- gsub("_.*", "", name)
    file$SeqId <- name
    anno_sub <- anno[which(anno$SeqId %in% name),]
    protein_gene <- as.character(anno_sub[,3])
    file$Protein_gene <- protein_gene
    file <- file[order(file$p),]
write.csv(file, paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_File_transfer/Fully_adjusted_model_files/", "MWAS_", name, "_", protein_gene, ".csv"), row.names = F)
}



# Batch 3

screen

R

library(data.table)
library(tidyverse) 
library(readxl)

# annotations file
anno <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Annotations_protein_file/annotation_formatted_for_paper.csv")

# Set location
path <- "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/04_FULL/batches/batch3/"
setwd(paste0(path))
L <- list.files(".", ".mlma")

# Select one test file
for (j in 1:length(L)){
file <- fread(L[j], header = T)
    file <- as.data.frame(file)
    name <- as.character(L[j])
    name <- gsub("_.*", "", name)
    file$SeqId <- name
    anno_sub <- anno[which(anno$SeqId %in% name),]
    protein_gene <- as.character(anno_sub[,3])
    file$Protein_gene <- protein_gene
    file <- file[order(file$p),]
write.csv(file, paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_File_transfer/Fully_adjusted_model_files/", "MWAS_", name, "_", protein_gene, ".csv"), row.names = F)
}




# Batch 4

screen

R

library(data.table)
library(tidyverse) 
library(readxl)

# annotations file
anno <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Annotations_protein_file/annotation_formatted_for_paper.csv")

# Set location
path <- "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/04_FULL/batches/batch4/"
setwd(paste0(path))
L <- list.files(".", ".mlma")

# Select one test file
for (j in 1:length(L)){
file <- fread(L[j], header = T)
    file <- as.data.frame(file)
    name <- as.character(L[j])
    name <- gsub("_.*", "", name)
    file$SeqId <- name
    anno_sub <- anno[which(anno$SeqId %in% name),]
    protein_gene <- as.character(anno_sub[,3])
    file$Protein_gene <- protein_gene
    file <- file[order(file$p),]
write.csv(file, paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_File_transfer/Fully_adjusted_model_files/", "MWAS_", name, "_", protein_gene, ".csv"), row.names = F)
}



# Batch 5

screen

R

library(data.table)
library(tidyverse) 
library(readxl)

# annotations file
anno <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Annotations_protein_file/annotation_formatted_for_paper.csv")

# Set location
path <- "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/04_FULL/batches/batch5/"
setwd(paste0(path))
L <- list.files(".", ".mlma")

# Select one test file
for (j in 1:length(L)){
file <- fread(L[j], header = T)
    file <- as.data.frame(file)
    name <- as.character(L[j])
    name <- gsub("_.*", "", name)
    file$SeqId <- name
    anno_sub <- anno[which(anno$SeqId %in% name),]
    protein_gene <- as.character(anno_sub[,3])
    file$Protein_gene <- protein_gene
    file <- file[order(file$p),]
write.csv(file, paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_File_transfer/Fully_adjusted_model_files/", "MWAS_", name, "_", protein_gene, ".csv"), row.names = F)
}




# Batch 6

screen

R

library(data.table)
library(tidyverse) 
library(readxl)

# annotations file
anno <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Annotations_protein_file/annotation_formatted_for_paper.csv")

# Set location
path <- "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/04_FULL/batches/batch6/"
setwd(paste0(path))
L <- list.files(".", ".mlma")

# Select one test file
for (j in 1:length(L)){
file <- fread(L[j], header = T)
    file <- as.data.frame(file)
    name <- as.character(L[j])
    name <- gsub("_.*", "", name)
    file$SeqId <- name
    anno_sub <- anno[which(anno$SeqId %in% name),]
    protein_gene <- as.character(anno_sub[,3])
    file$Protein_gene <- protein_gene
    file <- file[order(file$p),]
write.csv(file, paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_File_transfer/Fully_adjusted_model_files/", "MWAS_", name, "_", protein_gene, ".csv"), row.names = F)
}


# Batch 7

screen

R

library(data.table)
library(tidyverse) 
library(readxl)

# annotations file
anno <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Annotations_protein_file/annotation_formatted_for_paper.csv")

# Set location
path <- "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/04_FULL/batches/batch7/"
setwd(paste0(path))
L <- list.files(".", ".mlma")

# Select one test file
for (j in 1:length(L)){
file <- fread(L[j], header = T)
    file <- as.data.frame(file)
    name <- as.character(L[j])
    name <- gsub("_.*", "", name)
    file$SeqId <- name
    anno_sub <- anno[which(anno$SeqId %in% name),]
    protein_gene <- as.character(anno_sub[,3])
    file$Protein_gene <- protein_gene
    file <- file[order(file$p),]
write.csv(file, paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_File_transfer/Fully_adjusted_model_files/", "MWAS_", name, "_", protein_gene, ".csv"), row.names = F)
}



# Batch 8

screen

R

library(data.table)
library(tidyverse) 
library(readxl)

# annotations file
anno <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Annotations_protein_file/annotation_formatted_for_paper.csv")

# Set location
path <- "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/04_FULL/batches/batch8/"
setwd(paste0(path))
L <- list.files(".", ".mlma")

# Select one test file
for (j in 1:length(L)){
file <- fread(L[j], header = T)
    file <- as.data.frame(file)
    name <- as.character(L[j])
    name <- gsub("_.*", "", name)
    file$SeqId <- name
    anno_sub <- anno[which(anno$SeqId %in% name),]
    protein_gene <- as.character(anno_sub[,3])
    file$Protein_gene <- protein_gene
    file <- file[order(file$p),]
write.csv(file, paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_File_transfer/Fully_adjusted_model_files/", "MWAS_", name, "_", protein_gene, ".csv"), row.names = F)
}



# Batch 9

screen

R

library(data.table)
library(tidyverse) 
library(readxl)

# annotations file
anno <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Annotations_protein_file/annotation_formatted_for_paper.csv")

# Set location
path <- "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/04_FULL/batches/batch9/"
setwd(paste0(path))
L <- list.files(".", ".mlma")

# Select one test file
for (j in 1:length(L)){
file <- fread(L[j], header = T)
    file <- as.data.frame(file)
    name <- as.character(L[j])
    name <- gsub("_.*", "", name)
    file$SeqId <- name
    anno_sub <- anno[which(anno$SeqId %in% name),]
    protein_gene <- as.character(anno_sub[,3])
    file$Protein_gene <- protein_gene
    file <- file[order(file$p),]
write.csv(file, paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_File_transfer/Fully_adjusted_model_files/", "MWAS_", name, "_", protein_gene, ".csv"), row.names = F)
}



# Batch 10

screen

R

library(data.table)
library(tidyverse) 
library(readxl)

# annotations file
anno <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Annotations_protein_file/annotation_formatted_for_paper.csv")

# Set location
path <- "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/04_FULL/batches/batch10/"
setwd(paste0(path))
L <- list.files(".", ".mlma")

# Select one test file
for (j in 1:length(L)){
file <- fread(L[j], header = T)
    file <- as.data.frame(file)
    name <- as.character(L[j])
    name <- gsub("_.*", "", name)
    file$SeqId <- name
    anno_sub <- anno[which(anno$SeqId %in% name),]
    protein_gene <- as.character(anno_sub[,3])
    file$Protein_gene <- protein_gene
    file <- file[order(file$p),]
write.csv(file, paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_File_transfer/Fully_adjusted_model_files/", "MWAS_", name, "_", protein_gene, ".csv"), row.names = F)
}



# Batch 11

screen

R

library(data.table)
library(tidyverse) 
library(readxl)

# annotations file
anno <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Annotations_protein_file/annotation_formatted_for_paper.csv")

# Set location
path <- "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/04_FULL/batches/batch11/"
setwd(paste0(path))
L <- list.files(".", ".mlma")

# Select one test file
for (j in 1:length(L)){
file <- fread(L[j], header = T)
    file <- as.data.frame(file)
    name <- as.character(L[j])
    name <- gsub("_.*", "", name)
    file$SeqId <- name
    anno_sub <- anno[which(anno$SeqId %in% name),]
    protein_gene <- as.character(anno_sub[,3])
    file$Protein_gene <- protein_gene
    file <- file[order(file$p),]
write.csv(file, paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_File_transfer/Fully_adjusted_model_files/", "MWAS_", name, "_", protein_gene, ".csv"), row.names = F)
}




# Batch 12

screen

R

library(data.table)
library(tidyverse) 
library(readxl)

# annotations file
anno <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Annotations_protein_file/annotation_formatted_for_paper.csv")

# Set location
path <- "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/04_FULL/batches/batch12/"
setwd(paste0(path))
L <- list.files(".", ".mlma")

# Select one test file
for (j in 1:length(L)){
file <- fread(L[j], header = T)
    file <- as.data.frame(file)
    name <- as.character(L[j])
    name <- gsub("_.*", "", name)
    file$SeqId <- name
    anno_sub <- anno[which(anno$SeqId %in% name),]
    protein_gene <- as.character(anno_sub[,3])
    file$Protein_gene <- protein_gene
    file <- file[order(file$p),]
write.csv(file, paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_File_transfer/Fully_adjusted_model_files/", "MWAS_", name, "_", protein_gene, ".csv"), row.names = F)
}


# Check all relevant files are in the desired location

F <- list.files("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_File_transfer/Fully_adjusted_model_files/")

# Check that the SeqIds match those listed in the annotations file

anno <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Annotations_protein_file/annotation_formatted_for_paper.csv")

anno_sub <- anno[which(anno$SeqId %in% "9999-1"),] # read in okay - IRF6 is '9999-1'

names <- list()

for (i in 1:length(F)){
    name <- as.character(F[i])
    name <- gsub("MWAS_", "", name)
    name <- gsub("_.*", "", name)
    names[i] <- name
}

name_list <- do.call(rbind, names)
name_list <- as.data.frame(name_list)
names(name_list)[1] <- "SeqId"

anno2 <- anno[which(anno$SeqId %in% name_list$SeqId),]
dim(anno2) # all 4231 SeqIds are in the annotation file 

not <- anno[-which(anno$SeqId %in% name_list$SeqId),]