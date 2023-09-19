############################################################################################

### Lookup of cis trans effects for MWAS 

############################################################################################

cd /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/cis_trans/

screen

R

library(tidyverse)
library(readxl)

# Annotations
anno <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Annotations_protein_file/annotation_formatted_for_paper.csv")

# Load in the significant MWAS results with proteins listed to get unique proteins
basic <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Results_collation/Basic_EWAS_results_second_threshold_020221.csv")
WBC <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Results_collation/WBC_EWAS_results_second_threshold_020221.csv")
full <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Results_collation/FULL_EWAS_results_second_threshold_020221.csv")

###############################

## FULL MODEL FIRST 

cpgs <- full
unique <- unique(cpgs$SeqId)

length(unique)
# [1] 195

unique <- cpgs[c(1,10,11)] %>% unique() # xtracted into dataset with gene names 
tab <- unique
names(tab)[3] <- "gene" 

# get gene info set up from external source (Use Ensembl to Extract Relevant Features)
library(biomaRt)
ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
attributes = listAttributes(ensembl, page = "structure")
attributes[grep("transcript", attributes$description, ignore.case = TRUE), ]

# get list of proteins gene names 
list <- tab$gene

## Get Transcription Start Sites for soma proteins 
tss <- getBM(attributes = c("transcription_start_site", "chromosome_name",
                            "transcript_start", "transcript_end",
                            "strand",  "ensembl_gene_id",
                            "ensembl_transcript_id", "external_gene_name"),
             filters = "external_gene_name", "UCSC_RefGene_Group", values = c(list),
             mart = ensembl)

# see how many this worked for
length(unique(tss$external_gene_name)) 

# See if there are unique TSS extracted for these proteins 
dim(tss) 
length(unique(tss$transcription_start_site)) 

# which are the gene names that didnt extract
ov <- which(!list %in% tss$external_gene_name)
list2 <- list[ov]

# [1] "HDGFRP3" "CECR1"   "NT5C3L"  "CASC4"   "ADSSL1"

# HDGFL3
# ADA2 
# NT5C3B
# GOLM2
# ADSS1

# Try these aliases for the genes as per lookup on gene cards site 

list3 <- c("HDGFL3","NT5C3B","ADA2", "GOLM2", "ADSS1")

## Get Transcription Start Sites for remaining soma proteins 
tss2 <- getBM(attributes = c("transcription_start_site", "chromosome_name",
                            "transcript_start", "transcript_end",
                            "strand",  "ensembl_gene_id",
                            "ensembl_transcript_id", "external_gene_name"),
             filters = "external_gene_name", values = c(list3),
             mart = ensembl)

# see how many this worked for
length(unique(tss2$external_gene_name)) # 5

# Now we have a complete extracted set of tss information - join them together 

tss <- rbind(tss, tss2)

## set up code to create biomart dataframe with gene, start and end columns
biomart_df <- tss

names(biomart_df)[8] <- "gene"
names(biomart_df)[3] <- "start"
names(biomart_df)[4] <- "end"
names(biomart_df)[2] <- "chromosome_ensembl"

# Again, lets do a check to see if each protein has its own TSS extracted that is unique 

subset <- biomart_df[,c("transcription_start_site", "gene")]
subset <- unique(subset$gene)

## separate data into genes on positive strand and genes on negative strand

biomart_df_pos <- biomart_df[biomart_df$strand == 1,]
biomart_df_neg <- biomart_df[biomart_df$strand == -1,]

# sanity check to see if overlap is present/correct number of genes are present

a = unique(biomart_df_pos$gene)
b = unique(biomart_df_neg$gene)
which(a %in% b)
length(a)
length(b)

## create output data frames for genes on positive strand and those on negative strand to be r-bound later
   
   
    ##out_df1 == genes on positive strand
   
out_df1 <- matrix(nrow=length(unique(biomart_df_pos$gene)), ncol=2)
colnames(out_df1) <- c("start", "end")
rownames(out_df1) <- unique(biomart_df_pos$gene)


    ##out_df2 == genes on negative strand
   
out_df2 <- matrix(nrow=length(unique(biomart_df_neg$gene)), ncol=2)
colnames(out_df2) <- c("start", "end")
rownames(out_df2) <- unique(biomart_df_neg$gene)

## loop through gene names to fill in min 5'end (start) and max 3'end (end) as on positive strand you want furthest away 5'(tss) and biggest distance to 3' end
## fill out out_df1 based on the above

for(gene in unique(biomart_df_pos$gene)) {
tmp <- biomart_df_pos[which(biomart_df_pos$gene==gene), ]
out_df1[gene,"start"] <-  min(tmp$start)
out_df1[gene,"end"] <-  max(tmp$end)
}

## loop through gene names to fill in min 3'end (end) and max 5'end (start) as on negative strand you want furthest away 3'(tss (relative to positive strand)) and biggest distance to 5' end
## fill out out_df2 based on the above

for(gene in unique(biomart_df_neg$gene)) {
tmp <- biomart_df_neg[which(biomart_df_neg$gene==gene), ]
out_df2[gene,"start"] <-  max(tmp$end)
out_df2[gene,"end"] <-  min(tmp$start)
}

## row bind the positive and negative strand dataframes to give one out_df
out_df <- rbind(out_df1, out_df2)


## Tidy the files before saving for new row with protein names (out_df), then one for length of transcripts in biomart_df for future reference

    ## set as data frames
    out_df <- as.data.frame(out_df)
    biomart_df <- as.data.frame(biomart_df)

out_df$Gene <- row.names(out_df)
biomart_df$transcript_length <- abs(biomart_df$start - biomart_df$end)

# Check unique start sites here 
unique(out_df$start) # 143
# plot(out_df$start)

## Get Chromosome Number for Protein 
biomart_df$chromosome_ensembl2 <-  gsub("CHR_HSCHR", "", biomart_df$chromosome_ensembl)
biomart_df$chromosome_ensembl2 <- gsub("_.*", "", biomart_df$chromosome_ensembl2)
biomart_df <- biomart_df[,c(2,8,10)]
names(biomart_df) <- c("Chromosome_Biomarker", "Gene", "Chromosome_edit")
biomart_df = biomart_df[-which(duplicated(biomart_df$Gene)),]

# Join biomart_df file to the out_df file 
out_df <- merge(out_df, biomart_df, by = "Gene")

# # Check unique start sites here 
# unique(out_df$start) # 154
# plot(out_df$start)

names(out_df)[2] <- "Biomarker_TSS_Start"
names(out_df)[3] <- "Biomarker_TSS_End"
names(out_df)[1] <- "gene"

# Eyeball to make sure chr conversion has occured - fix any remaining formatting issues 
out_df$Chromosome_edit <- gsub("KIR", "", out_df$Chromosome_edit)
out_df$Chromosome_edit <- gsub("LRC", "", out_df$Chromosome_edit)
out_df$Chromosome_Biomarker <- NULL

# Replace the 4 trouble genes with the somamer name for joining as they are now extracted 
out_df$gene <- gsub("HDGFL3", "HDGFRP3", out_df$gene)
out_df$gene <- gsub("NT5C3B", "NT5C3L", out_df$gene)
out_df$gene <- gsub("ADA2", "CECR1", out_df$gene)
out_df$gene <- gsub("GOLM2", "CASC4", out_df$gene)
out_df$gene <- gsub("ADSS1", "ADSSL1", out_df$gene)

# [1] "HDGFRP3" "CECR1"   "NT5C3L"

# HDGFL3

# ADA2 

# NT5C3B


# row 159 - RBP5 is chr 12 - lookup online 
# row 60 - DEFB1 is chr 8 - lookup online 

out_df$Chromosome_edit[60] <- "8"
out_df$Chromosome_edit[159] <- "12"

# Merge back into main dataset with sig proteins 
tab2 <- left_join(tab, out_df, by = "gene") 

table(is.na(tab2))
# FALSE

write.csv(tab2, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/cis_trans/FULL_TSS_table_second_threshold.csv", row.names = F)


###############################

## WBC MODEL SECOND

cpgs <- WBC
unique <- unique(cpgs$SeqId)

length(unique)
# [1] 271

unique <- cpgs[c(1,10,11)] %>% unique() # xtracted into dataset with gene names 
tab <- unique
names(tab)[3] <- "gene" 

# get gene info set up from external source (Use Ensembl to Extract Relevant Features)
library(biomaRt)
ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
attributes = listAttributes(ensembl, page = "structure")
attributes[grep("transcript", attributes$description, ignore.case = TRUE), ]

# get list of proteins gene names 
list <- tab$gene

## Get Transcription Start Sites for soma proteins 
tss <- getBM(attributes = c("transcription_start_site", "chromosome_name",
                            "transcript_start", "transcript_end",
                            "strand",  "ensembl_gene_id",
                            "ensembl_transcript_id", "external_gene_name"),
             filters = "external_gene_name", "UCSC_RefGene_Group", values = c(list),
             mart = ensembl)

# see how many this worked for
length(unique(tss$external_gene_name)) 

# See if there are unique TSS extracted for these proteins 
dim(tss) 
length(unique(tss$transcription_start_site)) 

# which are the gene names that didnt extract
ov <- which(!list %in% tss$external_gene_name)
list2 <- list[ov]


# [1] "HDGFRP3"     "ALPPL2"      "CECR1"       "NT5C3L"      "ADSSL1"
# [6] "CASC4"       "HRSP12"      "PACAP"       "IL12B|IL23A"

# HDGFL3
# ALPG
# ADA2 
# NT5C3B
# ADSS1
# GOLM2
# RIDA
# ADCYAP1
# IL12B

# Try these aliases for the genes as per lookup on gene cards site 

list3 <- c("HDGFL3", "ALPG", "NT5C3B","ADA2", "ADSS1", "GOLM2", "RIDA", "ADCYAP1", "IL12B")

## Get Transcription Start Sites for remaining soma proteins 
tss2 <- getBM(attributes = c("transcription_start_site", "chromosome_name",
                            "transcript_start", "transcript_end",
                            "strand",  "ensembl_gene_id",
                            "ensembl_transcript_id", "external_gene_name"),
             filters = "external_gene_name", values = c(list3),
             mart = ensembl)

# see how many this worked for
length(unique(tss2$external_gene_name)) # 9

# Now we have a complete extracted set of tss information - join them together 

tss <- rbind(tss, tss2)

length(unique(tss$external_gene_name)) 

## set up code to create biomart dataframe with gene, start and end columns
biomart_df <- tss

names(biomart_df)[8] <- "gene"
names(biomart_df)[3] <- "start"
names(biomart_df)[4] <- "end"
names(biomart_df)[2] <- "chromosome_ensembl"

# Again, lets do a check to see if each protein has its own TSS extracted that is unique 

subset <- biomart_df[,c("transcription_start_site", "gene")]
subset <- unique(subset$gene)

## separate data into genes on positive strand and genes on negative strand

biomart_df_pos <- biomart_df[biomart_df$strand == 1,]
biomart_df_neg <- biomart_df[biomart_df$strand == -1,]

# sanity check to see if overlap is present/correct number of genes are present

a = unique(biomart_df_pos$gene)
b = unique(biomart_df_neg$gene)
which(a %in% b)
length(a)
length(b)

## create output data frames for genes on positive strand and those on negative strand to be r-bound later
   
   
    ##out_df1 == genes on positive strand
   
out_df1 <- matrix(nrow=length(unique(biomart_df_pos$gene)), ncol=2)
colnames(out_df1) <- c("start", "end")
rownames(out_df1) <- unique(biomart_df_pos$gene)


    ##out_df2 == genes on negative strand
   
out_df2 <- matrix(nrow=length(unique(biomart_df_neg$gene)), ncol=2)
colnames(out_df2) <- c("start", "end")
rownames(out_df2) <- unique(biomart_df_neg$gene)

## loop through gene names to fill in min 5'end (start) and max 3'end (end) as on positive strand you want furthest away 5'(tss) and biggest distance to 3' end
## fill out out_df1 based on the above

for(gene in unique(biomart_df_pos$gene)) {
tmp <- biomart_df_pos[which(biomart_df_pos$gene==gene), ]
out_df1[gene,"start"] <-  min(tmp$start)
out_df1[gene,"end"] <-  max(tmp$end)
}

## loop through gene names to fill in min 3'end (end) and max 5'end (start) as on negative strand you want furthest away 3'(tss (relative to positive strand)) and biggest distance to 5' end
## fill out out_df2 based on the above

for(gene in unique(biomart_df_neg$gene)) {
tmp <- biomart_df_neg[which(biomart_df_neg$gene==gene), ]
out_df2[gene,"start"] <-  max(tmp$end)
out_df2[gene,"end"] <-  min(tmp$start)
}

## row bind the positive and negative strand dataframes to give one out_df
out_df <- rbind(out_df1, out_df2)


## Tidy the files before saving for new row with protein names (out_df), then one for length of transcripts in biomart_df for future reference

    ## set as data frames
    out_df <- as.data.frame(out_df)
    biomart_df <- as.data.frame(biomart_df)

out_df$Gene <- row.names(out_df)
biomart_df$transcript_length <- abs(biomart_df$start - biomart_df$end)

# Check unique start sites here 
unique(out_df$start) # 143
# plot(out_df$start)

## Get Chromosome Number for Protein 
biomart_df$chromosome_ensembl2 <-  gsub("CHR_HSCHR", "", biomart_df$chromosome_ensembl)
biomart_df$chromosome_ensembl2 <- gsub("_.*", "", biomart_df$chromosome_ensembl2)
biomart_df <- biomart_df[,c(2,8,10)]
names(biomart_df) <- c("Chromosome_Biomarker", "Gene", "Chromosome_edit")
biomart_df = biomart_df[-which(duplicated(biomart_df$Gene)),]

# Join biomart_df file to the out_df file 
out_df <- merge(out_df, biomart_df, by = "Gene")

# # Check unique start sites here 
# unique(out_df$start) # 154
# plot(out_df$start)

names(out_df)[2] <- "Biomarker_TSS_Start"
names(out_df)[3] <- "Biomarker_TSS_End"
names(out_df)[1] <- "gene"

# Eyeball to make sure chr conversion has occured - fix any remaining formatting issues 
out_df$Chromosome_edit <- gsub("KIR", "", out_df$Chromosome_edit)
out_df$Chromosome_edit <- gsub("LRC", "", out_df$Chromosome_edit)
out_df$Chromosome_Biomarker <- NULL

# Replace the 4 trouble genes with the somamer name for joining as they are now extracted 
out_df$gene <- gsub("HDGFL3", "HDGFRP3", out_df$gene)
out_df$gene <- gsub("NT5C3B", "NT5C3L", out_df$gene)
out_df$gene <- gsub("ADA2", "CECR1", out_df$gene)
out_df$gene <- gsub("GOLM2", "CASC4", out_df$gene)
out_df$gene <- gsub("ALPG", "ALPPL2", out_df$gene)
out_df$gene <- gsub("ADSS1", "ADSSL1", out_df$gene)
out_df$gene <- gsub("RIDA", "HRSP12", out_df$gene)
out_df$gene <- gsub("ADCYAP1", "PACAP", out_df$gene)
out_df$gene <- gsub("IL12B", "IL12B|IL23A", out_df$gene)

# [1] "HDGFRP3"     "ALPPL2"      "CECR1"       "NT5C3L"      "ADSSL1"
# [6] "CASC4"       "HRSP12"      "PACAP"       "IL12B|IL23A"

# HDGFL3
# ALPG
# ADA2 
# NT5C3B
# ADSS1
# GOLM2
# RIDA
# ADCYAP1
# IL12B

# row 37 - C1Rl is chr 12 - lookup online 
# row 79 - DEFB1 is chr 8 - lookup online
# row 212 - RBP5 is chr 12 - lookup online
# row 226 - SEMA7A ischr 15 - lookup online 

out_df$Chromosome_edit[37] <- "12"
out_df$Chromosome_edit[79] <- "8"
out_df$Chromosome_edit[212] <- "12"
out_df$Chromosome_edit[226] <- "15"

# Merge back into main dataset with sig proteins 
tab2 <- left_join(tab, out_df, by = "gene") 

table(is.na(tab2))
# TRUE

# 184 assign to 246 to cover IL12B
tab2[184,c(4,5,6)] <- tab2[246,c(4,5,6)]

table(is.na(tab2))
# FALSE

write.csv(tab2, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/cis_trans/WBC_TSS_table_second_threshold.csv", row.names = F)


###############################

## BASIC MODEL THIRD

cpgs <- basic
unique <- unique(cpgs$SeqId)

length(unique)
# [1] 276

unique <- cpgs[c(1,10,11)] %>% unique() # xtracted into dataset with gene names 
tab <- unique
names(tab)[3] <- "gene" 

# get gene info set up from external source (Use Ensembl to Extract Relevant Features)
library(biomaRt)
ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
attributes = listAttributes(ensembl, page = "structure")
attributes[grep("transcript", attributes$description, ignore.case = TRUE), ]

# get list of proteins gene names 
list <- tab$gene

## Get Transcription Start Sites for soma proteins 
tss <- getBM(attributes = c("transcription_start_site", "chromosome_name",
                            "transcript_start", "transcript_end",
                            "strand",  "ensembl_gene_id",
                            "ensembl_transcript_id", "external_gene_name"),
             filters = "external_gene_name", "UCSC_RefGene_Group", values = c(list),
             mart = ensembl)

# see how many this worked for
length(unique(tss$external_gene_name)) # 188

# See if there are unique TSS extracted for these proteins 
dim(tss) 
length(unique(tss$transcription_start_site)) 

# which are the gene names that didnt extract
ov <- which(!list %in% tss$external_gene_name)
list2 <- list[ov]


 # [1] "ALPPL2"        "HDGFRP3"       "CECR1"         "LTA|LTB"
 # [5] "NT5C3L"        "ADSSL1"        "FAIM3"         "PACAP"
 # [9] "S100A8|S100A9" "CASC4"         "HRSP12"


# ALPG
# HDGFL3
# ADA2 
# LTA
# NT5C3B
# ADSS1
# FCMR
# ADCYAP1
# S100A8
# GOLM2
# RIDA

# Try these aliases for the genes as per lookup on gene cards site 

list3 <- c("ALPG","HDGFL3","ADA2", "LTA", "NT5C3B","ADSS1","FCMR", "ADCYAP1", "S100A8", "GOLM2", "RIDA")


## Get Transcription Start Sites for remaining soma proteins 
tss2 <- getBM(attributes = c("transcription_start_site", "chromosome_name",
                            "transcript_start", "transcript_end",
                            "strand",  "ensembl_gene_id",
                            "ensembl_transcript_id", "external_gene_name"),
             filters = "external_gene_name", values = c(list3),
             mart = ensembl)

# see how many this worked for
length(unique(tss2$external_gene_name)) # 11

# Now we have a complete extracted set of tss information - join them together 

tss <- rbind(tss, tss2)

length(unique(tss$external_gene_name)) 

## set up code to create biomart dataframe with gene, start and end columns
biomart_df <- tss

names(biomart_df)[8] <- "gene"
names(biomart_df)[3] <- "start"
names(biomart_df)[4] <- "end"
names(biomart_df)[2] <- "chromosome_ensembl"

# Again, lets do a check to see if each protein has its own TSS extracted that is unique 

subset <- biomart_df[,c("transcription_start_site", "gene")]
subset <- unique(subset$gene)

## separate data into genes on positive strand and genes on negative strand

biomart_df_pos <- biomart_df[biomart_df$strand == 1,]
biomart_df_neg <- biomart_df[biomart_df$strand == -1,]

# sanity check to see if overlap is present/correct number of genes are present

a = unique(biomart_df_pos$gene)
b = unique(biomart_df_neg$gene)
which(a %in% b)
length(a)
length(b)

## create output data frames for genes on positive strand and those on negative strand to be r-bound later
   
   
    ##out_df1 == genes on positive strand
   
out_df1 <- matrix(nrow=length(unique(biomart_df_pos$gene)), ncol=2)
colnames(out_df1) <- c("start", "end")
rownames(out_df1) <- unique(biomart_df_pos$gene)


    ##out_df2 == genes on negative strand
   
out_df2 <- matrix(nrow=length(unique(biomart_df_neg$gene)), ncol=2)
colnames(out_df2) <- c("start", "end")
rownames(out_df2) <- unique(biomart_df_neg$gene)

## loop through gene names to fill in min 5'end (start) and max 3'end (end) as on positive strand you want furthest away 5'(tss) and biggest distance to 3' end
## fill out out_df1 based on the above

for(gene in unique(biomart_df_pos$gene)) {
tmp <- biomart_df_pos[which(biomart_df_pos$gene==gene), ]
out_df1[gene,"start"] <-  min(tmp$start)
out_df1[gene,"end"] <-  max(tmp$end)
}

## loop through gene names to fill in min 3'end (end) and max 5'end (start) as on negative strand you want furthest away 3'(tss (relative to positive strand)) and biggest distance to 5' end
## fill out out_df2 based on the above

for(gene in unique(biomart_df_neg$gene)) {
tmp <- biomart_df_neg[which(biomart_df_neg$gene==gene), ]
out_df2[gene,"start"] <-  max(tmp$end)
out_df2[gene,"end"] <-  min(tmp$start)
}

## row bind the positive and negative strand dataframes to give one out_df
out_df <- rbind(out_df1, out_df2)


## Tidy the files before saving for new row with protein names (out_df), then one for length of transcripts in biomart_df for future reference

    ## set as data frames
    out_df <- as.data.frame(out_df)
    biomart_df <- as.data.frame(biomart_df)

out_df$Gene <- row.names(out_df)
biomart_df$transcript_length <- abs(biomart_df$start - biomart_df$end)

# Check unique start sites here 
unique(out_df$start) # 143
# plot(out_df$start)

## Get Chromosome Number for Protein 
biomart_df$chromosome_ensembl2 <-  gsub("CHR_HSCHR", "", biomart_df$chromosome_ensembl)
biomart_df$chromosome_ensembl2 <- gsub("_.*", "", biomart_df$chromosome_ensembl2)
biomart_df <- biomart_df[,c(2,8,10)]
names(biomart_df) <- c("Chromosome_Biomarker", "Gene", "Chromosome_edit")
biomart_df = biomart_df[-which(duplicated(biomart_df$Gene)),]

# Join biomart_df file to the out_df file 
out_df <- merge(out_df, biomart_df, by = "Gene")

# # Check unique start sites here 
# unique(out_df$start) # 154
# plot(out_df$start)

names(out_df)[2] <- "Biomarker_TSS_Start"
names(out_df)[3] <- "Biomarker_TSS_End"
names(out_df)[1] <- "gene"

# Eyeball to make sure chr conversion has occured - fix any remaining formatting issues 
out_df$Chromosome_edit <- gsub("KIR", "", out_df$Chromosome_edit)
out_df$Chromosome_edit <- gsub("LRC", "", out_df$Chromosome_edit)
out_df$Chromosome_Biomarker <- NULL

# Replace the 4 trouble genes with the somamer name for joining as they are now extracted 
out_df$gene <- gsub("HDGFL3", "HDGFRP3", out_df$gene)
out_df$gene <- gsub("NT5C3B", "NT5C3L", out_df$gene)
out_df$gene <- gsub("ADA2", "CECR1", out_df$gene)
out_df$gene <- gsub("GOLM2", "CASC4", out_df$gene)
out_df$gene <- gsub("ALPG", "ALPPL2", out_df$gene)
out_df$gene <- gsub("ADSS1", "ADSSL1", out_df$gene)
out_df$gene <- gsub("FCMR", "FAIM3", out_df$gene)
out_df$gene <- gsub("ADCYAP1", "PACAP", out_df$gene)
out_df$gene <- gsub("S100A8", "S100A8|S100A9", out_df$gene)
out_df$gene <- gsub("LTA", "LTA|LTB", out_df$gene)
out_df$gene <- gsub("RIDA", "HRSP12", out_df$gene)

# row 37 - CIRL is chr 12 - lookup online
# row 80 - DEFB1 is chr 8 - lookup online
# row 217 - RBP5 is chr 12 - lookup online 
# row 229 - SAA1 is chr 11 - lookup online 

out_df$Chromosome_edit[37] <- "12"
out_df$Chromosome_edit[80] <- "8"
out_df$Chromosome_edit[217] <- "12"
out_df$Chromosome_edit[229] <- "11"

# Merge back into main dataset with sig proteins 
tab2 <- left_join(tab, out_df, by = "gene") 

table(is.na(tab2))
# FALSE

write.csv(tab2, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/cis_trans/BASIC_TSS_table_second_threshold.csv", row.names = F)


#############################################################################################

### Now look at mapping CIS and TRANS associations for the cpgs that associated with the proteins 

# I will do this for basic, WBC and FULL models 

## FULL

# Read in the cpg positional info 
cpgs <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Results_collation/FULL_EWAS_results_second_threshold_020221.csv")

# Read in the protein gene annotations as generated above 
prot <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/cis_trans/FULL_TSS_table_second_threshold.csv")

# Merge gene positional info into the table with weights based on the gene names 
table <- left_join(cpgs, prot, by = c("SeqId"))
dim(table) # 2080 - this looks good 

# See whether we still have unique TSS sites 
length(unique(table$Biomarker_TSS_Start)) # 191 protein start sites 
length(unique(table$bp)) # 1837 cpg sites

# Align column names to cis trans code 

chrom <- table
names(chrom)[15] <- "gene_start"
names(chrom)[16] <- "gene_end"
names(chrom)[4] <- "Cpg_pos"

chrom$gene_start <- as.numeric(chrom$gene_start)
chrom$gene_end <- as.numeric(chrom$gene_end)

# Look at which matches are above and below in the range 
chrom$diff <- chrom$gene_start - chrom$Cpg_pos 
chrom$diff2 <- chrom$gene_start + chrom$Cpg_pos

# Change negative differences to absolute values 
chrom$diff <- abs(chrom$diff)
chrom$diff2 <- abs(chrom$diff2)

# match chromosome naming columns
names(chrom)[2] <- "chr_cpg"
names(chrom)[17] <- "Chr"

# Do everything in one go to assign CIS or TRANS 
chrom <- chrom %>% mutate(Effect = case_when(
  chr_cpg == Chr & (chrom$diff < 10000000 | chrom$diff2 < 10000000) ~ "CIS",
  chr_cpg == Chr & (!chrom$diff < 10000000 | chrom$diff < 10000000) ~ "TRANS",
  !chr_cpg == Chr ~ "TRANS"
))


# Make sure you get same results as decomposed version below (yes, it matches up)
test <- chrom %>% filter(Effect == "CIS") 
dim(test) # 451
test <- chrom %>% filter(Effect == "TRANS") 
dim(test) # 2477


# Write out table for use in next script with basic naming 
write.csv(chrom, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/cis_trans/FULL_cistrans_table_unformatted_naming_second_threshold.csv", row.names = F)

# Write the table for the suppl file
chrom <- chrom[c(3,2,5,4,6:9,1,11,10,12,15,16,17,20)]
names(chrom) <- c("CpG", "Chromsome of CpG", "Gene of CpG", "CpG Position", "Orientation",
  "Beta", "SE", "P", "SeqId", "Gene of Protein", "UniProt", "UniProt Full Name", "Gene Start", "Gene End", "Chromosome of Gene", "Association Type")
chrom <- chrom[order(chrom$P),]
write.csv(chrom, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/cis_trans/FULL_cistrans_table_formatted_naming_second_threshold.csv", row.names = F)

########################

## WBC

# Read in the cpg positional info 
cpgs <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Results_collation/WBC_EWAS_results_second_threshold_020221.csv")

# Read in the protein gene annotations as generated above 
prot <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/cis_trans/WBC_TSS_table_second_threshold.csv")

# Merge gene positional info into the table with weights based on the gene names 
table <- left_join(cpgs, prot, by = c("SeqId"))
dim(table) # 3213 - this looks good 

# See whether we still have unique TSS sites 
length(unique(table$Biomarker_TSS_Start)) # 262 protein start sites 
length(unique(table$bp)) # 1958 cpg sites

# Align column names to cis trans code 

chrom <- table
names(chrom)[15] <- "gene_start"
names(chrom)[16] <- "gene_end"
names(chrom)[4] <- "Cpg_pos"

chrom$gene_start <- as.numeric(chrom$gene_start)
chrom$gene_end <- as.numeric(chrom$gene_end)

# Look at which matches are above and below in the range 
chrom$diff <- chrom$gene_start - chrom$Cpg_pos 
chrom$diff2 <- chrom$gene_start + chrom$Cpg_pos

# Change negative differences to absolute values 
chrom$diff <- abs(chrom$diff)
chrom$diff2 <- abs(chrom$diff2)

# match chromosome naming columns
names(chrom)[2] <- "chr_cpg"
names(chrom)[17] <- "Chr"

# Do everything in one go to assign CIS or TRANS 
chrom <- chrom %>% mutate(Effect = case_when(
  chr_cpg == Chr & (chrom$diff < 10000000 | chrom$diff2 < 10000000) ~ "CIS",
  chr_cpg == Chr & (!chrom$diff < 10000000 | chrom$diff < 10000000) ~ "TRANS",
  !chr_cpg == Chr ~ "TRANS"
))


# Make sure you get same results as decomposed version below (yes, it matches up)
test <- chrom %>% filter(Effect == "CIS") 
dim(test) # 453
test <- chrom %>% filter(Effect == "TRANS") 
dim(test) # 2760


# Write out table for use in next script with basic naming 
write.csv(chrom, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/cis_trans/WBC_cistrans_table_unformatted_naming_second_threshold.csv", row.names = F)

# Write the table for the suppl file
chrom <- chrom[c(3,2,5,4,6:9,1,11,10,12,15,16,17,20)]
names(chrom) <- c("CpG", "Chromsome of CpG", "Gene of CpG", "CpG Position", "Orientation",
  "Beta", "SE", "P", "SeqId", "Gene of Protein", "UniProt", "UniProt Full Name", "Gene Start", "Gene End", "Chromosome of Gene", "Association Type")
chrom <- chrom[order(chrom$P),]
write.csv(chrom, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/cis_trans/WBC_cistrans_table_formatted_naming_second_threshold.csv", row.names = F)

########################

## BASIC

# Read in the cpg positional info 
cpgs <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Results_collation/Basic_EWAS_results_second_threshold_020221.csv")

# Read in the protein gene annotations as generated above 
prot <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/cis_trans/BASIC_TSS_table_second_threshold.csv")

# Merge gene positional info into the table with weights based on the gene names 
table <- left_join(cpgs, prot, by = c("SeqId"))
dim(table) # 238245 - this looks good 

# See whether we still have unique TSS sites 
length(unique(table$Biomarker_TSS_Start)) # 269 protein start sites 
length(unique(table$bp)) # 94767 cpg sites

# Align column names to cis trans code 

chrom <- table
names(chrom)[15] <- "gene_start"
names(chrom)[16] <- "gene_end"
names(chrom)[4] <- "Cpg_pos"

chrom$gene_start <- as.numeric(chrom$gene_start)
chrom$gene_end <- as.numeric(chrom$gene_end)

# Look at which matches are above and below in the range 
chrom$diff <- chrom$gene_start - chrom$Cpg_pos 
chrom$diff2 <- chrom$gene_start + chrom$Cpg_pos

# Change negative differences to absolute values 
chrom$diff <- abs(chrom$diff)
chrom$diff2 <- abs(chrom$diff2)

# match chromosome naming columns
names(chrom)[2] <- "chr_cpg"
names(chrom)[17] <- "Chr"

# Do everything in one go to assign CIS or TRANS 
chrom <- chrom %>% mutate(Effect = case_when(
  chr_cpg == Chr & (chrom$diff < 10000000 | chrom$diff2 < 10000000) ~ "CIS",
  chr_cpg == Chr & (!chrom$diff < 10000000 | chrom$diff < 10000000) ~ "TRANS",
  !chr_cpg == Chr ~ "TRANS"
))


# Make sure you get same results as decomposed version below (yes, it matches up)
test <- chrom %>% filter(Effect == "CIS") 
dim(test) # 2107
test <- chrom %>% filter(Effect == "TRANS") 
dim(test) # 236138


# Write out table for use in next script with basic naming 
write.csv(chrom, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/cis_trans/BASIC_cistrans_table_unformatted_naming_second_threshold.csv", row.names = F)

# Write the table for the suppl file
chrom <- chrom[c(3,2,5,4,6:9,1,11,10,12,15,16,17,20)]
names(chrom) <- c("CpG", "Chromsome of CpG", "Gene of CpG", "CpG Position", "Orientation",
  "Beta", "SE", "P", "SeqId", "Gene of Protein", "UniProt", "UniProt Full Name", "Gene Start", "Gene End", "Chromosome of Gene", "Association Type")
chrom <- chrom[order(chrom$P),]
write.csv(chrom, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/cis_trans/BASIC_cistrans_table_formatted_naming_second_threshold.csv", row.names = F)



###############################################################################

### UPDATE: REVISIONS

# Load in all of the above script for TSS cis/trans calculation, but set the threshold to 1mb rather than 10mb as a sensitivity 
# We will check how many pQTMs fall outside of this threshold

## FULL

# Read in the cpg positional info 
cpgs <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Results_collation/FULL_EWAS_results_020221.csv")

# Read in the protein gene annotations as generated above 
prot <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/cis_trans/FULL_TSS_table.csv")

# Merge gene positional info into the table with weights based on the gene names 
table <- left_join(cpgs, prot, by = c("SeqId"))
dim(table) # 2080 - this looks good 

# See whether we still have unique TSS sites 
length(unique(table$Biomarker_TSS_Start)) # 143 protein start sites 
length(unique(table$bp)) # 1305 cpg sites

# Align column names to cis trans code 

chrom <- table
names(chrom)[15] <- "gene_start"
names(chrom)[16] <- "gene_end"
names(chrom)[4] <- "Cpg_pos"

chrom$gene_start <- as.numeric(chrom$gene_start)
chrom$gene_end <- as.numeric(chrom$gene_end)

# Look at which matches are above and below in the range 
chrom$diff <- chrom$gene_start - chrom$Cpg_pos 
chrom$diff2 <- chrom$gene_start + chrom$Cpg_pos

# Change negative differences to absolute values 
chrom$diff <- abs(chrom$diff)
chrom$diff2 <- abs(chrom$diff2)

# match chromosome naming columns
names(chrom)[2] <- "chr_cpg"
names(chrom)[17] <- "Chr"

# Do everything in one go to assign CIS or TRANS 
chrom <- chrom %>% mutate(Effect = case_when(
  chr_cpg == Chr & (chrom$diff < 1000000 | chrom$diff2 < 1000000) ~ "CIS",
  chr_cpg == Chr & (!chrom$diff < 1000000 | chrom$diff < 1000000) ~ "TRANS",
  !chr_cpg == Chr ~ "TRANS"
))


# Make sure you get same results as decomposed version below (yes, it matches up)
test <- chrom %>% filter(Effect == "CIS") 
dim(test) # 350 - down to 318 
test <- chrom %>% filter(Effect == "TRANS") 
dim(test) # 1730 - down to 1762


# Write out table for use in next script with basic naming 
write.csv(chrom, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/cis_trans/1mb_FULL_cistrans_table_unformatted_naming.csv", row.names = F)

# Write the table for the suppl file
chrom <- chrom[c(3,2,5,4,6:9,1,11,10,12,15,16,17,20)]
names(chrom) <- c("CpG", "Chromsome of CpG", "Gene of CpG", "CpG Position", "Orientation",
  "Beta", "SE", "P", "SeqId", "Gene of Protein", "UniProt", "UniProt Full Name", "Gene Start", "Gene End", "Chromosome of Gene", "Association Type")
chrom <- chrom[order(chrom$P),]
write.csv(chrom, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/cis_trans/1mb_FULL_cistrans_table_formatted_naming.csv", row.names = F)


# Take a look at which cis associations (10mb) are altered to trans (in the 1mb sensitivity threshold) in full model
# 41 associations change overall
mb1 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/cis_trans/1mb_FULL_cistrans_table_unformatted_naming.csv")
mb10 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/cis_trans/FULL_cistrans_table_unformatted_naming.csv")

cis_mb1 <- mb1 %>% filter(Effect == "CIS")
cis_mb1$assoc <- paste0(cis_mb1$Probe, "_", cis_mb1$gene)

cis_mb10 <- mb10 %>% filter(Effect == "CIS")
cis_mb10$assoc <- paste0(cis_mb10$Probe, "_", cis_mb10$gene)


diff <- cis_mb10[-which(cis_mb10$assoc %in% cis_mb1$assoc),]

diff$assoc
#  [1] "cg04601780_ACE"    "cg20491697_COL6A1" "cg13155430_MX1"
#  [4] "cg14004602_CCL15"  "cg21549285_MX1"    "cg05300717_PRG3"
#  [7] "cg26312951_MX1"    "cg22862003_MX1"    "cg12482297_ACE"
# [10] "cg25087851_PRG3"   "cg10260205_COL6A1" "cg02487088_PRG3"
# [13] "cg21693883_COL6A1" "cg18367488_PRG3"   "cg07318398_NT5C3L"
# [16] "cg26657643_PRG3"   "cg11205545_PRG3"   "cg07636225_PRG3"
# [19] "cg16785077_MX1"    "cg24515533_NT5C3L" "cg14551643_AIF1L"
# [22] "cg03448301_PRG3"   "cg11884243_FCN2"   "cg18197277_PRG3"
# [25] "cg21203249_SEMA4D" "cg01556593_PRG3"   "cg14127188_FCN2"
# [28] "cg25116412_NT5C3L" "cg03776060_AIF1L"  "cg19502237_FCN2"
# [31] "cg21120609_COL6A1" "cg14631576_OGN"



########################

## WBC

# Read in the cpg positional info 
cpgs <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Results_collation/WBC_EWAS_results_020221.csv")

# Read in the protein gene annotations as generated above 
prot <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/cis_trans/WBC_TSS_table.csv")

# Merge gene positional info into the table with weights based on the gene names 
table <- left_join(cpgs, prot, by = c("SeqId"))
dim(table) # 2230 - this looks good 

# See whether we still have unique TSS sites 
length(unique(table$Biomarker_TSS_Start)) # 185 protein start sites 
length(unique(table$bp)) # 1362 cpg sites

# Align column names to cis trans code 

chrom <- table
names(chrom)[15] <- "gene_start"
names(chrom)[16] <- "gene_end"
names(chrom)[4] <- "Cpg_pos"

chrom$gene_start <- as.numeric(chrom$gene_start)
chrom$gene_end <- as.numeric(chrom$gene_end)

# Look at which matches are above and below in the range 
chrom$diff <- chrom$gene_start - chrom$Cpg_pos 
chrom$diff2 <- chrom$gene_start + chrom$Cpg_pos

# Change negative differences to absolute values 
chrom$diff <- abs(chrom$diff)
chrom$diff2 <- abs(chrom$diff2)

# match chromosome naming columns
names(chrom)[2] <- "chr_cpg"
names(chrom)[17] <- "Chr"

# Do everything in one go to assign CIS or TRANS 
chrom <- chrom %>% mutate(Effect = case_when(
  chr_cpg == Chr & (chrom$diff < 1000000 | chrom$diff2 < 1000000) ~ "CIS",
  chr_cpg == Chr & (!chrom$diff < 1000000 | chrom$diff < 1000000) ~ "TRANS",
  !chr_cpg == Chr ~ "TRANS"
))


# Make sure you get same results as decomposed version below (yes, it matches up)
test <- chrom %>% filter(Effect == "CIS") 
dim(test) # 360 - down to 327
test <- chrom %>% filter(Effect == "TRANS") 
dim(test) # 1870 - down to 1903


# Write out table for use in next script with basic naming 
write.csv(chrom, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/cis_trans/1mb_WBC_cistrans_table_unformatted_naming.csv", row.names = F)

# Write the table for the suppl file
chrom <- chrom[c(3,2,5,4,6:9,1,11,10,12,15,16,17,20)]
names(chrom) <- c("CpG", "Chromsome of CpG", "Gene of CpG", "CpG Position", "Orientation",
  "Beta", "SE", "P", "SeqId", "Gene of Protein", "UniProt", "UniProt Full Name", "Gene Start", "Gene End", "Chromosome of Gene", "Association Type")
chrom <- chrom[order(chrom$P),]
write.csv(chrom, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/cis_trans/1mb_WBC_cistrans_table_formatted_naming.csv", row.names = F)

########################

## BASIC

# Read in the cpg positional info 
cpgs <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Results_collation/Basic_EWAS_results_020221.csv")

# Read in the protein gene annotations as generated above 
prot <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/cis_trans/BASIC_TSS_table.csv")

# Merge gene positional info into the table with weights based on the gene names 
table <- left_join(cpgs, prot, by = c("SeqId"))
dim(table) # 158858 - this looks good 

# See whether we still have unique TSS sites 
length(unique(table$Biomarker_TSS_Start)) # 198 protein start sites 
length(unique(table$bp)) # 72874 cpg sites

# Align column names to cis trans code 

chrom <- table
names(chrom)[15] <- "gene_start"
names(chrom)[16] <- "gene_end"
names(chrom)[4] <- "Cpg_pos"

chrom$gene_start <- as.numeric(chrom$gene_start)
chrom$gene_end <- as.numeric(chrom$gene_end)

# Look at which matches are above and below in the range 
chrom$diff <- chrom$gene_start - chrom$Cpg_pos 
chrom$diff2 <- chrom$gene_start + chrom$Cpg_pos

# Change negative differences to absolute values 
chrom$diff <- abs(chrom$diff)
chrom$diff2 <- abs(chrom$diff2)

# match chromosome naming columns
names(chrom)[2] <- "chr_cpg"
names(chrom)[17] <- "Chr"

# Do everything in one go to assign CIS or TRANS 
chrom <- chrom %>% mutate(Effect = case_when(
  chr_cpg == Chr & (chrom$diff < 1000000 | chrom$diff2 < 1000000) ~ "CIS",
  chr_cpg == Chr & (!chrom$diff < 1000000 | chrom$diff < 1000000) ~ "TRANS",
  !chr_cpg == Chr ~ "TRANS"
))


# Make sure you get same results as decomposed version below (yes, it matches up)
test <- chrom %>% filter(Effect == "CIS") 
dim(test) # 1399 - down to 555
test <- chrom %>% filter(Effect == "TRANS") 
dim(test) # 157459 - up to 158303


# Write out table for use in next script with basic naming 
write.csv(chrom, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/cis_trans/1mb_BASIC_cistrans_table_unformatted_naming.csv", row.names = F)

# Write the table for the suppl file
chrom <- chrom[c(3,2,5,4,6:9,1,11,10,12,15,16,17,20)]
names(chrom) <- c("CpG", "Chromsome of CpG", "Gene of CpG", "CpG Position", "Orientation",
  "Beta", "SE", "P", "SeqId", "Gene of Protein", "UniProt", "UniProt Full Name", "Gene Start", "Gene End", "Chromosome of Gene", "Association Type")
chrom <- chrom[order(chrom$P),]
write.csv(chrom, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/cis_trans/1mb_BASIC_cistrans_table_formatted_naming.csv", row.names = F)



###############################################################################

### UPDATE: REVISIONS - BUT DO FOR SECOND THRESHOLD

# Load in all of the above script for TSS cis/trans calculation, but set the threshold to 1mb rather than 10mb as a sensitivity 
# We will check how many pQTMs fall outside of this threshold

## FULL

# Read in the cpg positional info 
cpgs <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Results_collation/FULL_EWAS_results_second_threshold_020221.csv")

# Read in the protein gene annotations as generated above 
prot <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/cis_trans/FULL_TSS_table_second_threshold.csv")

# Merge gene positional info into the table with weights based on the gene names 
table <- left_join(cpgs, prot, by = c("SeqId"))
dim(table) 

# See whether we still have unique TSS sites 
length(unique(table$Biomarker_TSS_Start)) 
length(unique(table$bp)) 

# Align column names to cis trans code 

chrom <- table
names(chrom)[15] <- "gene_start"
names(chrom)[16] <- "gene_end"
names(chrom)[4] <- "Cpg_pos"

chrom$gene_start <- as.numeric(chrom$gene_start)
chrom$gene_end <- as.numeric(chrom$gene_end)

# Look at which matches are above and below in the range 
chrom$diff <- chrom$gene_start - chrom$Cpg_pos 
chrom$diff2 <- chrom$gene_start + chrom$Cpg_pos

# Change negative differences to absolute values 
chrom$diff <- abs(chrom$diff)
chrom$diff2 <- abs(chrom$diff2)

# match chromosome naming columns
names(chrom)[2] <- "chr_cpg"
names(chrom)[17] <- "Chr"

# Do everything in one go to assign CIS or TRANS 
chrom <- chrom %>% mutate(Effect = case_when(
  chr_cpg == Chr & (chrom$diff < 1000000 | chrom$diff2 < 1000000) ~ "CIS",
  chr_cpg == Chr & (!chrom$diff < 1000000 | chrom$diff < 1000000) ~ "TRANS",
  !chr_cpg == Chr ~ "TRANS"
))


# Make sure you get same results as decomposed version below (yes, it matches up)
test <- chrom %>% filter(Effect == "CIS") 
dim(test) # 409
test <- chrom %>% filter(Effect == "TRANS") 
dim(test) # 2519


# # Write out table for use in next script with basic naming 
# write.csv(chrom, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/cis_trans/1mb_FULL_cistrans_table_unformatted_naming.csv", row.names = F)

# # Write the table for the suppl file
# chrom <- chrom[c(3,2,5,4,6:9,1,11,10,12,15,16,17,20)]
# names(chrom) <- c("CpG", "Chromsome of CpG", "Gene of CpG", "CpG Position", "Orientation",
#   "Beta", "SE", "P", "SeqId", "Gene of Protein", "UniProt", "UniProt Full Name", "Gene Start", "Gene End", "Chromosome of Gene", "Association Type")
# chrom <- chrom[order(chrom$P),]
# write.csv(chrom, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/cis_trans/1mb_FULL_cistrans_table_formatted_naming.csv", row.names = F)


# # Take a look at which cis associations (10mb) are altered to trans (in the 1mb sensitivity threshold) in full model
# # 41 associations change overall
# mb1 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/cis_trans/1mb_FULL_cistrans_table_unformatted_naming.csv")
# mb10 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/cis_trans/FULL_cistrans_table_unformatted_naming.csv")

# cis_mb1 <- mb1 %>% filter(Effect == "CIS")
# cis_mb1$assoc <- paste0(cis_mb1$Probe, "_", cis_mb1$gene)

# cis_mb10 <- mb10 %>% filter(Effect == "CIS")
# cis_mb10$assoc <- paste0(cis_mb10$Probe, "_", cis_mb10$gene)


# diff <- cis_mb10[-which(cis_mb10$assoc %in% cis_mb1$assoc),]

# diff$assoc
# #  [1] "cg04601780_ACE"    "cg20491697_COL6A1" "cg13155430_MX1"
# #  [4] "cg14004602_CCL15"  "cg21549285_MX1"    "cg05300717_PRG3"
# #  [7] "cg26312951_MX1"    "cg22862003_MX1"    "cg12482297_ACE"
# # [10] "cg25087851_PRG3"   "cg10260205_COL6A1" "cg02487088_PRG3"
# # [13] "cg21693883_COL6A1" "cg18367488_PRG3"   "cg07318398_NT5C3L"
# # [16] "cg26657643_PRG3"   "cg11205545_PRG3"   "cg07636225_PRG3"
# # [19] "cg16785077_MX1"    "cg24515533_NT5C3L" "cg14551643_AIF1L"
# # [22] "cg03448301_PRG3"   "cg11884243_FCN2"   "cg18197277_PRG3"
# # [25] "cg21203249_SEMA4D" "cg01556593_PRG3"   "cg14127188_FCN2"
# # [28] "cg25116412_NT5C3L" "cg03776060_AIF1L"  "cg19502237_FCN2"
# # [31] "cg21120609_COL6A1" "cg14631576_OGN"



########################

## WBC

# Read in the cpg positional info 
cpgs <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Results_collation/WBC_EWAS_results_second_threshold_020221.csv")

# Read in the protein gene annotations as generated above 
prot <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/cis_trans/WBC_TSS_table_second_threshold.csv")

# Merge gene positional info into the table with weights based on the gene names 
table <- left_join(cpgs, prot, by = c("SeqId"))
dim(table) 

# See whether we still have unique TSS sites 
length(unique(table$Biomarker_TSS_Start)) 
length(unique(table$bp)) 

# Align column names to cis trans code 

chrom <- table
names(chrom)[15] <- "gene_start"
names(chrom)[16] <- "gene_end"
names(chrom)[4] <- "Cpg_pos"

chrom$gene_start <- as.numeric(chrom$gene_start)
chrom$gene_end <- as.numeric(chrom$gene_end)

# Look at which matches are above and below in the range 
chrom$diff <- chrom$gene_start - chrom$Cpg_pos 
chrom$diff2 <- chrom$gene_start + chrom$Cpg_pos

# Change negative differences to absolute values 
chrom$diff <- abs(chrom$diff)
chrom$diff2 <- abs(chrom$diff2)

# match chromosome naming columns
names(chrom)[2] <- "chr_cpg"
names(chrom)[17] <- "Chr"

# Do everything in one go to assign CIS or TRANS 
chrom <- chrom %>% mutate(Effect = case_when(
  chr_cpg == Chr & (chrom$diff < 1000000 | chrom$diff2 < 1000000) ~ "CIS",
  chr_cpg == Chr & (!chrom$diff < 1000000 | chrom$diff < 1000000) ~ "TRANS",
  !chr_cpg == Chr ~ "TRANS"
))


# Make sure you get same results as decomposed version below (yes, it matches up)
test <- chrom %>% filter(Effect == "CIS") 
dim(test) # 413
test <- chrom %>% filter(Effect == "TRANS") 
dim(test) # 2800


# # Write out table for use in next script with basic naming 
# write.csv(chrom, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/cis_trans/1mb_WBC_cistrans_table_unformatted_naming.csv", row.names = F)

# # Write the table for the suppl file
# chrom <- chrom[c(3,2,5,4,6:9,1,11,10,12,15,16,17,20)]
# names(chrom) <- c("CpG", "Chromsome of CpG", "Gene of CpG", "CpG Position", "Orientation",
#   "Beta", "SE", "P", "SeqId", "Gene of Protein", "UniProt", "UniProt Full Name", "Gene Start", "Gene End", "Chromosome of Gene", "Association Type")
# chrom <- chrom[order(chrom$P),]
# write.csv(chrom, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/cis_trans/1mb_WBC_cistrans_table_formatted_naming.csv", row.names = F)

########################

## BASIC

# Read in the cpg positional info 
cpgs <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/00_Results_collation/Basic_EWAS_results_second_threshold_020221.csv")

# Read in the protein gene annotations as generated above 
prot <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/cis_trans/BASIC_TSS_table_second_threshold.csv")

# Merge gene positional info into the table with weights based on the gene names 
table <- left_join(cpgs, prot, by = c("SeqId"))
dim(table) # 158858 - this looks good 

# See whether we still have unique TSS sites 
length(unique(table$Biomarker_TSS_Start)) # 198 protein start sites 
length(unique(table$bp)) # 72874 cpg sites

# Align column names to cis trans code 

chrom <- table
names(chrom)[15] <- "gene_start"
names(chrom)[16] <- "gene_end"
names(chrom)[4] <- "Cpg_pos"

chrom$gene_start <- as.numeric(chrom$gene_start)
chrom$gene_end <- as.numeric(chrom$gene_end)

# Look at which matches are above and below in the range 
chrom$diff <- chrom$gene_start - chrom$Cpg_pos 
chrom$diff2 <- chrom$gene_start + chrom$Cpg_pos

# Change negative differences to absolute values 
chrom$diff <- abs(chrom$diff)
chrom$diff2 <- abs(chrom$diff2)

# match chromosome naming columns
names(chrom)[2] <- "chr_cpg"
names(chrom)[17] <- "Chr"

# Do everything in one go to assign CIS or TRANS 
chrom <- chrom %>% mutate(Effect = case_when(
  chr_cpg == Chr & (chrom$diff < 1000000 | chrom$diff2 < 1000000) ~ "CIS",
  chr_cpg == Chr & (!chrom$diff < 1000000 | chrom$diff < 1000000) ~ "TRANS",
  !chr_cpg == Chr ~ "TRANS"
))


# Make sure you get same results as decomposed version below (yes, it matches up)
test <- chrom %>% filter(Effect == "CIS") 
dim(test) # 752
test <- chrom %>% filter(Effect == "TRANS") 
dim(test) # 237493


# # Write out table for use in next script with basic naming 
# write.csv(chrom, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/cis_trans/1mb_BASIC_cistrans_table_unformatted_naming.csv", row.names = F)

# # Write the table for the suppl file
# chrom <- chrom[c(3,2,5,4,6:9,1,11,10,12,15,16,17,20)]
# names(chrom) <- c("CpG", "Chromsome of CpG", "Gene of CpG", "CpG Position", "Orientation",
#   "Beta", "SE", "P", "SeqId", "Gene of Protein", "UniProt", "UniProt Full Name", "Gene Start", "Gene End", "Chromosome of Gene", "Association Type")
# chrom <- chrom[order(chrom$P),]
# write.csv(chrom, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/cis_trans/1mb_BASIC_cistrans_table_formatted_naming.csv", row.names = F)


###############################################################################

### UPDATE: REVISIONS 250121

# Annotations for cpgs implicated in associations to be added using the UCSC reference database for EPIC array 
# Lookup promoter/genebody/TFsites - UCSC_RefGene_Group and Relation_to_CGI_Island

### FULL MODEL ONLY 

screen

R

# Read in the UCSC cpg annotations file 
anno <- readRDS("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Annotations_UCSC/annotations_UCSC_object_df_from_daniel.rds")

# Read in the annotated supplementary table with all assocations in fully adjusted EWAS
d <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/cis_trans/FULL_cistrans_table_formatted_naming_second_threshold.csv")
dim(d) # 2928

# Subset annotations file to cpgs involved in associations
anno <-  anno[which(rownames(anno) %in% d$CpG),]
dim(anno) # 1837 cpgs 

# Select columns of interest from the annotation file 
anno <- anno[,c("chr", "pos", "strand", "UCSC_RefGene_Group", "Islands_Name", "UCSC_RefGene_Name", 
    "Relation_to_Island", "OpenChromatin_NAME", "OpenChromatin_Evidence_Count")]
anno$CpG <- rownames(anno)

# Join in annotations to main table
d <- left_join(d, anno, by = "CpG")
write.csv(d, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/cis_trans/FULL_cistrans_table_formatted_naming_ANNOTATED.csv", row.names = F)

# Summarise ifnormation for the cis/trans associations
CIS <- d[which(d[,16] %in% "CIS"),]
dim(CIS) # 451
TRANS <- d[which(d[,16] %in% "TRANS"),] 
dim(TRANS) # 2477

table(CIS$Relation_to_Island)

 # Island N_Shelf N_Shore OpenSea S_Shelf S_Shore
 #     62      11      63     243      18      54

# proportion in either island, shore or shelf = (62 + 11 + 63 + 18 + 54) / 451 * 100 = 46%

table(TRANS$Relation_to_Island)

 # Island N_Shelf N_Shore OpenSea S_Shelf S_Shore
 #    110     108     240    1768     108     143

# proportion in either island, shore or shelf = (110 + 108 + 240 + 108 + 143) / 2477 * 100 = 29%


#############################################################################################














