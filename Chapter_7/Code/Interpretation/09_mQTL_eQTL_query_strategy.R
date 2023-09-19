##################################################################################################

### mQTL eQTL lookup 

##################################################################################################

screen

R

### Assess which CpGs were unique to EPIC array 

# Read in cpgs results file 
cpgs <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/cis_trans/slice_neuro_final_35_pQTMs_pQTLs_added.csv")

# Read in EPIC annotation file
anno <- readRDS("/Cluster_Filespace/Marioni_Group/Daniel/EPIC_AnnotationObject_df.rds")

# Subset annotation file to probes common to 450k and EPIC array
common_anno <- anno[which(anno$Methyl450_Loci == "TRUE"),]

# index which CpGs were in the 450K
length(which(cpgs$CpG %in% common_anno$Name)) # 19 are in 450K

# Create version that is additional to the 450k array 
anno_extra <- anno[-which(anno$Methyl450_Loci == "TRUE"),]

# > dim(anno)
# [1] 866836     46
# > dim(common_anno)
# [1] 453093     46
# > dim(anno_extra)
# [1] 413743     46

# Index whether any CpGs were EPIC only 
length(which(cpgs$CpG %in% anno_extra$Name)) # 16 are in EPIC 

# Save names of each set 
cpgs_450 <- cpgs[which(cpgs$CpG %in% common_anno$Name),]
cpgs_450$array <- "450K"

cpgs_EPIC <- cpgs[which(cpgs$CpG %in% anno_extra$Name),]
cpgs_EPIC$array <- "EPIC"


cpgs_450 <- cpgs_450[c(1,19)]
cpgs_EPIC <- cpgs_EPIC[c(1,19)]

# find unique cpgs in each 
unique(cpgs_450$CpG) # 16
unique(cpgs_EPIC$CpG) # 15


## Load requisite packages 
library(data.table)

## Read in CpG list 
cpgs1 = read.csv("/Cluster_Filespace/Marioni_Group/Rob/Somalogic/EWAS_Hits_Final_Dec.csv")

#  [1] "Gene_of_Somamer"       "Biomarker"             "Marker.x"
#  [4] "chr"                   "pos"                   "UCSC_RefGene_Name"
#  [7] "UCSC_RefGene_Group"    "Relation_to_Island"    "PIP"
# [10] "Biomarker.1"           "Beta"                  "Beta_LCI"
# [13] "Beta_HCI"              "Overlap"               "Chromosome_of_Somamer"
# [16] "Position_of_Somamer"   "Distance"              "Type"

cpgs <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/cis_trans/slice_neuro_final_35_pQTMs_pQTLs_added.csv")
# cpgs <- cpgs[c(10,11,1,2,4,3,5:9,12:18)]
names(cpgs)[10] <- "Gene_of_Somamer"
names(cpgs)[11] <- "Biomarker"
names(cpgs)[1] <- "Marker.x"

## Read in GoDMC mQTLs (at P < 5e-8)
mqtls = as.data.frame(fread("/Cluster_Filespace/Marioni_Group/Rob/GoDMC_mqtls_genomewide.csv"))

## Loop to extract mQTLs for CpGs 
## Set up output list 
out <- list() 

## Loop step 
for(i in 1:nrow(cpgs)){ 
  ## Extract CpG to query 
  tmp = cpgs[i,]
  ## Subset GoDMC to queried CpG - finding mQTLs for it 
  mqtls.tmp = mqtls[which(mqtls$cpg %in% tmp$Marker.x),]
  ## Store output 
  out[[i]] <- mqtls.tmp 
}
## collate mqtls 
out1 = as.data.frame(do.call("rbind", out))

## As a next step, we will format SNPs to chr:pos - helps for matching to eqtls etc.
## first extract chr and pos of mqtls as separate columns 
## CHR
out1$CHR <- sapply(strsplit(out1$snp, ":"), "[", 1)
out1$CHR <- as.numeric(gsub("chr", "", out1$CHR))
## POS 
out1$POS <- sapply(strsplit(out1$snp, ":"), "[", 2)
## then combine chr:pos 
out1$MarkerName <- paste(out1$CHR, out1$POS, sep = ":")

## annotate original file to denote mQTL status 
cpgs$mQTL <- "No"
cpgs[which(cpgs$Marker.x %in% out1$cpg),"mQTL"] <- "Yes"

## check how many have mQTLs - note, some cpgs are present across >1 protein 
table(cpgs$mQTL)

## double check number of unique CpGs match
table(unique(out1$cpg) == unique(cpgs[cpgs$mQTL %in% "Yes", "Marker.x"])) # TRUE

## Convert to rsid for clumping step 
# remotes::install_github("MRCIEU/ieugwasr")
library(ieugwasr)

## Read in rsid information 
anno = as.data.frame(fread("/Cluster_Filespace/Marioni_Group/Daniel/GWAS_AgeAccel/Cohort_Summary_Stats/EasyQC/Allele_Freq_and_Mapping_Info/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.rsid_map.gz"))
anno$MarkerName = paste0(anno$chr, ":", anno$pos)
## obtain rsids
out1$rsid = anno[match(out1$MarkerName, anno$MarkerName), "rsid"]
## clumping step 
## make CpG 'id' as this will force clumping to be done separately for each CpG 
names(out1)[1] <- "id"
## rename P 
out2 <- ld_clump(
  dat = out1,
  clump_kb = 10000,
  clump_r2 = 0.001,
  clump_p = 1,
  pop = "EUR",
  access_token = NULL,
  bfile = NULL,
  plink_bin = NULL
)

## convert 'id' back to CpG 
names(out2)[1] <- "cpg"
## merge with original CpG file to find which proteins the CpGs associated with (as we will use this in the eQTL step)
query = merge(out2, cpgs[,c("Gene_of_Somamer", "Marker.x")], by.x = "cpg", by.y = "Marker.x")

## some cpgs were associated with more than one protein - so query df is larger than out2 df
## we will ask whether the mQTLs for those CpGs, are also eQTLs for the protein (gene) that the CpG associated with

## Query whether those independent mQTLs are eQTLs 
## Read in eQTL Gen data 
eqtls = readRDS("/Cluster_Filespace/Marioni_Group/Rob/Neurology Proteins/eqtls/eQTLgen_data.rds")

# ## Convert summ stats to Beta and SE - would  need for coloc
# eqtls$BETA <- eqtls$Zscore/sqrt((2*eqtls$MAF)*(1-eqtls$MAF)*(eqtls$NrSamples + eqtls$Zscore^2))
# eqtls$SE <- 1/sqrt((2*eqtls$MAF)*(1-eqtls$MAF)*(eqtls$NrSamples + eqtls$Zscore^2))

## set up loop to ask whether mQTLs are also known eQTLs for gene that CpG associates with 
output <- as.data.frame(matrix(nrow = nrow(query), ncol = 4))
names(output) <- c("cpg", "mQTL", "gene", "eQTL_status")

for(i in 1:nrow(query)){ 
tmp = query[i,]
## extract mQTL for query (is it an eQTL for gene of interest)
snp = tmp[,"rsid"]
## extract gene of interest 
gene = tmp[,"Gene_of_Somamer"]
## extract cpg - for output table 
cpg = tmp[,"cpg"]
## Main step 
## subset eQTLGen to gene first (see if its there)
eqtl.tmp = eqtls[which(eqtls$GeneSymbol %in% gene),]
## if its not there, skip 
if(nrow(eqtl.tmp) == 0){ 
  output[i,1] <- as.character(cpg)
  output[i,2] <- as.character(snp)
  output[i,3] <- as.character(gene)
  output[i,4] <- "No"
  } else { 
## if it is there, we will further test whether mQTL is also an eQTL 
eqtl.tmp1 = eqtl.tmp[eqtl.tmp$SNP %in% snp,]
## if it is not there, eqtl.tmp1 will have no rows
if(nrow(eqtl.tmp1) == 0){ 
  output[i,1] <- as.character(cpg)
  output[i,2] <- as.character(snp)
  output[i,3] <- as.character(gene)
  output[i,4] <- "No"
} ## if it not there, eqtl.tmp1 will have entries 
else { 
  output[i,1] <- as.character(cpg)
  output[i,2] <- as.character(snp)
  output[i,3] <- as.character(gene)
  output[i,4] <- "Yes"}
  }
print(paste("Completing", paste0(as.character(i), "/", nrow(query))))
} 

## Subset to final mQTLs and eQTLs 
output.eqtls <- output[output$eQTL_status %in% "Yes", ]

# Wrtie out tables generated
write.csv(cpgs, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/mQTL_eQTL_lookup/mQTL_lookup_neuro_table_V2.csv", row.names = F)
write.csv(output.eqtls, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/mQTL_eQTL_lookup/mQTL_eQTL_lookup_V2.csv", row.names = F)

# Join eQTLs into main table
library(tidyverse)
cpgs_anno <- merge(cpgs, output.eqtls, by.x = "Marker.x", by.y = "cpg", all = T)
cpgs_anno <- cpgs_anno[order(cpgs_anno$P),]
write.csv(cpgs_anno, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/mQTL_eQTL_lookup/main_table_mQTL_eQTL_lookup_V2.csv", row.names = F)

# # Match order of cpgs to results file 
# cpgs_anno2 <- cpgs_anno[match(cpgs$Marker.x, cpgs_anno$Marker.x),]
# write.csv(cpgs_anno2, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/mQTL_eQTL_lookup/main_table_mQTL_eQTL_lookup_matched_subset_V2.csv", row.names = F)

# Write out query table with fullmQTLs 
write.csv(query, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/mQTL_eQTL_lookup/query.csv", row.names = F)


#######################

### ADD CpG INFORMATION TO THE SUMMARY TABLE FOR NEURO ASSOCIATIONS

# Read in the UCSC cpg annotations file 
anno <- readRDS("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Annotations_UCSC/annotations_UCSC_object_df_from_daniel.rds")

# Read in the table with mQTL/eQTL annotations 
neuro <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/mQTL_eQTL_lookup/main_table_mQTL_eQTL_lookup_V2.csv")
names(neuro)[1] <- "CpG"

# Subset annotations file to cpgs involved in associations
anno <-  anno[which(rownames(anno) %in% neuro$CpG),]
dim(anno) # 31 cpgs 

# Select columns of interest from the annotation file 
anno <- anno[,c("chr", "pos", "strand", "UCSC_RefGene_Group", "Islands_Name", "UCSC_RefGene_Name", 
    "Relation_to_Island", "OpenChromatin_NAME", "OpenChromatin_Evidence_Count")]
anno$CpG <- rownames(anno)

# Join in annotations to main table
neuro <- left_join(neuro, anno, by = "CpG")

write.csv(neuro, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/mQTL_eQTL_lookup/main_table_mQTL_eQTL_lookup_matched_subset_ANNO_ISLANDS.csv")

# Examine locational information

table(neuro$Relation_to_Island)

# N_Shelf N_Shore OpenSea S_Shelf S_Shore
#       1       4      33       1       1

table(neuro$UCSC_RefGene_Group)

#                                                      1
#                                                1stExon
#                                                      2
#                                          1stExon;5'UTR
#                                                      1
#                                                  3'UTR
#                                                      1
#                                                  5'UTR
#                                                      2
#                                      5'UTR;5'UTR;5'UTR
#                                                      1
#                                             5'UTR;Body
#                                                      1
#                                                   Body
#                                                     11
#                                              Body;Body
#                                                      3
#                                    Body;Body;Body;Body
#                                                      1
# Body;Body;Body;Body;Body;Body;Body;Body;Body;Body;Body
#                                                      1
#                                                TSS1500
#                                                      8
#                                                 TSS200
#                                                      7
