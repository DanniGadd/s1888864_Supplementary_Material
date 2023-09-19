############################################################################################
############################################################################################
################# PROTEINS OF INTEREST - STRADL ############################################
############################################################################################
############################################################################################

library(tidyverse)
library(readxl)

# Assessment of which proteins have been assessed previously in KORA/Rob studies

############################################################################################

### Read in STRADL datasets 

# Uniprot named
prot <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/prep/Phenotype_file_778.csv")
STRADL <- prot %>% as.data.frame()

############################################################################################

### Filter so no overlap with robs list of 250 for GWAS/EWAS rob has done 

############################################################################################

# Robs list of 250
rob1 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Proxies_stradl/275_Uniprots_in_STRADL_Annotated.csv")
match <- read_excel("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/rob_list_preprint_190821.xlsx")
names(match) <- c("Uniprot_ID", "SeqId")
list <- rob$Uniprot_ID

# Read in anno file with linkers from daniel / liu  
anno3 = read.csv("/Cluster_Filespace/Marioni_Group/Daniel/Somalogic_GWAS/annotation.csv", header=F)
anno3 = t(anno3)

# Get the linker file to s format that is seqIds and uniprot Ids for all STRADL proteins 
link <- anno3 %>% as.data.frame()
rownames(link) <- NULL 
link <- link[c(1,2)]
names(link) <- c("SeqId", "Uniprot_ID")

# Get the list from rob into a format where its the list of entries and uniprot IDs 
rob_link <- rob[c(1,2)]

# A little check to make sure we are picking up everything 
ov <- which(link$V2 %in% c("Q07954"))
link2 <- link[ov,]

# Match seqIds to the list from rob 
match <- left_join(rob_link, link, by = "Uniprot_ID")

# Chekc size
dim(match) # 333

# Save a copy of the rob list which has seqIds matched to it 
write.csv(match, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/prep/rob_list_indexed_to_seqIds.csv", row.names = F)

############################################################################################

### Now look at those which are also not in the KORA dataset used in their EWAS 

############################################################################################

### Work out which are not in the KORA dataset 

# Read in SomaLogic protein data in 778 with DNAm available - pre residualised for age/sex/PCs/plate
# prot1 <- readRDS("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/Validation/STRADL_778_residualised.rds")
prot <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/prep/Phenotype_file_778.csv")
names(prot)[2] <- "Stradl_id"
prot <- prot[-1] # we dont need GS id here but good to keep for reference 

stradl <- colnames(prot)[-1] 
stradl <- as.data.frame(stradl) 
names(stradl)[1] <- "UniProt"
dim(stradl) # 4235 

# Read in KORA list from EWAS study (1123)
KORA <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Proxies_stradl/KORA_proteins/Protein_info.csv")
# Take SeqId prior to _ as the seqid from KORA 
KORA$ID <- gsub("\\_.*","",KORA$SeqId)
KORA <- KORA[c(1,5,13)]

# Read in anno file with linkers from daniel / liu  
anno3 = read.csv("/Cluster_Filespace/Marioni_Group/Daniel/Somalogic_GWAS/annotation.csv", header=F)
anno3 = t(anno3)

# Take a look at how many seqIDs match between anno and KORA 
anno <- read.csv("/Cluster_Filespace/Marioni_Group/STRADL_proteins/annotation.csv", check.names = F)
names <- colnames(anno)
overlap <- which(KORA$ID %in% names) # 793
KORA2 <- KORA[overlap,]

# Save a copy
write.csv(KORA2, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/prep/KORA_list_793_indexed_to_seqIds.csv", row.names = F)


############################################################################################

### Get the list of proteins which are not in KORA or robs list 

############################################################################################

# Read in the previously generated matched file (in proxies_stradl scripts) for stradl with seqIds as colnames for reference 
seq <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Proxies_stradl/260121_Matching_cohorts/STRADL_778_residualised_matched.csv", check.names = F)

# Now read in the CORRECT protein file from the previous script (which has proper residualisation etc applied)
prot <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/prep/Phenotype_file_778.csv")
prot <- prot[-1] # can remove GS id to make sure proteins align

# As these are the same protein order (they were generated using the same raw file) we can apply column names between them 
colnames(prot)[2:4236] <- colnames(seq)[2:4236]
names(prot)[1] <- "ID"

# First remove overlap with robbie 
ov <- which(colnames(prot) %in% match$SeqId)
prot2 <- prot[,-ov]
dim(prot2) # 778 3928

# Save out the protein list which has neither KORA, nor Robs list included
write.csv(prot2, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/prep/proteins_minus_robs_3965_updated_3928_total.csv", row.names = F)

# Remove overlap with KORA
ov <- which(colnames(prot2) %in% KORA2$ID)
prot3 <- prot2[,-ov]
dim(prot3) # 3307

# Save out the protein list which has neither KORA, nor Robs list included
write.csv(prot3, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/prep/proteins_minus_robs_minus_KORA_3307_updated_3286_total.csv", row.names = F)


### The code below was used to identify proteins
### For now, as we are running on all of them, ill jump from this point straight to EWAS of the 3307 proteins 


#####################################################

# read in file with new proteins not in rob list or KORA 
proteins <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/prep/proteins_minus_robs_minus_KORA_3307_updated_3286_total.csv", check.names = F)

# read in the pQTM results file 
all <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/cis_trans/slice_neuro_final_35_pQTMs.csv")
names(all)[1] <- "Probe"

# create an assoc column 
all$assoc <- paste(all$SeqId, all$Probe, sep = "_")

# Now read in the pQTMs from Shazas paper and do similar column
library(readxl)
library(tidyverse)
zag <- read_excel("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_240521_results/formatted_tables/Zaghlool_41467_2019_13831_MOESM5_ESM_EDITED.xlsx")
zag <- as.data.frame(zag)
zag$proteinID <- sub("\\_.", "", zag$proteinID) # get seqIds for matching 
list <- zag$proteinID
zag$assoc <- paste(zag$proteinID, zag$cpg, sep = "_")

# Now read in robs pQTMs and do similar column
# rob <- read_excel("Y:/Danni/01_A_place_for_general_files/Rob_soma_paper_review/Robs_list_proteins_CpGs.xlsx")
# rob1 <- read_excel("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Interpretation/interpretation_090621/Robs_list_proteins_CpGs.xlsx")
rob <- read_excel("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/rob_bayes_pQTMs.xlsx")
rob <- as.data.frame(rob)
names(rob) <- c("CpG", "UniProt", "Gene")

# Get seqIds matched to robs proteins 
mapping <- read_excel("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Interpretation/interpretation_090621/Rob_soma_mapping.xlsx")
mapping <- as.data.frame(mapping)
names(mapping) <- c("UniProt", "SeqId")

rob <- left_join(rob, mapping, by = "UniProt")
rob$assoc <- paste(rob$SeqId, rob$CpG, sep = "_")

# Subset so that no overlapping associations remain in my results file from current lists - for the 89 neuro first to check no overlaps 
ov <- all[which(all$assoc %in% zag$assoc),] # none 
res <- all[-which(all$assoc %in% zag$assoc),] # none

# Now look at which are overlapping with robs 
overlap <- res[which(res$assoc %in% rob$assoc),] # none

# # Now get those that did not overlap:
result <- res[-which(all$assoc %in% rob$assoc),] 

### Now do the same lookup but across the fully adjusted results in general, not just neuro results 
full <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/cis_trans/FULL_cistrans_table_unformatted_naming_second_threshold.csv")
full$assoc <- paste(full$SeqId, full$Probe, sep = "_") # 2797 to subset from 

# ov2 <- full[which(full$assoc %in% rob$assoc),] # 10 
ov <- full[which(full$assoc %in% zag$assoc),] # 31 overlapped 

# 26 pQTMs overlapping with zaghlool et al
#                   assoc
# 130  4148-49_cg04983687
# 194  5124-69_cg22910295
# 205  4500-50_cg10773601
# 302  5124-69_cg03650189
# 328  5124-69_cg21994045
# 351  4148-49_cg00582671
# 389  5124-69_cg15011409
# 406  3311-27_cg07839457
# 410  4500-50_cg16651537
# 423  4148-49_cg04497992
# 515  4148-49_cg12227660
# 582  3292-75_cg07839457
# 634  4148-49_cg25087851
# 643  4148-49_cg19434937
# 682  3311-27_cg16411857
# 899  5028-59_cg07839457
# 1207 4148-49_cg08799394
# 1649 4148-49_cg18805734
# 1673 4148-49_cg18927901
# 1727 4148-49_cg12614529
# 1833 4148-49_cg21919729
# 1862 4148-49_cg23990557
# 2305 4148-49_cg16335893
# 2338 4148-49_cg23196129
# 2440 3485-28_cg07839457
# 2830  3038-9_cg07839457

sub1 <- full[-which(full$assoc %in% zag$assoc),] # 2928 to 2902

# Now look at robs overlap across all 
ov2 <- sub1[which(sub1$assoc %in% rob$assoc),] # 10 

#      gene_start  gene_end Chr     diff     diff2 Effect               assoc
# 10    113105788 113120685  13   667370 226878946    CIS  3184-25_cg19086603
# 42     60402883  60347134  15   288712 121094478    CIS 13700-10_cg27554954
# 51    117232525 117204337  11   161669 234303381    CIS  4459-68_cg07173352
# 99     20012668  19992052   2   199200  40224536    CIS 19361-78_cg24416238
# 114    60402883  60347134  15   249143 121054909    CIS 13700-10_cg12166217
# 349     3532446   3503612   4 27964503  35029395  TRANS  3640-14_cg12259379
# 1307   41163186  41158506   6 18856050 101182422  TRANS  16300-4_cg02521229
# 1512  150113372 150053291   5 93090350 207136394  TRANS 13682-47_cg07839457
# 1593   94612384  94624055  14   494543 189719311    CIS  4153-11_cg15876198
# 1995   71998613  72005715  11 14975591 129021635  TRANS  3073-51_cg07839457

sub2 <- sub1[-which(sub1$assoc %in% rob$assoc),] # 2892
# updated to 2892 novel, with 26 from shaza and durther 10 from robs - 36 total (2928 overall)

# Now check to see how many of the novel pQTMs were due to not being measured previously in Zaghlool at al
# read in 793 matches across ours and KORA study 
matches <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/prep/KORA_list_793_indexed_to_seqIds.csv")
# Subset our novel results to remove any previously measured proteins
sub3 <- sub2[-which(sub2$SeqId %in% matches$ID),] # 1783
before <- sub2[which(sub2$SeqId %in% matches$ID),] # 1109
# so 1783 of the 2928 were due to not being measured before
length(unique(sub3$SeqId)) # 150 additional seqids
length(unique(sub3$gene)) # 148 proteins 

#   [1] "MDGA1"     "EBI3"      "CTNNB1"    "LILRA3"    "GSTM3"     "B3GNT8"
#   [7] "NRXN1"     "LRP11"     "CRISP2"    "COPS7B"    "IL26"      "ICAM4"
#  [13] "QSOX2"     "CLEC4C"    "GSTM4"     "GSTZ1"     "PRG3"      "DECR2"
#  [19] "ACE"       "SHMT1"     "COL6A1"    "CTRB2"     "CHI3L1"    "KLB"
#  [25] "FAM3D"     "MX1"       "CRELD1"    "CSAG1"     "NID2"      "ERAP2"
#  [31] "CCL15"     "PILRA"     "HDGFRP3"   "ITIH3"     "GMPR"      "GSTM1"
#  [37] "TMPRSS6"   "PSAT1"     "ALPI"      "MATN3"     "HSPB1"     "LEFTY2"
#  [43] "PROK2"     "APOL3"     "SCUBE1"    "TESC"      "CRYZ"      "DAP"
#  [49] "MMAB"      "S100A7"    "BOLA1"     "GOLM1"     "LAG3"      "ANXA2"
#  [55] "DDAH1"     "ISG15"     "VNN2"      "GBP1"      "FCGR3A"    "RBL2"
#  [61] "SSU72"     "IL18R1"    "SIGLEC14"  "CECR1"     "COL2A1"    "CLPS"
#  [67] "SDF2"      "FAHD1"     "NT5C3L"    "SIGLEC5"   "ITPA"      "FGL1"
#  [73] "KLK10"     "GZMK"      "AIF1L"     "LHCGR"     "FCN2"      "OLFM2"
#  [79] "F8"        "GATM"      "SEMA4D"    "SIGLEC12"  "INHBC"     "LILRB3"
#  [85] "STXBP6"    "ACP6"      "A1BG"      "PRSS57"    "SMPD1"     "ECI2"
#  [91] "LY75"      "IGLL1"     "KIR2DS2"   "ALDOB"     "DPEP1"     "CRTAM"
#  [97] "TSTD1"     "DDX58"     "TCL1A"     "FABP2"     "B3GNT2"    "OGN"
# [103] "UGDH"      "IAPP"      "AXIN2"     "MBD1"      "RPN1"      "AOC1"
# [109] "HEXB"      "IL15RA"    "CASC4"     "PLXND1"    "SERPINA10" "CFHR1"
# [115] "TTC9B"     "APOA5"     "ACADM"     "CTSH"      "CRHBP"     "ENGASE"
# [121] "IL12B"     "ISLR2"     "IFIT3"     "SHANK3"    "GSTT1"     "CTSF"
# [127] "TLR3"      "PRG2"      "CSF1R"     "APOC3"     "CPQ"       "SIRPB1"
# [133] "AMY2A"     "LIPN"      "ADSSL1"    "STAT1"     "CD300A"    "ICAM5"
# [139] "CREG1"     "ADGRF5"    "APOA2"     "GM2A"      "A4GALT"    "PSMB1"
# [145] "DEFB1"     "RBP5"      "IL2RB"     "LRRC32"


# If we remove PRG3, how many associations that are novel and unmeasured remain 
PRG3 <- sub3[which(sub3$gene %in% "PRG3"),] # 1116
sub4 <- sub3[-which(sub3$gene %in% "PRG3"),] # 667 

# so 1109 were due to protein measured before
length(unique(before$gene)) # 41 proteins
before_sub <- before[-which(before$gene %in% "PAPPA"),] # 136 associations remaining - with 40 proteins 

#  [1] "HGFAC"    "PDK1"     "MICA"     "PCSK7"    "TPSB2"    "CLEC11A"
#  [7] "ERAP1"    "MMP1"     "CXCL6"    "LY9"      "CSF3"     "MAP2K2"
# [13] "PRSS2"    "CD163"    "FCGR2B"   "CD48"     "CST7"     "SIGLEC6"
# [19] "GP1BA"    "FCGR3B"   "GZMA"     "ACY1"     "FASLG"    "CXCL11"
# [25] "CD5L"     "EDAR"     "POR"      "RETN"     "VCAM1"    "ECM1"
# [31] "F7"       "ICAM5"    "SLITRK5"  "LGALS3"   "CST5"     "IMPDH1"
# [37] "TNFRSF1B" "MPL"      "LRPAP1"   "B2M"


