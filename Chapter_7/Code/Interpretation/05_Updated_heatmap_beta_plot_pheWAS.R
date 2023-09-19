
###################################################################################################################

### DO A CROSS-NEURO PHENOTYPE PLOT

###################################################################################################################


### COLLATE PHeWAS FILES to lists per phenotype group

screen 

R

library(tidyverse)
library(readxl)

### APOE 
apoe <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/APOE/APOE_processed/APOE_JOINT_annotated_THR.csv")
apoe$type <- "APOE"
list <- apoe[which(apoe$Status == "pass"),]
list1 <- list[c("SeqId", "phenotype", "type", "Gene.Name.Name", "beta", "SE", "Pcalc", "UniProt.Full.Name")] # 16 proteins for APOE 

### IMAGING 
brain <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/00_Results_collation/imaging/brain_accel_result_annotated_THR.csv")
brain$type <- "imaging"
brain$phenotype <- "Brain Age Acceleration"
brain <- brain[which(brain$Status == "pass"),]
list1 <- rbind(list1, brain[c("SeqId", "phenotype", "type", "Gene.Name.Name", "beta", "SE", "Pcalc", "UniProt.Full.Name")])

GGM <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/00_Results_collation/imaging/Global_GM_Volume_result_annotated_THR.csv")
GGM$type <- "imaging"
GGM$phenotype <- "Global Grey Matter Volume"
GGM <- GGM[which(GGM$Status == "pass"),]
list1 <- rbind(list1, GGM[c("SeqId", "phenotype", "type", "Gene.Name.Name", "beta", "SE", "Pcalc", "UniProt.Full.Name")])

GFA <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/00_Results_collation/imaging/gFA_result_annotated_THR.csv")
GFA$type <- "imaging"
GFA$phenotype <- "General Fractional Anisotropy"
GFA <- GFA[which(GFA$Status == "pass"),]
list1 <- rbind(list1, GFA[c("SeqId", "phenotype", "type", "Gene.Name.Name", "beta", "SE", "Pcalc", "UniProt.Full.Name")])

GMD <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/00_Results_collation/imaging/gMD_result_updated_gMD_annotated_THR.csv")
GMD$type <- "imaging"
GMD$phenotype <- "General Mean Diffusivity"
GMD <- GMD[which(GMD$Status == "pass"),]
list1 <- rbind(list1, GMD[c("SeqId", "phenotype", "type", "Gene.Name.Name", "beta", "SE", "Pcalc", "UniProt.Full.Name")])

FAZ <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/00_Results_collation/imaging/Fazekas_Score_Total_result_annotated_THR.csv")
FAZ$type <- "imaging"
FAZ$phenotype <- "Fazekas White Matter Hyperintensity Score"
FAZ <- FAZ[which(FAZ$Status == "pass"),]
list1 <- rbind(list1, FAZ[c("SeqId", "phenotype", "type", "Gene.Name.Name", "beta", "SE", "Pcalc", "UniProt.Full.Name")])

WM <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/00_Results_collation/imaging/WMH_Volume_Total_result_updated_WMHV_with_ICV_and_site_and_editor_annotated_THR.csv")
WM$type <- "imaging"
WM <- WM[which(WM$Status == "pass"),]
list1 <- rbind(list1, WM[c("SeqId", "phenotype", "type", "Gene.Name.Name", "beta", "SE", "Pcalc", "UniProt.Full.Name")])

WBV <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/00_Results_collation/imaging/WBV_No_Ventricles_result_annotated_THR.csv")
WBV$type <- "imaging"
WBV$phenotype <- "Whole Brain Volume"
WBV <- WBV[which(WBV$Status == "pass"),]
list1 <- rbind(list1, WBV[c("SeqId", "phenotype", "type", "Gene.Name.Name", "beta", "SE", "Pcalc", "UniProt.Full.Name")])


### COG SCORES
g <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/00_Results_collation/cognitive/g_result_annotated_THR.csv")
g$type <- "cognitive"
g$phenotype <- "General Cognitive Ability"
g <- g[which(g$Status == "pass"),]
list1 <- rbind(list1, g[c("SeqId", "phenotype", "type", "Gene.Name.Name", "beta", "SE", "Pcalc", "UniProt.Full.Name")])

gf <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/00_Results_collation/cognitive/gf_result_annotated_THR.csv")
gf$type <- "cognitive"
gf$phenotype <- "General Fluid Cognitive Ability"
gf <- gf[which(gf$Status == "pass"),]
list1 <- rbind(list1, gf[c("SeqId", "phenotype", "type", "Gene.Name.Name", "beta", "SE", "Pcalc", "UniProt.Full.Name")])

digit_symbol <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/00_Results_collation/cognitive/digit_symbol_result_annotated_THR.csv")
digit_symbol$type <- "cognitive"
digit_symbol$phenotype <- "Processing Speed"
digit_symbol <- digit_symbol[which(digit_symbol$Status == "pass"),]
list1 <- rbind(list1, digit_symbol[c("SeqId", "phenotype", "type", "Gene.Name.Name", "beta", "SE", "Pcalc", "UniProt.Full.Name")])

LM <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/00_Results_collation/cognitive/LM_result_annotated_THR.csv")
LM$type <- "cognitive"
LM$phenotype <- "Logical Memory"
LM <- LM[which(LM$Status == "pass"),]
list1 <- rbind(list1, LM[c("SeqId", "phenotype", "type", "Gene.Name.Name", "beta", "SE", "Pcalc", "UniProt.Full.Name")])

mr_correct <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/00_Results_collation/cognitive/mr_correct_result_annotated_THR.csv")
mr_correct$type <- "cognitive"
mr_correct$phenotype <- "Non-Verbal Reasoning"
mr_correct <- mr_correct[which(mr_correct$Status == "pass"),]
list1 <- rbind(list1, mr_correct[c("SeqId", "phenotype", "type", "Gene.Name.Name", "beta", "SE", "Pcalc", "UniProt.Full.Name")])

verbal_total <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/00_Results_collation/cognitive/verbal_total_result_annotated_THR.csv")
verbal_total$type <- "cognitive"
verbal_total$phenotype <- "Verbal Reasoning"
verbal_total <- verbal_total[which(verbal_total$Status == "pass"),]
list1 <- rbind(list1, verbal_total[c("SeqId", "phenotype", "type", "Gene.Name.Name", "beta", "SE", "Pcalc", "UniProt.Full.Name")])

vocabulary <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/00_Results_collation/cognitive/vocabulary_result_annotated_THR.csv")
vocabulary$type <- "cognitive"
vocabulary$phenotype <- "Vocabulary"
vocabulary <- vocabulary[which(vocabulary$Status == "pass"),]
list1 <- rbind(list1, vocabulary[c("SeqId", "phenotype", "type", "Gene.Name.Name", "beta", "SE", "Pcalc", "UniProt.Full.Name")])


# Work out how many assocs for each phenotype and modality 
# apoe is already 11 as seen above 

dim(list1) # 405 total associations for phenos 
unique(list1[,4]) # 191 proteins unique to these associations

apoel <- list1 %>% filter(type == "APOE")
dim(apoel) # 14 APOE 

cog <- list1 %>% filter(type == "cognitive")
dim(cog) # 296 cognitive 

im <- list1 %>% filter(type == "imaging")
dim(im) # 95 imaging


## Next, read the files in and chop down to key info for plotting 
apoe <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/APOE/APOE_processed/APOE_JOINT_annotated_THR.csv")
brain <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/00_Results_collation/imaging/brain_accel_result_annotated_THR.csv")
GGM <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/00_Results_collation/imaging/Global_GM_Volume_result_annotated_THR.csv")
GFA <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/00_Results_collation/imaging/gFA_result_annotated_THR.csv")
GMD <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/00_Results_collation/imaging/gMD_result_updated_gMD_annotated_THR.csv")
FAZ <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/00_Results_collation/imaging/Fazekas_Score_Total_result_annotated_THR.csv")
WM <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/00_Results_collation/imaging/WMH_Volume_Total_result_updated_WMHV_with_ICV_and_site_and_editor_annotated_THR.csv")
WBV <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/00_Results_collation/imaging/WBV_No_Ventricles_result_annotated_THR.csv")
g <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/00_Results_collation/cognitive/g_result_annotated_THR.csv")
gf <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/00_Results_collation/cognitive/gf_result_annotated_THR.csv")
digit_symbol <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/00_Results_collation/cognitive/digit_symbol_result_annotated_THR.csv")
LM <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/00_Results_collation/cognitive/LM_result_annotated_THR.csv")
mr_correct <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/00_Results_collation/cognitive/mr_correct_result_annotated_THR.csv")
verbal_total <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/00_Results_collation/cognitive/verbal_total_result_annotated_THR.csv")
vocabulary <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/00_Results_collation/cognitive/vocabulary_result_annotated_THR.csv")

## Plot for APOE, cognitive and brain imaging separately 
library(ggplot2)

# Cognitive plot 

g <- g[c(1,2,3,5,7)]
g$phenotype <- "General Cognitive Ability"
gf <- gf[c(1,2,3,5,7)]
gf$phenotype <- "General Fluid Cognitive Ability"
digit_symbol <- digit_symbol[c(1,2,3,5,7)]
digit_symbol$phenotype <- "Processing Speed"
LM <- LM[c(1,2,3,5,7)]
LM$phenotype <- "Logical Memory"
mr_correct <- mr_correct[c(1,2,3,5,7)]
mr_correct$phenotype <- "Non-Verbal Reasoning"
verbal_total <- verbal_total[c(1,2,3,5,7)]
verbal_total$phenotype <- "Verbal Reasoning"
vocabulary <- vocabulary[c(1,2,3,5,7)]
vocabulary$phenotype <- "Vocabulary"

bind <- rbind(verbal_total, vocabulary)
bind <- rbind(bind, LM)
bind <- rbind(bind, mr_correct)
bind <- rbind(bind, digit_symbol)
bind <- rbind(bind, gf)
bind <- rbind(bind, g)

# Restrict to the proteins of interest across all phenotypes from above 
bind <- bind[which(bind$SeqId %in% cog$SeqId),] # list1 has 497 associations in protein PheWAS

# Add stars for significance (at threshold)
bind$stars <- cut(bind$Pcalc, breaks=c(-Inf, 0.0003496503, 0.01, 0.05, Inf), label=c("*", "", "", ""))

# Make genes unique 
names(bind)[2] <- "gene"

# Assign to variable for stitching 
cogpplot <- ggplot(aes(x=gene, y=phenotype, fill=beta), data=bind) +
 geom_tile() + scale_fill_gradient2(low="#2C7BB6", mid="white", high="#D7191C") + 
 geom_text(aes(label=stars), color="black", size=5) +
theme(legend.position = 'None', axis.text.x = element_text(size=8, angle=45, hjust = 0.9),
          axis.text.y = element_text(size=10, angle=0)) + xlab("") + ylab("")  + guides(fill=guide_legend(title="Beta effect")) +
scale_y_discrete(limits = c("Verbal Reasoning", "Vocabulary", "Logical Memory", "Non-Verbal Reasoning", 
  "General Fluid Cognitive Ability", "General Cognitive Ability", "Processing Speed"))

# APOE 

apoe <- apoe[c(1,2,3,5,7)]
apoe$phenotype <- "APOE"
bind <- apoe

# Restrict to the proteins of interest across all phenotypes from above 
bind <- bind[which(bind$SeqId %in% apoel$SeqId),] # list1 has 497 associations in protein PheWAS

# Add stars for significance (at threshold)
bind$stars <- cut(bind$Pcalc, breaks=c(-Inf, 0.0003496503, 0.01, 0.05, Inf), label=c("*", "", "", ""))

# Make genes unique 
names(bind)[2] <- "gene"

# Assign 
apoeplot <- ggplot(aes(x=gene, y=phenotype, fill=beta), data=bind) +
 geom_tile() + scale_fill_gradient2(low="#2C7BB6", mid="white", high="#D7191C") + 
 geom_text(aes(label=stars), color="black", size=5) +
theme(legend.position = 'None', axis.text.x = element_text(size=8, angle=45, hjust = 0.9),
          axis.text.y = element_text(size=10, angle=0)) + xlab("") + ylab("")  + guides(fill=guide_legend(title="Beta effect")) +
scale_y_discrete(limits = c("APOE"))


# Imaging
brain <- brain[c(1,2,3,5,7)]
brain$phenotype <- "Relative Brain Age"
GGM <- GGM[c(1,2,3,5,7)]
GGM$phenotype <- "Global Grey Matter Volume"
GFA <- GFA[c(1,2,3,5,7)]
GFA$phenotype <- "General Fractional Anisotropy"
GMD <- GMD[c(1,2,3,5,7)]
GMD$phenotype <- "General Mean Diffusivity"
FAZ <- FAZ[c(1,2,3,5,7)]
FAZ$phenotype <- "Fazekas White Matter Hyperintensity Score"
WM <- WM[c(1,2,3,5,7)]
WM$phenotype <- "White Matter Hyperintensity Volume"
WBV <- WBV[c(1,2,3,5,7)]
WBV$phenotype <- "Whole Brain Volume"

bind <- brain
bind <- rbind(bind, GGM)
bind <- rbind(bind, GFA)
bind <- rbind(bind, GMD)
bind <- rbind(bind, WM)
bind <- rbind(bind, FAZ)
bind <- rbind(bind, WBV)

# Restrict to the proteins of interest across all phenotypes from above 
bind <- bind[which(bind$SeqId %in% im$SeqId),] # list1 has 497 associations in protein PheWAS

# Add stars for significance (at threshold)
bind$stars <- cut(bind$Pcalc, breaks=c(-Inf, 0.0003496503, 0.01, 0.05, Inf), label=c("*", "", "", ""))

# Make genes unique 
names(bind)[2] <- "gene"

# Assign
implot <- ggplot(aes(x=gene, y=phenotype, fill=beta), data=bind) +
 geom_tile() + scale_fill_gradient2(low="#2C7BB6", mid="white", high="#D7191C") + 
 geom_text(aes(label=stars), color="black", size=5) +
theme(legend.position = 'None', axis.text.x = element_text(size=8, angle=45, hjust = 0.9),
          axis.text.y = element_text(size=10, angle=0)) + xlab("") + ylab("")  + guides(fill=guide_legend(title="Beta effect")) +
scale_y_discrete(limits = c("Fazekas White Matter Hyperintensity Score", "General Mean Diffusivity", "Whole Brain Volume", "Global Grey Matter Volume", 
  "General Fractional Anisotropy", "Relative Brain Age", "White Matter Hyperintensity Volume"))

### PATCHWORK 
library(patchwork)
pdf("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/heatmap_phewas/heatmap_thresholded_joint_individual_plots.pdf", width = 25, height = 9)
implot / cogpplot / apoeplot +  plot_layout(heights = unit(c(7, 7, 1), c('cm', 'cm', 'cm')))
dev.off()

###################################################################################################

### Do version with common proteins across imaging and cognitive 

###################################################################################################

screen 

R

## Next, read the files in and chop down to key info for plotting 
apoe <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/APOE/APOE_processed/APOE_JOINT_annotated_THR.csv")
brain <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/00_Results_collation/imaging/brain_accel_result_annotated_THR.csv")
GGM <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/00_Results_collation/imaging/Global_GM_Volume_result_annotated_THR.csv")
GFA <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/00_Results_collation/imaging/gFA_result_annotated_THR.csv")
GMD <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/00_Results_collation/imaging/gMD_result_updated_gMD_annotated_THR.csv")
FAZ <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/00_Results_collation/imaging/Fazekas_Score_Total_result_annotated_THR.csv")
WM <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/00_Results_collation/imaging/WMH_Volume_Total_result_updated_WMHV_with_ICV_and_site_and_editor_annotated_THR.csv")
WBV <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/00_Results_collation/imaging/WBV_No_Ventricles_result_annotated_THR.csv")
g <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/00_Results_collation/cognitive/g_result_annotated_THR.csv")
gf <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/00_Results_collation/cognitive/gf_result_annotated_THR.csv")
digit_symbol <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/00_Results_collation/cognitive/digit_symbol_result_annotated_THR.csv")
LM <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/00_Results_collation/cognitive/LM_result_annotated_THR.csv")
mr_correct <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/00_Results_collation/cognitive/mr_correct_result_annotated_THR.csv")
verbal_total <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/00_Results_collation/cognitive/verbal_total_result_annotated_THR.csv")
vocabulary <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/00_Results_collation/cognitive/vocabulary_result_annotated_THR.csv")

## Plot for APOE, cognitive and brain imaging separately 
library(ggplot2)

# Cognitive plot 

g <- g[c(1,2,3,5,7)]
g$phenotype <- "General Cognitive Ability"
gf <- gf[c(1,2,3,5,7)]
gf$phenotype <- "General Fluid Cognitive Ability"
digit_symbol <- digit_symbol[c(1,2,3,5,7)]
digit_symbol$phenotype <- "Processing Speed"
LM <- LM[c(1,2,3,5,7)]
LM$phenotype <- "Logical Memory"
mr_correct <- mr_correct[c(1,2,3,5,7)]
mr_correct$phenotype <- "Non-Verbal Reasoning"
verbal_total <- verbal_total[c(1,2,3,5,7)]
verbal_total$phenotype <- "Verbal Reasoning"
vocabulary <- vocabulary[c(1,2,3,5,7)]
vocabulary$phenotype <- "Vocabulary"

bind <- rbind(verbal_total, vocabulary)
bind <- rbind(bind, LM)
bind <- rbind(bind, mr_correct)
bind <- rbind(bind, digit_symbol)
bind <- rbind(bind, gf)
bind <- rbind(bind, g)

apoe <- apoe[c(1,2,3,5,7)]
apoe$phenotype <- "APOE haplotype"

bind <- rbind(bind, apoe)

brain <- brain[c(1,2,3,5,7)]
brain$phenotype <- "Relative Brain Age"
GGM <- GGM[c(1,2,3,5,7)]
GGM$phenotype <- "Global Grey Matter Volume"
GFA <- GFA[c(1,2,3,5,7)]
GFA$phenotype <- "General Fractional Anisotropy"
GMD <- GMD[c(1,2,3,5,7)]
GMD$phenotype <- "General Mean Diffusivity"
FAZ <- FAZ[c(1,2,3,5,7)]
FAZ$phenotype <- "Fazekas White Matter Hyperintensity Score"
WM <- WM[c(1,2,3,5,7)]
WM$phenotype <- "White Matter Hyperintensity Volume"
WBV <- WBV[c(1,2,3,5,7)]
WBV$phenotype <- "Whole Brain Volume"

bind <- rbind(bind, brain)
bind <- rbind(bind, GGM)
bind <- rbind(bind, GFA)
bind <- rbind(bind, GMD)
bind <- rbind(bind, WM)
bind <- rbind(bind, FAZ)
bind <- rbind(bind, WBV)

# Read in the common associations table 
common <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/common_proteins/common_proteins.csv")
subset <- common$SeqId %>% as.data.frame()
listZ <- c("17671-58", "2797-56", "4337-49") %>% as.data.frame()
subset <- rbind(subset, listZ)
# Subset bind file to SeqIds in the common assocs file 
bind <- bind[which(bind$SeqId %in% subset[,1]),]

# Index the common associations 
bind$stars <- cut(bind$Pcalc, breaks=c(-Inf, 0.0003496503, 0.01, 0.05, Inf), label=c("*", "", "", ""))

# Make genes unique 
names(bind)[2] <- "gene"

length(unique(bind$SeqId)) # 26 somamers 
length(unique(bind$gene)) # 25 proteins (22 cog/imaging and 3 APOE/cog) 

# Write source data
write.csv(bind, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Source_data/Fig4b_source_heatmap.csv", row.names = F)


# Plot 
library(ggplot2)

pdf("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/heatmap_phewas/common_somamers.pdf", width = 14.5, height = 6)
ggplot(aes(x=gene, y=phenotype, fill=beta), data=bind) +
 geom_tile() + scale_fill_gradient2(low="#2C7BB6", mid="white", high="#D7191C") + 
 geom_text(aes(label=stars), color="black", size=7) +
theme(legend.position = 'None', axis.text.x = element_text(size=14, angle=45, hjust = 0.9),
          axis.text.y = element_text(size=16, angle=0)) + xlab("") + ylab("") + guides(fill=guide_legend(title="Beta effect")) +
scale_y_discrete(limits = c("Fazekas White Matter Hyperintensity Score", "Whole Brain Volume", "General Mean Diffusivity",
     "Verbal Reasoning", "Logical Memory", "Vocabulary", "APOE haplotype", "Global Grey Matter Volume", "White Matter Hyperintensity Volume", 
     "General Fractional Anisotropy", "Relative Brain Age", "Non-Verbal Reasoning", 
  "General Fluid Cognitive Ability", "General Cognitive Ability", "Processing Speed")) 
dev.off()

# Count how many asterisks for associations
table(bind$stars) # 77

####################################

### CREATE PLOT SHOWING TOTAL NUMBER OF ASSOCIATIONS PER 15 TRAITS 

screen 

R

library(MetBrewer)
library(readxl)
library(ggplot2)
library(ggthemes)
library(readxl)
library(tidyverse)

table <- read_excel("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/heatmap_phewas/phewas_table_counts.xlsx")
table <- as.data.frame(table)
table <- table[c(1:3)]
names(table) <- c("Phenotype", "Marker", "Count")


#                                          Phenotype           Marker Count
# 1                               Relative Brain Age   Imaging marker    24
# 2                    General Fractional Anisotropy   Imaging marker    22
# 3               White Matter Hyperintensity Volume   Imaging marker    19
# 4                        Global Grey Matter Volume   Imaging marker    15
# 5                               Whole Brain Volume   Imaging marker     6
# 6                         General Mean Diffusivity   Imaging marker     6
# 7                              Fazekas score (WMH)   Imaging marker     3
# 8                                   APOE haplotype   APOE haplotype    14
# 9                  Processing speed - Digit Symbol Cognitive marker   102
# 10                       General Cognitive Ability Cognitive marker    73
# 11                 General Fluid Cognitive Ability Cognitive marker    54
# 12 Non-Verbal Reasoning - Matrix Reasining Correct Cognitive marker    38
# 13                                Vocabulary Score Cognitive marker    13
# 15                                  Logical Memory Cognitive marker     9
# 14                 Verbal Reasoning - Verbal Total Cognitive marker     7


table$Phenotype <- c("Relative Brain Age", "General Fractional Anisotropy", "White Matter Hyperintensity Volume", 
  "Global Grey Matter Volume", "Whole Brain Volume", "General Mean Diffusivity", "Fazekas White Matter Hyperintensity Score", "APOE status", 
  "Processing Speed", "General Cognitive Ability", "General Fluid Cognitive Ability", "Non-Verbal Reasoning", "Vocabulary", "Verbal Reasoning",
  "Logical Memory")

plot_data <- table %>% 
  group_by(Phenotype) %>% 
  arrange(Count) %>% 
  ungroup() %>% 
  mutate(order = row_number())

# Write source data
write.csv(plot_data, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Source_data/Fig4a_source_barplot.csv", row.names = F)


pdf("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/heatmap_phewas/plot_bar_V2.pdf", width = 14.5, height = 6)
ggplot(plot_data, aes(order, Count, fill = Marker)) +
  geom_col() +
  scale_x_continuous(breaks = plot_data$order, labels = plot_data$Phenotype) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  coord_flip() +
scale_fill_manual(breaks = c("Imaging marker", "Cognitive marker", "APOE haplotype"),
values = c("cadetblue1", "lightseagreen", "royalblue3")) + theme_minimal() +
theme(panel.grid.major.y = element_blank()) + theme(legend.title = element_blank()) + theme(legend.position = 'None') +
  labs(y = "Number of associations with protein levels",
       x = "") +
  theme(axis.text.x = element_text(size=16),
          axis.text.y = element_text(size=16, angle=0), axis.title = element_text(size=16)) 
dev.off()

