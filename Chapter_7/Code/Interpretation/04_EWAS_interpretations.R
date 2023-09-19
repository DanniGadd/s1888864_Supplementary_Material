############################################################################################

### EWAS INTERPRETATION

############################################################################################

# Create annotated CpG count table with EWAS catalogue
# Plot pleiotropic signals as miami plot to show visually 

############################################################################################################

### EWAS CpG count summary table

############################################################################################################

screen 

R

library(tidyverse)
library(readxl)

# cpgs = read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/cis_trans_090621/no_eGFR_cistrans_table_unformatted_naming.csv") 
# cpgs = read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/cis_trans/FULL_cistrans_table_formatted_naming_ANNOTATED.csv")
cpgs = read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/cis_trans/FULL_cistrans_table_formatted_naming_ANNOTATED.csv") 
list1 <- unique(cpgs$SeqId) # 195 somamers

# Do a quick count of the total cpgs for each protein marker and the total cpg signatures across proteins 
EWAS_phen <- cpgs
names(EWAS_phen)[10] <- "gene"

# Of the proteome, which proteins had the most amount of CpG signals?
	Gene <- EWAS_phen[c("SeqId", "gene")]
	count <- Gene %>% group_by(SeqId) %>% count(SeqId)
	count <- as.data.frame(count)
	count <- count[order(-count$n),]
	count <- left_join(count, Gene, by = "SeqId")
	count <- unique(count)
	count <- count[c(1,3,2)]
	names(count)[3] <- "N"

count$gene <- as.character(count$gene)
count$N <- as.integer(count$N)

# Which proteins had pQTLs in Sun et al 
list <- read_excel("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_pQTL_extraction_GS/pQTLs_seqId_edited.xlsx")
list <- as.data.frame(list)
list <- list[-1,]

# 1981 pQTLs from Sun et al - for 1561 unique SeqIds 
# unique SNPs = 

count2 <- left_join(count, list, by = "SeqId")

# Work out which proteins had SNPs in our extraction (837 SNPs)
# read in SNPs extracted 
SNPs <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_pQTL_extraction_GS/chr_extraction_imputed/run_all_chr_joint_to_single_file/pQTLs_170521.csv", check.names = F)

# read in list of SNPs with SeqIds for proteins in the list for extraction
table <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Interpretation/interpretation_090621/pQTLs_in_extraction_list.csv")

# subset the table to only those SNPs we have extracted for protein regression available 
names <- colnames(SNPs)[-1]
names <- as.data.frame(names)
names(names)[1] <- "pQTL"

join <- left_join(names, table, by = "pQTL") # we have 1690 associations with biomarkers accross the 837 SNPs extracted 
names(join)[2] <- "SeqId"

# Add this info on what SNPs were available to the main table of EWAS protein counts 
count3 <- left_join(count2, join, by = "SeqId")

# Format table 
count <- count3[c(1,2,3,9,10,11,14,15,16)]
names(count) <- c("SeqId", "Gene of Protein", "Number of CpG Associations", "Sun et al Sentinel Variant", "Sun et al Chromosome", "Sun et al Position", "Sun et al Effect Allele", "Sun et al Other Allele", "pQTL extracted in STRADL")

write.csv(count, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/miami/proteins_annotated_suppl.csv")

###################

### EWAS CPG CAT ANNO

library(tidyverse)
library(readxl)

cpgs = read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/cis_trans/FULL_cistrans_table_formatted_naming_ANNOTATED.csv") 
list1 <- unique(cpgs$SeqId)

# Do a quick count of the total cpgs for each protein marker and the total cpg signatures across proteins 
EWAS_phen <- cpgs
names(EWAS_phen)[1] <- "Probe"

	# Which were the most common cpgs in the proteome EWAS signature?
	probe <- EWAS_phen[c("Probe")]
	count <- probe %>% group_by(Probe) %>% count(Probe)
	count <- as.data.frame(count)
	count <- count[order(-count$n),]
	names(count)[2] <- "No_times_selected_across_proteins"

# Read in catalogue file 
EWAS <- read.delim("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins/outputs/Glmnet_re_runs_for_all_analyses/EWAS_catalogue/EWAS_Catalog_03-07-2019.txt")

results <- count

# Set a place for the annotations 
results$Annotation <- "X"

# Get annotations for the cpgs of interest 
dim(count) # 1837 

for (i in 1:1837){
	cpg <- results[i,1]
	anno <- EWAS
	anno_cpg <- anno[which(anno$CpG %in% cpg),]
	anno_cpg <- anno_cpg[which(anno_cpg$P < 3.6e-8),]
	trait <- anno_cpg$Trait %>% unique()
	str <- str_c(trait, collapse = ", ")
	results[i,3] <- str
}

# Now we want to add a column which lists which of the proteins each cpg associates with 
names(EWAS_phen)[10] <- "Short_name"

unique(EWAS_phen$Probe) # 1837 unique probes

for (i in 1:2895){ # get dimensons of the pred results file 
	cpg <- as.character(results[i,1])
	data <- EWAS_phen
	data_cpg <- data[which(data$Probe %in% cpg),]
	list <- data_cpg$Short_name %>% unique()
	str <- str_c(list, collapse = ", ")
	results[i,"EpiScore"] <- str
}

# Add in gene of CpG and trans/cis info 
names(EWAS_phen)[3] <- "CpG_gene"
add <- EWAS_phen[,c(3,1)]

results <- left_join(results, add, by = "Probe")

# Save off results file with annotations for full EWAS 
results <- unique(results)
write.csv(results, paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/miami/cpgs_annotated_suppl.csv"), row.names = F)


##################################################################################################

### MIAMI PLOT PLEIOTROPY

screen 

R

library(tidyverse)
library(readxl)

# Fully-adjusted MWAS associations
cpgs <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/cis_trans/FULL_cistrans_table_formatted_naming_ANNOTATED.csv")

# get cpgs table ready for plotting 
cpg <- cpgs[c("Chromsome.of.CpG", "CpG.Position", "P", "CpG")]
cpg$Data <- "EWAS_v1"

# get protein table ready for plotting
prot <- cpgs[c("Chromosome.of.Gene", "Gene.Start", "P", "Gene.of.Protein")]
prot$Data <- "EWAS_v2"

names(cpg) <- c("CHR", "BP", "P", "SNP", "Data")
names(prot) <- c("CHR", "BP", "P", "SNP", "Data")

# numeric for chr of cpgs
cpg$CHR <- as.numeric(cpg$CHR)

# Rename X to chr 23 and convert to numeric 
ex <- prot[which(prot$CHR %in% "X"),] # work out how many X proteins excluded 

# > dim(ex)
# [1] 11  5
# > table(ex$SNP)

# CSAG1    F8
#     8     3

# we will present the associations for chr 1-22 as there are no CpGs on X chromosomes - only two proteins given here
prot <- prot[-which(prot$CHR %in% "X"),] 
prot$CHR <- as.numeric(prot$CHR)

# Join datasets for miami
miami_dat = rbind(cpg[,c("CHR", "BP", "SNP", "Data", "P")], prot[,c("CHR", "BP", "SNP", "Data", "P")])

### TRY GG MIAMI 
library(devtools)
install_github("juliedwhite/miamiplot")
library(miamiplot)

### Create labels for upper and lower plots separately 
# Assign labels to miami_dat file for input 
names(miami_dat)[1:3] <- c("chr", "pos", "rsid")
# Access the data being plotted
plot_data <- prep_miami_data(data = miami_dat, split_by = "Data", 
                             split_at = "EWAS_v2", p = "P")

# Create a labelling system for the top 
studyA_labels <- plot_data$upper 
studyA_labels <- subset(studyA_labels, rsid %in% c("PRG3", "PAPPA", "GZMK", "MDGA1", "MX1", "ISG15", "GBP1", "NRXN1", "COPS7B", "IL26", "CRISP2", "CRYZ", "PROK2","PCSK7",
	"LRP11", "PILRA", "HDGFRP3", "QSOX2", "GZMA"))
require(data.table)
group <- as.data.table(studyA_labels)
group_labels <- group[group[, .I[which.max(P)], by=rsid]$V1]
group_labels <- group_labels[,c("rel_pos", "logged_p", "rsid")]
names(group_labels) <- c("rel_pos", "logged_p", "label")

# Create a labelling system for the top 
studyB_labels <- plot_data$lower
studyB_labels <- subset(studyB_labels, rsid %in% c("cg07839457", "cg08122652", "cg06690548","cg16411857","cg21160290","cg22535403","cg24267699",
"cg11294350","cg11698867"))
require(data.table)
group2 <- as.data.table(studyB_labels)
group_labels2 <- group2[group2[, .I[which.max(P)], by=rsid]$V1]
group_labels2 <- group_labels2[,c("rel_pos", "logged_p", "rsid")]
names(group_labels2) <- c("rel_pos", "logged_p", "label")

# Create a highlighting system to colour key results for upper
studyA_labels <- plot_data$upper 
studyA_labels <- subset(studyA_labels, rsid %in% c("PRG3", "PAPPA", "GZMK", "MDGA1", "MX1", "ISG15", "GBP1", "NRXN1", "COPS7B", "IL26", "CRISP2", "CRYZ", "PROK2","PCSK7",
	"LRP11", "PILRA", "HDGFRP3", "QSOX2", "GZMA"))
studyA_labels$color <- "darkcyan"

# Create a highlighting system to colour key results for lower
studyB_labels <- plot_data$lower
studyB_labels <- subset(studyB_labels, rsid %in% c("cg07839457", "cg08122652", "cg06690548","cg16411857","cg21160290","cg22535403","cg24267699",
"cg11294350","cg11698867"))
studyB_labels$color <- "darkcyan"

# Miami in ggmiami 
pdf("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/miami/ggmiami_plot_labelled.pdf", width = 10, height = 5)
ggmiami(miami_dat, split_by = "Data", split_at = "EWAS_v2", suggestive_line = NULL, genome_line = NULL, chr_colors = c("darkblue", "darkblue"), chr = "chr", pos = "pos", p = "P", upper_ylab = "-log10(p)",
  lower_ylab = "-log10(p)", upper_labels_df = group_labels, lower_labels_df = group_labels2, 
  upper_highlight = studyA_labels$rsid, upper_highlight_col = "rsid", lower_highlight = studyB_labels$rsid, lower_highlight_col = "rsid",
  upper_highlight_color = studyA_labels$color, lower_highlight_color = studyB_labels$color)
dev.off()

# Write source data
write.csv(miami_dat, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Source_data/Fig3a_source_circos.csv", row.names = F)

