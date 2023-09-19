
### EWAS protein cis and trans plot 

# Read in the table that has our CpG-protein associations, with the chromosome and positional information extracted for all relevant variables

library(circlize)
library(tidyverse)
cpgs = read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/cis_trans/FULL_cistrans_table_formatted_naming_ANNOTATED.csv") # 2928
names(cpgs)[10] <- "gene"

# Take a look at the associations for the remaining proteins (580)
c <- table(cpgs$SeqId) %>% as.data.frame()
c <- c[order(-c$Freq),]

# Join in gene info for SeqIds 
gene <- cpgs[c(9,10)]
names(c)[1] <- "SeqId"
d <- left_join(c, gene, by = "SeqId")

library(dplyr)
e <- distinct(d)

# save out summary of pleiotropic associations
write.csv(e, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/cis_trans_plot/full_ewas_summary_proteins_top_second_threshold.csv")

## Subset full results to remove the 2 highly pleiotropic proteins 
cpgs <- cpgs[-which(cpgs$gene == "PRG3"),] 
dim(cpgs) # 1812

cpgs <- cpgs[-which(cpgs$gene == "PAPPA"),] 
dim(cpgs) # 825


# First see how many cis/trans unique sites we have 
names(cpgs)[16] <- "Effect"
cis = cpgs[cpgs$Effect %in% "CIS", ]
dim(cis) # 434

trans = cpgs[cpgs$Effect %in% "TRANS", ]
dim(trans) # 391


# Prep naming
names(cpgs)[3] <- "Gene_of_Hit"
cpgs[which(cpgs$Gene_of_Hit %in% ""),"Gene_of_Hit"] <- "Unannotated"
names(cpgs)[2] <- "Chromosome_of_Hit"
names(cpgs)[4] <- "Position_of_Hit"
names(cpgs)[15] <- "Chromosome_of_Somamer"
names(cpgs)[18] <- "Position_of_Somamer"

# Read in chromosome length (taken from Rob's original cis/trans plot)
gen <- cpgs
chr = read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/cis_trans_circus/CHR_Length.csv")
chr_hit <- chr
names(chr)[2] <- "Max_Chromosome_Biomarker"
names(chr_hit)[2] <- "Max_Chromosome_Hit"

gen = merge(gen,chr_hit,by.x="Chromosome_of_Hit", by.y="CHR")
gen = merge(gen,chr,by.x="Chromosome_of_Hit", by.y="CHR")

gen$Chromosome_of_Hit <- as.character(gen$Chromosome_of_Hit)
gen[gen$Chromosome_of_Hit %in% "X", "Chromosome_of_Hit"] <- "23"
gen$Chromosome_of_Somamer = as.character(gen$Chromosome_of_Somamer)
gen[gen$Chromosome_of_Somamer %in% "X", "Chromosome_of_Somamer"] <- "23"
gen$Chromosome_of_Hit = as.numeric(gen$Chromosome_of_Hit)
gen$Chromosome_of_Somamer = as.numeric(gen$Chromosome_of_Somamer)

gen$Relative_Hit_Position <- gen$Position_of_Hit/gen$Max_Chromosome_Hit
gen$Relative_Gene_Position <- gen$Position_of_Somamer/gen$Max_Chromosome_Biomarker

gen$Relative_Hit_Position <- gen$Relative_Hit_Position + gen$Chromosome_of_Hit
gen$Relative_Gene_Position <- gen$Relative_Gene_Position + gen$Chromosome_of_Somamer

# Write source data
write.csv(gen, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Source_data/Fig2c_source_cistrans.csv", row.names = F)


# Plot relative positions of hits 

library(ggplot2)
pdf("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/cis_trans_plot/cistrans_plot1_second_threshold.pdf", width = 8.9, height = 7.2)
p = ggplot(gen, aes(Relative_Hit_Position, Relative_Gene_Position)) + scale_color_manual(values = c("CIS" = "purple2", "TRANS" = "seagreen3"))
q = p + geom_jitter(aes(colour = as.factor(Effect)), size = 0.7) + xlab("CpG Position") + ylab("Protein Position")
q = q + scale_x_continuous(limits = c(1,24), breaks = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24), labels = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X", "")) + scale_y_continuous(limits = c(1,24), breaks = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24), labels = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X", "")) + theme_bw()
q = q + theme(legend.title = element_blank()) + theme(legend.text = element_text(face = "italic", size = 18)) + theme(panel.background = element_rect(colour = "black")) + geom_abline(intercept = 0, slope = 1, linetype = "dashed", size = 0.1, color = "darkslategray4")
q = q + theme(axis.text=element_text(size=14)) 
q = q + theme(axis.title.x.bottom = element_text(margin = margin(14, 0, 0, 0)))
q = q + theme(axis.title.y.left = element_text(margin = margin(0, 14, 0, 0)))
q + theme(axis.title=element_text(size=18)) 
dev.off()

