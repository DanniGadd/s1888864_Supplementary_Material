### PLOT FOR COMMON PROTEINS PCA 

cd /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/common_proteins/

screen

R

library(tidyverse)
library(ggplot2)
library(readxl)
library(psych)
library(ggcorrplot)
library(cowplot)
library(tidyverse)

# Read in the protein file 
prot <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS/01_Phenotype_collation/prot_file_220221.csv", check.names = F)

# Read in "common" file from the updated beta plots script for the plot with 25 proteins 
# List the top 25 proteins (26 somamers)
common <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/common_proteins/common_proteins_V2.csv")
list <- unique(common$SeqId)

# subset prot file to include just these somamers 
prot2 <- prot[,which(colnames(prot) %in% list)]

# Assign variable to joint 
joint <- prot2

# # 4876-32 - F9 
# # 5307-12 - F9 

# joint <- joint[,-which(colnames(joint) %in% "5307-12")]

# Scale data and get PC scores
scores_pca <- principal(scale(joint), rotate="none", nfactors=ncol(joint))$scores
variance_pca <- principal(scale(joint), rotate="none", nfactors=ncol(joint))$Vaccounted
pca <- principal(scale(joint), rotate="none", nfactors=ncol(joint))

scores <- as.data.frame(scores_pca)
var1 <- as.data.frame(variance_pca)

var <- var1[1,] # get variance 
cum <- var1[5,] # get cumulative variance

var <- gather(var) # gather
var$col <- ifelse(var$value >= 1, "darkgrey", "orange")
var$num <- as.integer(1:ncol(joint))
var$mes <- as.numeric(var$value)

cum <- gather(cum) # gather
cum$num <- as.integer(1:ncol(joint))
cum$mes <- as.numeric(cum$value)

cor = cor(joint)

## Read in seq-id conversion file 
anno <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Annotations_protein_file/annotation_formatted_for_paper.csv")

## subset seq-ids
anno1 = anno[which(anno$SeqId %in% colnames(cor)),] 
## match up their order 
ids = colnames(cor)[which(colnames(cor) %in% anno1$SeqId)] 
anno1 = anno1[match(ids, anno1$SeqId),]

## check they match
table(anno1$SeqId == colnames(cor)[which(colnames(cor) %in% anno1$SeqId)] )

## replace seq-ids with gene names 
names(anno1)[3] <- "Name"
colnames(cor)[which(colnames(cor) %in% anno1$SeqId)] <- as.character(anno1$Name)
row.names(cor)[which(row.names(cor) %in% anno1$SeqId)] <- as.character(anno1$Name)

# a <- if(ncol(joint) > 12) { 
a <- ggplot(var, aes(num, mes, col = col)) +
geom_point(size = 4) + theme(legend.position = "none", axis.text=element_text(size=13), axis.title=element_text(size=13)) +
xlab("Principal component") + ylab("Eigenvalue") + geom_hline(yintercept=1, color = "darkgrey", size=1) + 
scale_color_manual(values=c("#E69F00", "#999999")) 


b <- ggplot(cum, aes(num, mes)) +
geom_bar(stat = "identity", fill = "steelblue") + theme(legend.position = "none", axis.text=element_text(size=13), axis.title=element_text(size=13)) +
xlab("Principal component") + ylab("Cumulative proportion") + ylim(0,1) 

c <- ggcorrplot(cor, 
           hc.order = TRUE,
           type = "lower")


c <- c + theme(
  panel.background = element_rect(fill = "white", colour = "white",
                                size = 2, linetype = "solid"),
  panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                colour = "white"), 
  panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                colour = "white")
  )


p1 = plot_grid(c,a,b, nrow = 1, labels = c("a", "b", "c"), rel_widths = c(0.85,0.5,0.5))


pdf("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/common_proteins/FIGURE_lower_22_proteins.pdf", height =6 , width = 17)
p1
dev.off()





