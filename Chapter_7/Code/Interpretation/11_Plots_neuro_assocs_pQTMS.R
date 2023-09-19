
########################################################################

### READ IN NEURO CPGS AND PLOT ASSOCIATIONS WITH PROTEINS

screen 

R

slice <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/cis_trans/slice_neuro_final_41_pQTMs_pQTLs_added.csv")

prot <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS/01_Phenotype_collation/prot_file_270121.csv", check.names = F)


# Join protei data to methylation data 
prot$GS_id <- as.character(prot$GS_id)

neuro <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/plots_neuro_assocs/plots/neuro_joint_protein_cpgs.csv", check.names = F)

# Plot associations of interest 
table <- slice[c("CpG", "SeqId", "Gene.of.Protein")]

library(ggpubr)
library(ggplot2)

plots <- list()

for(i in 1:length(table$CpG)){
    CpG <- as.character(table[i,1])
    protein <- as.character(table[i,2])
    gene <- as.character(table[i,3])
    data <- neuro[,CpG] %>% as.data.frame()
    data2 <- neuro[,protein] %>% as.data.frame()
    dataset <- cbind(data,data2)
    names(dataset) <- c("CpG", "Protein")
     p <- ggplot(dataset, aes(x=CpG,y=Protein)) +
     geom_point(alpha=0.5) +
     labs(x= CpG, y=gene)
    plots[[i]] <- p
}

# Save off
pdf(file = paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/plots_neuro_assocs/plots/neuro_pQTMs.pdf"))
for (i in 1:41) {
    print(plots[[i]])
}
dev.off()

