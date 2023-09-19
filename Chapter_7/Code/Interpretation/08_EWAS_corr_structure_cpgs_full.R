############################################################################################

### CORRELATION STRUCTURE FOR CPGS IN FULL MODEL AND ASSOC PLOTS FOR PQTMS

############################################################################################

cd /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/corr_structure/

screen

R

library(data.table)

### LOAD DNAm as per fully-adjusted models
## Read in methylation files and transpose so CpGs are columns and IDs are rows
library(data.table) 
w2 = readRDS("/Cluster_Filespace/Marioni_Group/Rob/Somalogic/wave2-STRADL-mvals.rds")
w3 = readRDS("/Cluster_Filespace/Marioni_Group/Rob/Somalogic/wave3-STRADL-mvals.rds")
w2 = w2[which(row.names(w2) %in% row.names(w3)),]
w3 = w3[which(row.names(w3) %in% row.names(w2)),]

## Combine w2 and w3 
meth = cbind(w2,w3)

## Add in GS ids instead of DNAm_ids 
stradl_ids <- read.table("/Cluster_Filespace/Marioni_Group/Rob/Somalogic/STRADL_DNAm_target_REM_17April2020.txt",header =T)
id.778 <- read.csv("/Cluster_Filespace/Marioni_Group/Rob/Somalogic/IDs_778.csv")

## Subset to the people in the ID file 
stradl_ids <- stradl_ids[which(stradl_ids$GS_id %in% id.778$x),]
ids = stradl_ids$DNAm_id
meth <- meth[,match(ids, colnames(meth))]
table(stradl_ids$DNAm_id == colnames(meth))
colnames(meth) <- stradl_ids$GS_id

## Prepare phenotypes file for covariates that we will use to pre-correct data 
stradl_ids <- stradl_ids[,c("GS_id", "age_stradl", "sex", "wave", "Stradl_id")]

# Add batch 
w2_batch <- read.csv("/Cluster_Filespace/Marioni_Group/Rob/Somalogic/wave2_batch.csv") 
w2_batch$Batch[w2_batch$Batch == 1]<- "w1_1"
w2_batch$Batch[w2_batch$Batch == 2]<- "w1_2"
w2_batch$Batch[w2_batch$Batch == 3]<- "w1_3"
w2_batch$Batch[w2_batch$Batch == 4]<- "w1_4"
w2_batch$Batch[w2_batch$Batch == 5]<- "w1_5"
w2_batch$Batch[w2_batch$Batch == 6]<- "w1_6"
w3_batch <- read.csv("/Cluster_Filespace/Marioni_Group/Rob/Somalogic/wave_3_original_samplesheet.csv") 
w2_batch <- w2_batch[,c("Sample_Name","Batch")]
names(w2_batch) <- c("GS_id", "Batch")
w3_batch <- w3_batch[,c("Sample_Name","Batch_all")]
names(w3_batch) <- c("GS_id", "Batch")
w3_batch$GS_id <- gsub("ST", "", w3_batch$GS_id)
batch = rbind(w2_batch, w3_batch)
batch = batch[-which(duplicated(batch$GS_id)),]

phenotypes= merge(stradl_ids, batch, by = "GS_id")

# Join in SCID depression file
dep <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/Depression_case_controls_200521/prep_files/SCID_joint_to_GP_hospital_combined_SCID_only_270121.csv")
dep <- dep[c(1,2)]
names(dep)[2] <- "combined"
dep$combined <- as.character(dep$combined)
names(dep)[1] <- "Stradl_id"
# which(phenotypes$Stradl_id %in% dep$Stradl_id) 
phenotypes = merge(phenotypes, dep, by = "Stradl_id")

# Subset to 774 individuals with complete depression information
sample <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/774_IDs_meth_order.csv")
phenotypes <- phenotypes[which(phenotypes$Stradl_id %in% sample$Stradl_id),]

# Add WBC info into phenotypes
w2cells <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/w2_cells_joint.csv")
w3cells <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/w3_cells_joint.csv")
merged <- rbind(w2cells, w3cells)
names(merged)[1] <- "GS_id"

phenotypes = merge(phenotypes, merged, by = "GS_id")

# Join in epismoker info for those in STRADL to linker file, then to cov file 
epi <- readRDS("/Cluster_Filespace/Marioni_Group/Danni/stradl_epismoker.rds")
names(epi)[1] <- "DNAm_id"
target <- read.table("/Cluster_Filespace/Marioni_Group/Rob/Somalogic/STRADL_DNAm_target_REM_17April2020.txt",header =T)
target <- target[c(2,3)]
names(target)[2] <- "DNAm_id"
library(tidyverse)
epi <- left_join(epi, target, by = "DNAm_id")

phenotypes = merge(phenotypes, epi, by = "GS_id")

# Add BMI
BMI <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS/01_Phenotype_collation/prot_file_270121.csv")
BMI <- BMI[c("GS_id", "BMI")]

phenotypes = merge(phenotypes, BMI, by = "GS_id")

# Subset meth file to 774 individuals with complete information
meth <- meth[,which(colnames(meth) %in% phenotypes$GS_id)]

# Match ID order between phenotype and meth files 
ids= colnames(meth)
phenotypes = phenotypes[match(ids, as.character(phenotypes$GS_id)),]
table(colnames(meth) == phenotypes$GS_id)

# Transpose methylation data
meth = t(meth)

# Resdualise meth matrix 
for(i in 1:(ncol(meth))){
  print(i)
  meth[,i] <- resid(lm(meth[,i] ~ phenotypes$age_stradl + as.factor(phenotypes$sex) + as.factor(phenotypes$Batch) + as.factor(phenotypes$wave) + as.factor(phenotypes$combined) + phenotypes$Bcell + phenotypes$CD4T + phenotypes$CD8T + phenotypes$Gran + phenotypes$NK + phenotypes$Mono + phenotypes$Eos + phenotypes$smokingScore + phenotypes$BMI, na.action = na.exclude))
}

meth <- as.data.frame(meth)

## Combined IDS of individuals (from the relatedness subtyping script for GRM)
osca_dat <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/774_IDs_meth_order.csv")

# match meth file to the order in the FID IID 774 file 
ids = osca_dat$GS_id
meth <- meth[match(ids, row.names(meth)),]
table(osca_dat$GS_id == row.names(meth))

## Combine IDs for FID and IID and rearrange columns so it goes FID, IID, cpg1, cpg2... 
tmp <- meth
tmp$IID <- osca_dat$IID
tmp$FID <- osca_dat$FID
tmp <- tmp[,c(ncol(tmp), ncol(tmp)-1, 1:(ncol(tmp)-2))]
tmp$FID <- as.character(tmp$FID)
tmp$IID <- as.character(tmp$IID)
identical(tmp$IID, as.character(osca_dat$IID)) # true 

##############################################################

full <- tmp

library(tidyverse)
library(readxl)
library(data.table)

# Fully adjusted results
cpgs <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/cis_trans/FULL_cistrans_table_formatted_naming_ANNOTATED.csv")
names(cpgs)[10] <- "gene" # 2928

# Create variables for PAPPA and PRG3
PRG3 <- cpgs[which(cpgs$gene == "PRG3"),] 
dim(PRG3) # 1116

PAPPA <- cpgs[which(cpgs$gene == "PAPPA"),] 
dim(PAPPA) # 987

## Subset full results to remove the 2 highly pleiotropic proteins 
cpgs <- cpgs[-which(cpgs$gene == "PRG3"),] 
cpgs <- cpgs[-which(cpgs$gene == "PAPPA"),] 
dim(cpgs) # 825
length(unique(cpgs$CpG)) # 423 CpGs in 580 associations - updated to 597 CpGs in 825 associations

# Subset DNAm to three groups as above
meth_cpgs <- full[,which(colnames(full) %in% cpgs$CpG)]
meth_PRG3 <- full[,which(colnames(full) %in% PRG3$CpG)]
meth_PAPPA <- full[,which(colnames(full) %in% PAPPA$CpG)]

# > dim(meth_cpgs)
# [1] 774 597
# > dim(meth_PAPPA)
# [1] 774 987
# > dim(meth_PRG3)
# [1]  774 1116

# Extract all as a file too 
# Fully adjusted results
cpgs2 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/cis_trans/FULL_cistrans_table_formatted_naming_ANNOTATED.csv")
names(cpgs2)[10] <- "gene" # 2928
meth_all <- full[,which(colnames(full) %in% cpgs2$CpG)]

# > dim(meth_all)
# [1]  774 1837

# Top 100 P PAPPA
PAPPA_top <- PAPPA[c(1:100),]

# Top 100 P PRG3 
PRG3_top <- PRG3[c(1:100),]
meth_PRG3_top <- full[,which(colnames(full) %in% PRG3_top$CpG)]
meth_PAPPA_top <- full[,which(colnames(full) %in% PAPPA_top$CpG)]


##############################################################################

### PLOT pQTM associations 

slice <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/cis_trans/slice_neuro_final_35_pQTMs_pQTLs_added.csv")
prot <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS/01_Phenotype_collation/prot_file_220221.csv", check.names = F)

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
     geom_point(alpha=0.5, color = "darkblue") +
     labs(x= CpG, y=gene) + theme_classic() + 
  theme(legend.title = element_blank(),
    legend.position = "none", 
    axis.title.x = element_text(size = 27), 
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.y = element_text(size = 27),
    plot.title = element_text(size = 27),
    axis.title.x.bottom = element_text(margin = margin(14, 0, 0, 0)))
    plots[[i]] <- p
}

# library(gridExtra)
ggsave("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/corr_structure_MWAS/arrange.pdf", arrangeGrob(grobs = plots), device = "pdf", height = 20, width = 25)



#############################################################

# Run a PCA on the signals in each subgroup 

## CPGS

library(factoextra)
res.pca <- prcomp(meth_cpgs, scale = TRUE)

# names(res.pca)
# # [1] "sdev"     "rotation" "center"   "scale"    "x"
loadings <- res.pca$rotation # this provides the loadings

loadings[1:5,1:4]

#                 PC1          PC2           PC3          PC4
# 10000-28 0.01578677  0.006109059  0.0086136447 -0.010829816
# 10001-7  0.01712065 -0.020603307 -0.0022774886  0.010428056
# 10003-15 0.01163660  0.014874150  0.0088900875  0.012455513
# 10006-25 0.01947500 -0.008470357 -0.0055983740 -0.003819188
# 10008-43 0.01578426  0.003989797 -0.0008580165 -0.017475235

# the matrix x has the principal component score vectors
dim(res.pca$x)
# [1] 1065 1065

# Compute variance explained by componenets using standard deviations
std_dev <- res.pca$sdev
pr_var <- std_dev^2
prop_varex <- pr_var/sum(pr_var)
head(prop_varex)
# [1] 0.469395583 0.066381867 0.026601245 0.018127898 0.014454832 0.009869478

# Plot cumulative variance explained
pdf("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/corr_structure_MWAS/PCA_cpgs.pdf")
plot(cumsum(prop_varex), xlab = "Principal Component",
              ylab = "Cumulative Proportion of Variance Explained",
              type = "b")
dev.off()

# Eigenvalues
eig.val <- get_eigenvalue(res.pca)
# eig.val
write.csv(eig.val, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/corr_structure_MWAS/cpgs_eig_values_prcomp.csv")

# Create corrplot 
library(foreign) 
library(dplyr)
library(psych)

# ggcorrplot
library(ggcorrplot)
corr <- cor(meth_cpgs)
pdf("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/corr_structure_MWAS/ggcorrplot_cpgs.pdf", width = 50, height = 50)
ggcorrplot(corr, hc.order = TRUE, type = "lower",
     outline.col = "white") + theme(text = element_text(size = 0.2)) 
dev.off()

# Create plots for cum var and eigen values 
var <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/corr_structure_MWAS/cpgs_eig_values_prcomp.csv")
var$num <- 1:597
names(var)[2] <- "mes"
# var[1,1] <- 300
var$col <- ifelse(var$mes >= 1, "darkgrey", "orange")
dim(var[which(var$mes >= 1),]) # 150 eigenvalues greater than or equal to 1 
# Eigen values 
pdf("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/corr_structure_MWAS/updated_130322_extra_plots/CPGS_597_eigen_vals.pdf")
ggplot(var, aes(num, mes, col = col)) +
geom_point(size = 1) + theme(legend.position = "none", axis.text=element_text(size=13), axis.title=element_text(size=13)) +
xlab("Principal component") + ylab("Eigenvalue") + geom_hline(yintercept=1, color = "darkgrey", size=1) + 
scale_color_manual(values=c("#E69F00", "#999999")) 
dev.off()

# first eigenvalue was 1987.89029 - but set to 300 for purposes of visualisation
names(var)[4] <- "cum"
# Cumulative variance 
pdf("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/corr_structure_MWAS/updated_130322_extra_plots/CPGS_597_cum_var.pdf")
ggplot(var, aes(num, cum)) +
geom_bar(stat = "identity", fill = "steelblue") + theme(legend.position = "none", axis.text=element_text(size=13), axis.title=element_text(size=13)) +
xlab("Principal component") + ylab("Cumulative proportion") 
dev.off()

########################################################

## PAPPA

library(factoextra)
res.pca <- prcomp(meth_PAPPA, scale = TRUE)

# names(res.pca)
# # [1] "sdev"     "rotation" "center"   "scale"    "x"
loadings <- res.pca$rotation # this provides the loadings

# the matrix x has the principal component score vectors
dim(res.pca$x)
# [1] 1065 1065

# Compute variance explained by componenets using standard deviations
std_dev <- res.pca$sdev
pr_var <- std_dev^2
prop_varex <- pr_var/sum(pr_var)
head(prop_varex)
# [1] 0.469395583 0.066381867 0.026601245 0.018127898 0.014454832 0.009869478

# Plot cumulative variance explained
pdf("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/corr_structure_MWAS/PCA_PAPPA.pdf")
plot(cumsum(prop_varex), xlab = "Principal Component",
              ylab = "Cumulative Proportion of Variance Explained",
              type = "b")
dev.off()

# Eigenvalues
eig.val <- get_eigenvalue(res.pca)
eig.val
write.csv(eig.val, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/corr_structure_MWAS/PAPPA_eig_values_prcomp.csv")


# Create plots for cum var and eigen values 
var <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/corr_structure_MWAS/PAPPA_eig_values_prcomp.csv")
var$num <- 1:774
names(var)[2] <- "mes"
# var[1,1] <- 300
var$col <- ifelse(var$mes >= 1, "darkgrey", "orange")
dim(var[which(var$mes >= 1),]) # 150 eigenvalues greater than or equal to 1 
# Eigen values 
pdf("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/corr_structure_MWAS/updated_130322_extra_plots/PAPPA_eigen_vals.pdf")
ggplot(var, aes(num, mes, col = col)) +
geom_point(size = 1) + theme(legend.position = "none", axis.text=element_text(size=13), axis.title=element_text(size=13)) +
xlab("Principal component") + ylab("Eigenvalue") + geom_hline(yintercept=1, color = "darkgrey", size=1) + 
scale_color_manual(values=c("#E69F00", "#999999")) 
dev.off()

# first eigenvalue was 1987.89029 - but set to 300 for purposes of visualisation

names(var)[4] <- "cum"
# Cumulative variance 
pdf("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/corr_structure_MWAS/updated_130322_extra_plots/PAPPA_cum_var.pdf")
ggplot(var, aes(num, cum)) +
geom_bar(stat = "identity", fill = "steelblue") + theme(legend.position = "none", axis.text=element_text(size=13), axis.title=element_text(size=13)) +
xlab("Principal component") + ylab("Cumulative proportion") 
dev.off()

########################################################

## PRG3

library(factoextra)
res.pca <- prcomp(meth_PRG3, scale = TRUE)

# names(res.pca)
# # [1] "sdev"     "rotation" "center"   "scale"    "x"
loadings <- res.pca$rotation # this provides the loadings

# the matrix x has the principal component score vectors
dim(res.pca$x)
# [1] 1065 1065

# Compute variance explained by componenets using standard deviations
std_dev <- res.pca$sdev
pr_var <- std_dev^2
prop_varex <- pr_var/sum(pr_var)
head(prop_varex)
# [1] 0.469395583 0.066381867 0.026601245 0.018127898 0.014454832 0.009869478

# Plot cumulative variance explained
pdf("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/corr_structure_MWAS/PCA_PRG3.pdf")
plot(cumsum(prop_varex), xlab = "Principal Component",
              ylab = "Cumulative Proportion of Variance Explained",
              type = "b")
dev.off()

# Eigenvalues
eig.val <- get_eigenvalue(res.pca)
eig.val
write.csv(eig.val, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/corr_structure_MWAS/PRG3_eig_values_prcomp.csv")


# Create plots for cum var and eigen values 
var <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/corr_structure_MWAS/PRG3_eig_values_prcomp.csv")
var$num <- 1:774
names(var)[2] <- "mes"
# var[1,1] <- 300
var$col <- ifelse(var$mes >= 1, "darkgrey", "orange")
dim(var[which(var$mes >= 1),]) # 150 eigenvalues greater than or equal to 1 
# Eigen values 
pdf("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/corr_structure_MWAS/updated_130322_extra_plots/PRG3_eigen_vals.pdf")
ggplot(var, aes(num, mes, col = col)) +
geom_point(size = 1) + theme(legend.position = "none", axis.text=element_text(size=13), axis.title=element_text(size=13)) +
xlab("Principal component") + ylab("Eigenvalue") + geom_hline(yintercept=1, color = "darkgrey", size=1) + 
scale_color_manual(values=c("#E69F00", "#999999")) 
dev.off()

# first eigenvalue was 1987.89029 - but set to 300 for purposes of visualisation

names(var)[4] <- "cum"
# Cumulative variance 
pdf("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/corr_structure_MWAS/updated_130322_extra_plots/PRG3_cum_var.pdf")
ggplot(var, aes(num, cum)) +
geom_bar(stat = "identity", fill = "steelblue") + theme(legend.position = "none", axis.text=element_text(size=13), axis.title=element_text(size=13)) +
xlab("Principal component") + ylab("Cumulative proportion") 
dev.off()


########################################################

## ALL

library(factoextra)
res.pca <- prcomp(meth_all, scale = TRUE)

# names(res.pca)
# # [1] "sdev"     "rotation" "center"   "scale"    "x"
loadings <- res.pca$rotation # this provides the loadings

# the matrix x has the principal component score vectors
dim(res.pca$x)
# [1] 1065 1065

# Compute variance explained by componenets using standard deviations
std_dev <- res.pca$sdev
pr_var <- std_dev^2
prop_varex <- pr_var/sum(pr_var)
head(prop_varex)
# [1] 0.469395583 0.066381867 0.026601245 0.018127898 0.014454832 0.009869478

# Plot cumulative variance explained
pdf("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/corr_structure_MWAS/PCA_all_1837.pdf")
plot(cumsum(prop_varex), xlab = "Principal Component",
              ylab = "Cumulative Proportion of Variance Explained",
              type = "b")
dev.off()

# Eigenvalues
eig.val <- get_eigenvalue(res.pca)
eig.val
write.csv(eig.val, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/corr_structure_MWAS/ALL_eig_values_prcomp.csv")


# Create plots for cum var and eigen values 
var <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/corr_structure_MWAS/ALL_eig_values_prcomp.csv")
var$num <- 1:774
names(var)[2] <- "mes"
# var[1,1] <- 300
var$col <- ifelse(var$mes >= 1, "darkgrey", "orange")
dim(var[which(var$mes >= 1),]) # 150 eigenvalues greater than or equal to 1 
# Eigen values 
pdf("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/corr_structure_MWAS/updated_130322_extra_plots/ALL_eigen_vals.pdf")
ggplot(var, aes(num, mes, col = col)) +
geom_point(size = 1) + theme(legend.position = "none", axis.text=element_text(size=13), axis.title=element_text(size=13)) +
xlab("Principal component") + ylab("Eigenvalue") + geom_hline(yintercept=1, color = "darkgrey", size=1) + 
scale_color_manual(values=c("#E69F00", "#999999")) 
dev.off()

# first eigenvalue was 1987.89029 - but set to 300 for purposes of visualisation

names(var)[4] <- "cum"
# Cumulative variance 
pdf("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/corr_structure_MWAS/updated_130322_extra_plots/ALL_cum_var.pdf")
ggplot(var, aes(num, cum)) +
geom_bar(stat = "identity", fill = "steelblue") + theme(legend.position = "none", axis.text=element_text(size=13), axis.title=element_text(size=13)) +
xlab("Principal component") + ylab("Cumulative proportion") 
dev.off()


# ########################################################

# ### TAKE TOP 100 FROM PAPPA AND PRG3 and correlate with main cpgs 

# # # ggcorrplot
# # library(ggcorrplot)
# # corr <- cor(meth_PAPPA_top)
# # pdf("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/corr_structure_MWAS/ggcorrplot_top_100_PAPPA.pdf", width = 60, height = 60)
# # ggcorrplot(corr,
# #      outline.col = "white") + theme(text = element_text(size = 0.2)) 
# # dev.off()


# # # ggcorrplot
# # library(ggcorrplot)
# # corr <- cor(meth_PRG3_top)
# # pdf("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/corr_structure_MWAS/ggcorrplot_top_100_PRG3.pdf", width = 60, height = 60)
# # ggcorrplot(corr,
# #      outline.col = "white") + theme(text = element_text(size = 0.2)) 
# # dev.off()


# # ggcorrplot
# library(ggcorrplot)
# corr <- cor(meth_cpgs)
# pdf("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/corr_structure_MWAS/ggcorrplot_597_unique_cpgs.pdf", width = 25, height = 25)
# ggcorrplot(corr, type = "upper",
#      outline.col = "white", hc.order = TRUE) +
#   theme(axis.text.x=element_text(size=1),
#         axis.text.y=element_text(size=1), legend.position = "none")
# dev.off()

# # ggcorrplot
# library(ggcorrplot)
# corr <- cor(meth_PAPPA)
# pdf("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/corr_structure_MWAS/ggcorrplot_PAPPA_all.pdf", width = 25, height = 25)
# ggcorrplot(corr, type = "upper",
#      outline.col = "white", hc.order = TRUE) +
#   theme(axis.text.x=element_text(size=1),
#         axis.text.y=element_text(size=1), legend.position = "none")
# dev.off()

# # ggcorrplot
# library(ggcorrplot)
# corr <- cor(meth_PRG3)
# pdf("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/corr_structure_MWAS/ggcorrplot_PRG3_all.pdf", width = 25, height = 25)
# ggcorrplot(corr, type = "upper",
#      outline.col = "white", hc.order = TRUE) +
#   theme(axis.text.x=element_text(size=1),
#         axis.text.y=element_text(size=1), legend.position = "none")
# dev.off()

# # ggcorrplot
# library(ggcorrplot)
# corr <- cor(meth_all)
# pdf("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/corr_structure_MWAS/ggcorrplot_ALL_cpgs.pdf", width = 25, height = 25)
# ggcorrplot(corr, type = "upper",
#      outline.col = "white", hc.order = TRUE) +
#   theme(axis.text.x=element_text(size=1),
#         axis.text.y=element_text(size=1), legend.position = "none")
# dev.off()
