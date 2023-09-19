##############################################################################################################

### PROCESS PheWAS models

##############################################################################################################

## Assess the independent signals in the 4235 somamers for inclusion in PheWAS/MWAS
## This will be used to decide the threshold for multuple testing adjustment in PheWAS vs FDR levels

# prcomp: https://www.analyticsvidhya.com/blog/2016/03/pca-practical-guide-principal-component-analysis-python/
# prcomp: http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/118-principal-component-analysis-in-r-prcomp-vs-princomp/

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
prot <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS/01_Phenotype_collation/prot_file_270121.csv", check.names = F)

# Isolate just the protein columns of interest for PCA
prot <- prot[c(33:4267)]

# Try prcomp
library(factoextra)
res.pca <- prcomp(prot, scale = TRUE)

names(res.pca)
# [1] "sdev"     "rotation" "center"   "scale"    "x"

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
plot(cumsum(prop_varex), xlab = "Principal Component",
              ylab = "Cumulative Proportion of Variance Explained",
              type = "b")

# Eigenvalues
eig.val <- get_eigenvalue(res.pca)
eig.val
write.csv(eig.val, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS/00_PCA_proteins/eig_values_prcomp.csv")


# Results for Variables
res.var <- get_pca_var(res.pca)
res.var$coord          # Coordinates
res.var$contrib        # Contributions to the PCs
res.var$cos2           # Quality of representation 
# Results for individuals
res.ind <- get_pca_ind(res.pca)
res.ind$coord          # Coordinates
res.ind$contrib        # Contributions to the PCs
res.ind$cos2           # Quality of representation 

# The function princomp() uses the spectral decomposition approach. 
# The functions prcomp() and PCA()[FactoMineR] use the singular value decomposition (SVD).

# Spectral decomposition which examines the covariances / correlations between variables
# Singular value decomposition which examines the covariances / correlations between individuals

### PLOT EIGENVALUES AND CUMULATIVE VARIANCE WITH CORR PLOT 

var <- eig.val
var$num <- 1:1065
names(var)[1] <- "mes"
var[1,1] <- 300
var$col <- ifelse(var$mes >= 1, "darkgrey", "orange")
dim(var[which(var$mes >= 1),]) # 483 eigenvalues greater than or equal to 1 
# Eigen values 
a <- ggplot(var, aes(num, mes, col = col)) +
geom_point(size = 1) + theme(legend.position = "none", axis.text=element_text(size=13), axis.title=element_text(size=13)) +
xlab("Principal component") + ylab("Eigenvalue") + geom_hline(yintercept=1, color = "darkgrey", size=1) + 
scale_color_manual(values=c("#E69F00", "#999999")) 

pdf("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/corr_structure/4235_proteins_eigen_vals.pdf")
ggplot(var, aes(num, mes, col = col)) +
geom_point(size = 1) + theme(legend.position = "none", axis.text=element_text(size=13), axis.title=element_text(size=13)) +
xlab("Principal component") + ylab("Eigenvalue") + geom_hline(yintercept=1, color = "darkgrey", size=1) + 
scale_color_manual(values=c("#E69F00", "#999999")) 
dev.off()

# first eigenvalue was 1987.89029 - but set to 300 for purposes of visualisation

names(var)[3] <- "cum"
# Cumulative variance 
b <- ggplot(var, aes(num, cum)) +
geom_bar(stat = "identity", fill = "steelblue") + theme(legend.position = "none", axis.text=element_text(size=13), axis.title=element_text(size=13)) +
xlab("Principal component") + ylab("Cumulative proportion") 

pdf("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/corr_structure/4235_proteins_cum_var.pdf")
ggplot(var, aes(num, cum)) +
geom_bar(stat = "identity", fill = "steelblue") + theme(legend.position = "none", axis.text=element_text(size=13), axis.title=element_text(size=13)) +
xlab("Principal component") + ylab("Cumulative proportion") 
dev.off()

# ggcorrplot
library(ggcorrplot)
corr <- cor(prot)
c <- ggcorrplot(corr, hc.order = TRUE, type = "lower",
     outline.col = "white") + theme(text = element_text(size = 0.2)) 

pdf("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/corr_structure/4235_proteins_corr_plot.pdf", width = 50, height = 50)
ggcorrplot(corr, hc.order = TRUE, type = "lower",
     outline.col = "white") + theme(text = element_text(size = 0.2)) 
dev.off()


###################################################
