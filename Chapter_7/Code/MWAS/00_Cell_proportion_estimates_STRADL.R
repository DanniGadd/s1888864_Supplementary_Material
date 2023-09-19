###############################################################################################

### ESTIMATES FOR CELL PROPORTIONS USING MEFFIL

###############################################################################################

## We are redoing all cell estimates for both sets using meffil and will use these in
## MWAS analysis to adjust for WBCs, including eosinophils.

# # Install meffil from command line 
# cd /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/meffil/

# wget https://github.com/perishky/meffil/archive/master.zip
#  unzip master.zip
#  mv meffil-master meffil
#  R CMD INSTALL meffil

library(data.table) 
library(tidyverse)

###############################################################################################

# Load original WBC estimates used
w2cells <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/stradl_w2_cells.csv")
w3cells <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/stradl_w3_cells.csv")

# Read in methylation files
w2 = readRDS("/Cluster_Filespace/Marioni_Group/Rob/Somalogic/wave2-STRADL-mvals.rds")
w3 = readRDS("/Cluster_Filespace/Marioni_Group/Rob/Somalogic/wave3-STRADL-mvals.rds")
w2 = w2[which(row.names(w2) %in% row.names(w3)),]
w3 = w3[which(row.names(w3) %in% row.names(w2)),]

# Subset to WBC estimate groups used previously 
w2 <- w2[,which(colnames(w2) %in% w2cells$DNAm_id)]
w3 <- w3[,which(colnames(w3) %in% w3cells$DNAm_id)]

# Remove 3 individuals that do not have depression data available (NAs)
SCID <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/Depression_case_controls_200521/prep_files/SCID_joint_to_GP_hospital_combined_SCID_only_270121.csv")
SCID <- SCID[c(1,2)]
SCID <- na.omit(SCID)
names(SCID)[1] <- "Stradl_id"
target <- read.table("/Cluster_Filespace/Marioni_Group/Rob/Somalogic/STRADL_DNAm_target_REM_17April2020.txt",header =T)
SCID <- left_join(SCID, target, by = "Stradl_id")

w2 <- w2[,which(colnames(w2) %in% SCID$DNAm_id)]
w3 <- w3[,which(colnames(w3) %in% SCID$DNAm_id)]

# Save out IDs for use in EWAS
names1 <- colnames(w2) %>% as.data.frame()
names2 <- colnames(w3) %>% as.data.frame()
names <- rbind(names1, names2)
names(names)[1] <- "DNAm_id"
target <- target[c(1:3)]
names <- left_join(names, target, by = "DNAm_id")
write.csv(names, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/774_ids_MWAS.csv", row.names = F)

# > dim(w2)
# [1] 772619    476
# > dim(w3)
# [1] 772619    298

# Create estimates using the two wave files - run by Daniel

# DLM
# > dim(w2)
# [1] 772619    479
# > dim(w3)
# [1] 772619    299

## Apply meffil estimate 
# Matrix of methylation levels (rows = CpG sites, columns = subjects) - use this as the input 
# Do each set individually 
require(meffil)
est_w2 <- meffil.estimate.cell.counts.from.betas(w2, cell.type.reference="blood gse35069 complete")
est_w3 <- meffil.estimate.cell.counts.from.betas(w3, cell.type.reference="blood gse35069 complete")

# Save out 
# write.csv(est_w2, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/meffil/w2__cell_estimates.csv", row.names = F)
# write.csv(est_w3, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/meffil/w3_299_cell_estimates.csv", row.names = F)

# DLM: I don't have permission to write to the subdirectories so I've written to Danni folder. Please copy across to above directories
write.csv(est_w2, "/Cluster_Filespace/Marioni_Group/Danni/w2_cell_estimates.csv", row.names = F)
write.csv(est_w3, "/Cluster_Filespace/Marioni_Group/Danni/w3_299_cell_estimates.csv", row.names = F)

est_w2_2 <- meffil.estimate.cell.counts.from.betas(w2, cell.type.reference="blood gse35069") saveRDS("/Cluster_Filespace/Marioni_Group/Danni/w2_cell_estimates_6types.rds")
est_w3_2 <- meffil.estimate.cell.counts.from.betas(w3, cell.type.reference="blood gse35069")saveRDS("/Cluster_Filespace/Marioni_Group/Danni/w3_cell_estimates_6types.rds")

# Daniel - has generated both sets of reference panel estimates and saved them out for use 

###############################################################################################

### COMPARE ESTIMATES

# Used previously:
# phenotypes$Bcell + phenotypes$CD4T + phenotypes$CD8T + phenotypes$Gran + phenotypes$NK
# Updated:
# "Bcell" "CD4T"  "CD8T"  "Eos"   "Mono"  "Neu"   "NK"

library(tidyverse)
library(ggcorrplot)

# New estimates (generated through meffil with first reference panel)
w2_6types <- readRDS("/Cluster_Filespace/Marioni_Group/Danni/w2_cell_estimates_6types.rds")
w3_6types <- readRDS("/Cluster_Filespace/Marioni_Group/Danni/w3_cell_estimates_6types.rds")
w2_6types <- as.data.frame(w2_6types)
w3_6types <- as.data.frame(w3_6types)

# Look at whether rows sum to 1
# w2_6types <- w2_6types[,-1]
# w3_6types <- w3_6types[,-1]
# rowSums(w2_6types)
# rowSums(w3_6types) # these dont quite sum to 1 either

# New estimates (generated through meffil with second reference panel) and subset to 774
w2 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/w2_cell_estimates.csv")
w3 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/w3_299_cell_estimates.csv")
w2 <- w2[which(w2$X %in% rownames(w2_6types)),]
w3 <- w3[which(w3$X %in% rownames(w3_6types)),]

# Look at whether rows sum to 1
# w2 <- w2[,-1]
# w3 <- w3[,-1]
# rowSums(w2)
# rowSums(w3) # these dont quite sum to 1

## How do the Bcell, CD4T and CD8T correlate between meffil versions 
identical(w2$X, rownames(w2_6types)) # TRUE
cor.test(w2$Bcell, w2_6types$Bcell) # 0.97
cor.test(w2$CD4T, w2_6types$CD4T) # 0.99
cor.test(w2$CD8T, w2_6types$CD8T) # 0.97

# I will therefore take the set we used previously (Bcell, CD4T, CD8T, Gran, NK, Mono)
# And i will add eosinophils and neutrophils to it from the separate reference panel 
# Check what the correlations look like 

w2 <- w2[c(1,5,7)]
w3 <- w3[c(1,5,7)]

w2_6types$X <- rownames(w2_6types)
w3_6types$X <- rownames(w3_6types)

# Save out and plot correlations for cell counts to use in study 

target <- read.table("/Cluster_Filespace/Marioni_Group/Rob/Somalogic/STRADL_DNAm_target_REM_17April2020.txt",header =T)
target <- target[c(2,3)]
names(target)[2] <- "X"

w2 <- left_join(w2, w2_6types, by = "X")
w2 <- left_join(w2, target, by = "X")
w2 <- w2[c(10,2:9)]
write.csv(w2, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/w2_cells_joint.csv", row.names = F)

w2cor <- w2[-1]
cor <- cor(w2cor)
pdf("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/w2_WBCs_corrplot.pdf")
ggcorrplot(cor, type = "lower",
   lab = TRUE)
dev.off()


w3 <- left_join(w3, w3_6types, by = "X")
w3 <- left_join(w3, target, by = "X")
w3 <- w3[c(10,2:9)]
write.csv(w3, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/w3_cells_joint.csv", row.names = F)


w3cor <- w3[-1]
cor2 <- cor(w3cor)
pdf("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/w3_WBCs_corrplot.pdf")
ggcorrplot(cor2, type = "lower",
   lab = TRUE)
dev.off()



# ## How do the updated meffils correlate with the previous versions used

# # Load original WBC estimates used
# w2cells <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/stradl_w2_cells.csv")
# w3cells <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/stradl_w3_cells.csv")

# w2cells <- w2cells[which(w2cells$DNAm_id %in% w2$X),]
# w3cells <- w3cells[which(w3cells$DNAm_id %in% w3$X),]

# # Match order 

# w2cells <- w2cells[match(w2$X, w2cells$DNAm_id),]

# cor.test(w2cells$CD4T, w2$CD4T) # 0.91
# cor.test(w2cells$CD8T, w2$CD8T) # 0.82


# w3cells <- w3cells[match(w3$X, w3cells$DNAm_id),]

# cor.test(w3cells$CD4T, w3$CD4T) # 0.85
# cor.test(w3cells$CD8T, w3$CD8T) # 0.87



###############################################################################################

# Paste session info below to track versions of packages used
# sessionInfo()

# R version 4.0.3 (2020-10-10)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: openSUSE Leap 15.1

# Matrix products: default
# BLAS:   /usr/local/lib64/R/lib/libRblas.so
# LAPACK: /usr/local/lib64/R/lib/libRlapack.so

# locale:
#  [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C
#  [3] LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8
#  [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8
#  [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C
#  [9] LC_ADDRESS=C               LC_TELEPHONE=C
# [11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C

# attached base packages:
# [1] parallel  stats     graphics  grDevices utils     datasets  methods
# [8] base

# other attached packages:
#  [1] meffil_1.3.1          preprocessCore_1.52.0 SmartSVA_0.1.3
#  [4] RSpectra_0.16-0       isva_1.9              JADE_2.0-3
#  [7] qvalue_2.22.0         gdsfmt_1.26.1         statmod_1.4.35
# [10] quadprog_1.5-8        DNAcopy_1.64.0        fastICA_1.2-3
# [13] lme4_1.1-25           Matrix_1.2-18         multcomp_1.4-14
# [16] TH.data_1.0-10        survival_3.2-7        mvtnorm_1.1-1
# [19] matrixStats_0.57.0    markdown_1.1          gridExtra_2.3
# [22] Cairo_1.5-12.2        knitr_1.30            reshape2_1.4.4
# [25] plyr_1.8.6            ggplot2_3.3.5         sva_3.38.0
# [28] BiocParallel_1.24.0   genefilter_1.72.0     mgcv_1.8-33
# [31] nlme_3.1-150          limma_3.46.0          sandwich_3.0-0
# [34] lmtest_0.9-39         zoo_1.8-8             MASS_7.3-53
# [37] illuminaio_0.32.0

# loaded via a namespace (and not attached):
#  [1] Biobase_2.50.0       httr_1.4.2           edgeR_3.32.1
#  [4] bit64_4.0.5          splines_4.0.3        assertthat_0.2.1
#  [7] askpass_1.1          stats4_4.0.3         blob_1.2.1
# [10] pillar_1.6.1         RSQLite_2.2.1        lattice_0.20-41
# [13] glue_1.4.2           digest_0.6.27        minqa_1.2.4
# [16] colorspace_2.0-1     XML_3.99-0.5         pkgconfig_2.0.3
# [19] purrr_0.3.4          xtable_1.8-4         scales_1.1.1
# [22] tibble_3.1.2         openssl_1.4.3        annotate_1.68.0
# [25] generics_0.1.0       IRanges_2.24.0       ellipsis_0.3.2
# [28] withr_2.4.2          BiocGenerics_0.36.0  magrittr_2.0.1
# [31] crayon_1.4.1         memoise_1.1.0        fansi_0.5.0
# [34] tools_4.0.3          lifecycle_1.0.0      stringr_1.4.0
# [37] S4Vectors_0.28.0     munsell_0.5.0        locfit_1.5-9.4
# [40] cluster_2.1.0        AnnotationDbi_1.52.0 base64_2.0
# [43] compiler_4.0.3       rlang_0.4.11         nloptr_1.2.2.2
# [46] grid_4.0.3           boot_1.3-25          gtable_0.3.0
# [49] codetools_0.2-16     DBI_1.1.0            R6_2.5.0
# [52] dplyr_1.0.7          bit_4.0.4            utf8_1.2.1
# [55] clue_0.3-60          stringi_1.6.2        Rcpp_1.0.6
# [58] vctrs_0.3.8          tidyselect_1.1.1     xfun_0.18

###############################################################################################