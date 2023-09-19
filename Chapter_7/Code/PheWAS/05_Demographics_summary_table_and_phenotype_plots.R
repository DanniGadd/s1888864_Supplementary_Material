###############################################################################################

### PLOT STRADL PHEWAS PHENOTYPES 
### CREATE THE DEMOGRAPHICS TABLE SUMMARY FOR THE STRADL COHORT 

# Use n=1095 and n=778 group IDs to create 2 sets of output summaries 

# REVISION UPDATE: Revisit the supplmentary table 1 information to update as requested by reviewers
# This includes further additions (depression, estimated cell proportions, BMI categories, smoking status categories)
# It also includes further clarifications on IQR for scores (cog and brain imaging), and distribution assessments
# Updated summary information will be collated for both analysis groups 
# Also, create distribution plots for cognitive scores and imaging variables (histograms)

###############################################################################################

### PLOT PHENOTYPES

screen

R

library(tidyverse)
library(ggplot2)

## Cognitive scores 

# Read in the processed cognitive data from the script daniel shared with me - composite gf and g scores and other scores with outliers > 3.5 sd from mean removed 
comp <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS/01_Phenotype_collation/prot_file_220221.csv", check.names = F)

# Make plots for cognitive scores 

# Digit symbol 
d <- comp[,c("GS_id", "digit_symbol")]
d <- na.omit(d)
d$digit_symbol <- as.numeric(d$digit_symbol)
mean_value <- mean(d$digit_symbol)
bw <- 2 * IQR(d$digit_symbol) / length(d$digit_symbol)^(1/3)

plot1 <- ggplot(d, aes(x=digit_symbol)) + 
  geom_histogram(color="black", fill="azure3", binwidth = bw) + geom_vline(xintercept = mean_value,      
             col = "dodgerblue4",
             lwd = 2) + xlab("Processing Speed") + ylab("Count") +
 theme_classic() + theme(
    axis.text = element_text(size = 20), 
    axis.title = element_text(size = 20))

pdf("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS/01_Phenotype_collation/Plots/digit.pdf")
plot1
dev.off()

# Logical memory 
d <- comp[,c("GS_id", "LM")]
d <- na.omit(d)
d$LM <- as.numeric(d$LM)
mean_value <- mean(d$LM)
bw <- 2 * IQR(d$LM) / length(d$LM)^(1/3)

plot2 <- ggplot(d, aes(x=LM)) + 
  geom_histogram(color="black", fill="azure3", binwidth = bw) + geom_vline(xintercept = mean_value,      
             col = "dodgerblue4",
             lwd = 2) + xlab("Logical Memory") + ylab("Count") +
 theme_classic() + theme(
    axis.text = element_text(size = 20), 
    axis.title = element_text(size = 20))

pdf("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS/01_Phenotype_collation/Plots/LM.pdf")
plot2
dev.off()

# Verbal total (verbal reasoning)
d <- comp[,c("GS_id", "verbal_total")]
d <- na.omit(d)
d$verbal_total <- as.numeric(d$verbal_total)
mean_value <- mean(d$verbal_total)
bw <- 2 * IQR(d$verbal_total) / length(d$verbal_total)^(1/3)

plot3 <- ggplot(d, aes(x=verbal_total)) + 
  geom_histogram(color="black", fill="azure3", binwidth = bw) + geom_vline(xintercept = mean_value,      
             col = "dodgerblue4",
             lwd = 2) + xlab("Verbal Reasoning") + ylab("Count") +
 theme_classic() + theme(
    axis.text = element_text(size = 20), 
    axis.title = element_text(size = 20))

pdf("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS/01_Phenotype_collation/Plots/verbal_reasoning.pdf")
plot3
dev.off()


# Matrix reasoning (Non-Verbal reasoning)
d <- comp[,c("GS_id", "mr_correct")]
d <- na.omit(d)
d$verbal_total <- as.numeric(d$mr_correct)
mean_value <- mean(d$mr_correct)
bw <- 2 * IQR(d$mr_correct) / length(d$mr_correct)^(1/3)

plot4 <- ggplot(d, aes(x=mr_correct)) + 
  geom_histogram(color="black", fill="azure3", binwidth = bw) + geom_vline(xintercept = mean_value,      
             col = "dodgerblue4",
             lwd = 2) + xlab("Non-Verbal Reasoning") + ylab("Count") +
 theme_classic() + theme(
    axis.text = element_text(size = 20), 
    axis.title = element_text(size = 20))

pdf("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS/01_Phenotype_collation/Plots/nonverbal_reasoning.pdf")
plot4
dev.off()

# Vocbulary
d <- comp[,c("GS_id", "vocabulary")]
d <- na.omit(d)
d$verbal_total <- as.numeric(d$vocabulary)
mean_value <- mean(d$vocabulary)
bw <- 2 * IQR(d$vocabulary) / length(d$vocabulary)^(1/3)

plot5 <- ggplot(d, aes(x=vocabulary)) + 
  geom_histogram(color="black", fill="azure3", binwidth = bw) + geom_vline(xintercept = mean_value,      
             col = "dodgerblue4",
             lwd = 2) + xlab("Vocabulary") + ylab("Count") +
 theme_classic() + theme(
    axis.text = element_text(size = 20), 
    axis.title = element_text(size = 20))

pdf("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS/01_Phenotype_collation/Plots/vocabulary.pdf")
plot5
dev.off()

# g
d <- comp[,c("GS_id", "g")]
d <- na.omit(d)
d$g <- as.numeric(d$g)
mean_value <- mean(d$g)
bw <- 2 * IQR(d$g) / length(d$g)^(1/3)

plot6 <- ggplot(d, aes(x=g)) + 
  geom_histogram(color="black", fill="azure3", binwidth = bw) + geom_vline(xintercept = mean_value,      
             col = "dodgerblue4",
             lwd = 2) + xlab("General Cognitive Ability") + ylab("Count") +
 theme_classic() + theme(
    axis.text = element_text(size = 20), 
    axis.title = element_text(size = 20))

pdf("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS/01_Phenotype_collation/Plots/g.pdf")
plot6
dev.off()


# gf 
d <- comp[,c("GS_id", "gf")]
d <- na.omit(d)
d$gf <- as.numeric(d$gf)
mean_value <- mean(d$gf)
bw <- 2 * IQR(d$gf) / length(d$gf)^(1/3)

plot7 <- ggplot(d, aes(x=gf)) + 
  geom_histogram(color="black", fill="azure3", binwidth = bw) + geom_vline(xintercept = mean_value,      
             col = "dodgerblue4",
             lwd = 2) + xlab("General Fluid Cognitive Ability") + ylab("Count") +
 theme_classic() + theme(
    axis.text = element_text(size = 20), 
    axis.title = element_text(size = 20))

pdf("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS/01_Phenotype_collation/Plots/gf.pdf")
plot7
dev.off()

# join together as suppl plot 

library(patchwork)

pdf("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS/01_Phenotype_collation/Plots/joint_cognitive_V2.pdf", width = 15, height = 10)
plot1 + plot2 + plot3 + plot4 + plot5 + plot6 + plot7
dev.off()


## Imaging 

# Global_GM_Volume
d <- comp[,c("GS_id", "Global_GM_Volume")]
d <- na.omit(d)
d$Global_GM_Volume <- as.numeric(d$Global_GM_Volume)
mean_value <- mean(d$Global_GM_Volume)
bw <- 2 * IQR(d$Global_GM_Volume) / length(d$Global_GM_Volume)^(1/3)

plot1 <- ggplot(d, aes(x=Global_GM_Volume)) + 
  geom_histogram(color="black", fill="azure3", binwidth = bw) + geom_vline(xintercept = mean_value,      
             col = "dodgerblue4",
             lwd = 2) + xlab("Global Grey Matter Volume") + ylab("Count") +
 theme_classic() + theme(
    axis.text = element_text(size = 20), 
    axis.title = element_text(size = 20))

pdf("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS/01_Phenotype_collation/Plots/global_gm_volume.pdf")
plot1
dev.off()


# WMHV
d <- comp[,c("GS_id", "WMH_Volume_Total")]
d <- na.omit(d)
d$WMH_Volume_Total <- as.numeric(d$WMH_Volume_Total)
mean_value <- mean(d$WMH_Volume_Total)
bw <- 2 * IQR(d$WMH_Volume_Total) / length(d$WMH_Volume_Total)^(1/3)

plot2 <- ggplot(d, aes(x=WMH_Volume_Total)) + 
  geom_histogram(color="black", fill="azure3", binwidth = bw) + geom_vline(xintercept = mean_value,      
             col = "dodgerblue4",
             lwd = 2) + xlab("White Matter Hyperintensity Volume") + ylab("Count") +
 theme_classic() + theme(
    axis.text = element_text(size = 20), 
    axis.title = element_text(size = 20))

pdf("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS/01_Phenotype_collation/Plots/WMHV.pdf")
plot2
dev.off()

# WBV
d <- comp[,c("GS_id", "WBV_No_Ventricles")]
d <- na.omit(d)
d$WBV_No_Ventricles <- as.numeric(d$WBV_No_Ventricles)
mean_value <- mean(d$WBV_No_Ventricles)
bw <- 2 * IQR(d$WBV_No_Ventricles) / length(d$WBV_No_Ventricles)^(1/3)

plot3 <- ggplot(d, aes(x=WBV_No_Ventricles)) + 
  geom_histogram(color="black", fill="azure3", binwidth = bw) + geom_vline(xintercept = mean_value,      
             col = "dodgerblue4",
             lwd = 2) + xlab("Whole Brain Volume") + ylab("Count") +
 theme_classic() + theme(
    axis.text = element_text(size = 20), 
    axis.title = element_text(size = 20))

pdf("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS/01_Phenotype_collation/Plots/WBV_no_vent.pdf")
plot3
dev.off()

# GFA
d <- comp[,c("GS_id", "gFA")]
d <- na.omit(d)
d$gFA <- as.numeric(d$gFA)
mean_value <- mean(d$gFA)
bw <- 2 * IQR(d$gFA) / length(d$gFA)^(1/3)

plot4 <- ggplot(d, aes(x=gFA)) + 
  geom_histogram(color="black", fill="azure3", binwidth = bw) + geom_vline(xintercept = mean_value,      
             col = "dodgerblue4",
             lwd = 2) + xlab("General Fractional Anisotropy") + ylab("Count") +
 theme_classic() + theme(
    axis.text = element_text(size = 20), 
    axis.title = element_text(size = 20))

pdf("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS/01_Phenotype_collation/Plots/gFA.pdf")
plot4
dev.off()


# GMD
d <- comp[,c("GS_id", "gMD")]
d <- na.omit(d)
d$gMD <- as.numeric(d$gMD)
mean_value <- mean(d$gMD)
bw <- 2 * IQR(d$gMD) / length(d$gMD)^(1/3)

plot5 <- ggplot(d, aes(x=gMD)) + 
  geom_histogram(color="black", fill="azure3", binwidth = bw) + geom_vline(xintercept = mean_value,      
             col = "dodgerblue4",
             lwd = 2) + xlab("General Mean Diffusivity") + ylab("Count") +
 theme_classic() + theme(
    axis.text = element_text(size = 20), 
    axis.title = element_text(size = 20))

pdf("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS/01_Phenotype_collation/Plots/gMD.pdf")
plot5
dev.off()

# Relative brain age 
d <- comp[,c("GS_id", "brain_accel")]
d <- na.omit(d)
d$brain_accel <- as.numeric(d$brain_accel)
mean_value <- mean(d$brain_accel)
bw <- 2 * IQR(d$brain_accel) / length(d$brain_accel)^(1/3)

plot6 <- ggplot(d, aes(x=brain_accel)) + 
  geom_histogram(color="black", fill="azure3", binwidth = bw) + geom_vline(xintercept = mean_value,      
             col = "dodgerblue4",
             lwd = 2) + xlab("Relative Brain Age") + ylab("Count") +
 theme_classic() + theme(
    axis.text = element_text(size = 20), 
    axis.title = element_text(size = 20))

pdf("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS/01_Phenotype_collation/Plots/Relative_brain_age.pdf")
plot6
dev.off()


# Fazekas WMH
d <- comp[,c("GS_id", "Fazekas_Score_Total")]
d <- na.omit(d)
d$Fazekas_Score_Total <- as.numeric(d$Fazekas_Score_Total)
mean_value <- mean(d$Fazekas_Score_Total)
bw <- 2 * IQR(d$Fazekas_Score_Total) / length(d$Fazekas_Score_Total)^(1/3)

plot7 <- ggplot(d, aes(x=Fazekas_Score_Total)) + 
  geom_histogram(color="black", fill="azure3", binwidth = bw) + geom_vline(xintercept = mean_value,      
             col = "dodgerblue4",
             lwd = 2) + xlab("Fazekas White Matter Hyperintensity Score") + ylab("Count") +
 theme_classic() + theme(
    axis.text = element_text(size = 20), 
    axis.title = element_text(size = 20))

pdf("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS/01_Phenotype_collation/Plots/Faz_WMH.pdf")
plot7
dev.off()

# join together as suppl plot 

library(patchwork)

pdf("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS/01_Phenotype_collation/Plots/joint_imaging_V2.pdf", width = 16, height = 10)
plot1 + plot2 + plot3 + plot4 + plot5 + plot6 + plot7
dev.off()

###############################################################################################

### SUMMARISE PHENOTYPES 

screen

R

library(tidyverse)
library(ggplot2)

# Data from PheWAS used for the plots above 
comp <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS/01_Phenotype_collation/prot_file_220221.csv", check.names = F)

# Load Meth IDs for EWAS to allow for creation of two subsets (1065 and 778)
order <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/matching_data_788.csv")

# order <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/MWAS/774_ids_MWAS.csv")
# target <- read.table("/Cluster_Filespace/Marioni_Group/Rob/Somalogic/STRADL_DNAm_target_REM_17April2020.txt",header =T)
# target <- target[c(2,3)]
# names(target)[2] <- "DNAm_id"

# Calculate 2 subsets (778 vs 1065)
MWAS <- comp[which(comp$GS_id %in% order$GS_id),]

## AGE

# N available: 778 MWAS, 1065 comp
table(is.na(comp$st_age)) # FALSE
mean(comp$st_age, na.rm = T)
sd(comp$st_age, na.rm = T)

# > mean(comp$st_age, na.rm = T)
# [1] 59.92958
# > sd(comp$st_age, na.rm = T)
# [1] 9.561029

table(is.na(MWAS$st_age)) # FALSE
mean(MWAS$st_age, na.rm = T)
sd(MWAS$st_age, na.rm = T)

# > mean(MWAS$st_age, na.rm = T)
# [1] 60.16967
# > sd(MWAS$st_age, na.rm = T)
# [1] 8.817731


## SEX
table(is.na(comp$sex)) # FALSE
table(comp$sex)
#   F   M
# 630 435

table(is.na(MWAS$sex)) # FALSE
table(MWAS$sex)
#   F   M
# 439 339

## COGNITIVE 

# [4283] "logical_mem_1"                 "logical_mem_2"
# [4285] "digit_symbol"                  "verbal_C"
# [4287] "verbal_F"                      "verbal_L"
# [4289] "verbal_total"                  "vocabulary"
# [4291] "mr_correct"                    "mr_errors"
# [4293] "mr_time"                       "LM"
# [4295] "g"                             "gf"

# List phenotypes 
num <- c("digit_symbol", "LM",  "mr_correct", "vocabulary",  "verbal_total", "g", "gf")

# Create a numeric phenotype table  
output <- matrix(nrow = 1*length(num), ncol = 15)
output <- as.data.frame(output)
names(output) <- c("Phenotype", "N", "P", "Mean_1065", "(SD)_1065", "IQR", "range_lower", "range_upper", "N778", "P778", "Mean_778", "(SD)_778", "IQR778", "range778_lower", "range778_upper")

for (i in 1:length(num)){
  # for the 1065
  trait <- as.character(num[i])
  dataset <- comp[,which(colnames(comp) %in% c("GS_id", trait))]
  names(dataset)[2] <- "V2"
  dataset$V2 <- as.numeric(dataset$V2)
  dataset <- na.omit(dataset)
  mean <- round(mean(dataset$V2), digits = 1)
  sd <- round(sd(dataset$V2), digits = 1) 
  N <- dim(dataset)[1]
  P <- ((dim(dataset)[1]) / 1065 )* 100
  IQR <- IQR(dataset$V2)[1]
  range_lower <- range(dataset$V2)[1]
  range_upper <- range(dataset$V2)[2]

  # for the 778 
  dataset2 <- MWAS[,which(colnames(MWAS) %in% c("GS_id", trait))]
  names(dataset2)[2] <- "V2"
  dataset2$V2 <- as.numeric(dataset2$V2)
  dataset2 <- na.omit(dataset2)
  N2 <- dim(dataset2)[1]
  P2 <- ((dim(dataset2)[1]) / 778 )* 100
  mean2 <- round(mean(dataset2$V2), digits = 1)
  sd2 <- round(sd(dataset2$V2), digits = 1) 
  IQR2 <- IQR(dataset2$V2)[1]
  range_lower2 <- range(dataset2$V2)[1]
  range_upper2 <- range(dataset2$V2)[2]

  # round 
  P <- round(P, digits = 1)
  P2 <- round(P2, digits = 1)
  IQR <- round(IQR, digits = 1)
  IQR2 <- round(IQR2, digits = 1)
  range_lower <- round(range_lower, digits = 1)
  range_upper <- round(range_upper, digits = 1)
  range_lower2 <- round(range_lower2, digits = 1)
  range_upper2 <- round(range_upper2, digits = 1)

  # output
  output[i,1] <- trait
  output[i,2] <- N
  output[i,3] <- P
  output[i,4] <- mean
  output[i,5] <- sd
  output[i,6] <- IQR
  output[i,7] <- range_lower
  output[i,8] <- range_upper
  output[i,9] <- N2
  output[i,10] <- P2
  output[i,11] <- mean2
  output[i,12] <- sd2
  output[i,13] <- IQR2
  output[i,14] <- range_lower2
  output[i,15] <- range_upper2
}

write.csv(output, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/demo_suppl_table/cognitive.csv", row.names = F)

#########################################################

## IMAGING

# "Site"                          "Edited"
# [4301] "Batch"                         "Estimated_ICV"
# [4303] "WBV_With_Ventricles"           "WBV_No_Ventricles"
# [4305] "Global_GM_Volume"              "Cerebrum_WM_Volume"
# [4307] "Cerebellum_WM_Volume_Left"     "Cerebellum_WM_Volume_Right"
# [4309] "gFA"                           "gMD"
# [4311] "Reviewer"                      "Kappa_Final"
# [4313] "WMH_Volume_Total"              "WMH_Count"
# [4315] "Brain_age"                     "Fazekas_Score_Periventricular"
# [4317] "Fazekas_Score_Deep"            "Fazekas_Score_Total"
# [4319] "brain_accel"


# List phenotypes - do continuous imaging variables 
num <- c("Brain_age", "brain_accel", "Global_GM_Volume", "WBV_No_Ventricles", "gFA", "gMD", "WMH_Volume_Total", "Estimated_ICV")

# Create a numeric phenotype table  
output <- matrix(nrow = 1*length(num), ncol = 15)
output <- as.data.frame(output)
names(output) <- c("Phenotype", "N", "P", "Mean_1065", "(SD)_1065",  "range_lower", "range_upper", "IQR", "N778", "P778", "Mean_778", "(SD)_778",  "range778_lower", "range778_upper", "IQR778")

for (i in 1:length(num)){
  # for the 1065
  trait <- as.character(num[i])
  dataset <- comp[,which(colnames(comp) %in% c("GS_id", trait))]
  names(dataset)[2] <- "V2"
  dataset$V2 <- as.numeric(dataset$V2)
  dataset <- na.omit(dataset)
  mean <- mean(dataset$V2)
  sd <- sd(dataset$V2)
  N <- dim(dataset)[1]
  P <- ((dim(dataset)[1]) / 1065 )* 100
  IQR <- IQR(dataset$V2)[1]
  range_lower <- range(dataset$V2)[1]
  range_upper <- range(dataset$V2)[2]

  # for the 778 
  dataset2 <- MWAS[,which(colnames(MWAS) %in% c("GS_id", trait))]
  names(dataset2)[2] <- "V2"
  dataset2$V2 <- as.numeric(dataset2$V2)
  dataset2 <- na.omit(dataset2)
  N2 <- dim(dataset2)[1]
  P2 <- ((dim(dataset2)[1]) / 778 )* 100
  mean2 <- mean(dataset2$V2)
  sd2 <- sd(dataset2$V2)
  IQR2 <- IQR(dataset2$V2)[1]
  range_lower2 <- range(dataset2$V2)[1]
  range_upper2 <- range(dataset2$V2)[2]

  # # round 
  # P <- round(P, digits = 1)
  # P2 <- round(P2, digits = 1)
  # IQR <- round(IQR, digits = 1)
  # IQR2 <- round(IQR2, digits = 1)
  # mean <- round(mean, digits = 1)
  # mean2 <- round(mean2, digits = 1)
  # sd2 <- round(sd2, digits = 1)
  # sd <- round(sd, digits = 1)
  # range_lower <- round(range_lower, digits = 1)
  # range_upper <- round(range_upper, digits = 1)
  # range_lower2 <- round(range_lower2, digits = 1)
  # range_upper2 <- round(range_upper2, digits = 1)

  # output
  output[i,1] <- trait
  output[i,2] <- N
  output[i,3] <- P
  output[i,4] <- mean
  output[i,5] <- sd
  output[i,6] <- range_lower
  output[i,7] <- range_upper
  output[i,8] <- IQR
  output[i,9] <- N2
  output[i,10] <- P2
  output[i,11] <- mean2
  output[i,12] <- sd2
  output[i,13] <- range_lower2
  output[i,14] <- range_upper2
  output[i,15] <- IQR2
}

write.csv(output, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/demo_suppl_table/imaging_V2_unrounded.csv", row.names = F)


## FAZ scores 

table(comp$Fazekas_Score_Total)

#   1   2   3   4   5
# 159 503 194  69  15

table(MWAS$Fazekas_Score_Total)

#  1   2   3   4   5
# 107 379 152  50  12

## APOE

table(is.na(comp$apoe))

table(comp$apoe)

# e2e2 e2e3 e2e4 e3e3 e3e4 e4e4
#    5  121   22  633  234   35

table(MWAS$apoe)

# e2e2 e2e3 e2e4 e3e3 e3e4 e4e4
#    4   89   14  460  171   28


# DEPRESSION
table(comp$DiagnosisGiven)

table(MWAS$DiagnosisGiven)


# BMI 

table(is.na(comp$BMI)) # FALSE
mean(comp$BMI, na.rm = T)
sd(comp$BMI, na.rm = T)
IQR(comp$BMI)[1]
range(comp$BMI)[1]
range(comp$BMI)[2]

mean(MWAS$BMI, na.rm = T)
sd(MWAS$BMI, na.rm = T)
IQR(MWAS$BMI)[1]
range(MWAS$BMI)[1]
range(MWAS$BMI)[2]

### Read in EWAS phenotypes file to do the following markers:
MWAS <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/demo_suppl_table/MWAS_phenotypes_file.csv")

# EPISMOKER
mean(MWAS$smokingScore, na.rm = T)
sd(MWAS$smokingScore, na.rm = T)
IQR(MWAS$smokingScore)[1]
range(MWAS$smokingScore)[1]
range(MWAS$smokingScore)[2]

## CELL ESTIMATES

# "Eos"          "Neu"          "Bcell"
# [11] "CD4T"         "CD8T"         "Gran"         "Mono"         "NK"

# List phenotypes - do continuous imaging variables 
num <- c("CD4T", "CD8T", "Gran", "Bcell", "Mono", "NK", "Eos", "Neu")

# Create a numeric phenotype table  
output <- matrix(nrow = 1*length(num), ncol = 15)
output <- as.data.frame(output)
names(output) <- c("Phenotype", "N778", "P778", "Mean_778", "(SD)_778",  "range778_lower", "range778_upper", "IQR778")

for (i in 1:length(num)){
  # for the 1065
  trait <- as.character(num[i])

  # for the 778 
  dataset2 <- MWAS[,which(colnames(MWAS) %in% c("GS_id", trait))]
  names(dataset2)[2] <- "V2"
  dataset2$V2 <- as.numeric(dataset2$V2)
  dataset2 <- na.omit(dataset2)
  N2 <- dim(dataset2)[1]
  P2 <- ((dim(dataset2)[1]) / 778 )* 100
  mean2 <- mean(dataset2$V2)
  sd2 <- sd(dataset2$V2)
  IQR2 <- IQR(dataset2$V2)[1]
  range_lower2 <- range(dataset2$V2)[1]
  range_upper2 <- range(dataset2$V2)[2]

  # round 
  P2 <- round(P2, digits = 2)
  IQR2 <- round(IQR2, digits = 2)
  mean2 <- round(mean2, digits = 2)
  sd2 <- round(sd2, digits = 2)
  range_lower2 <- round(range_lower2, digits = 2)
  range_upper2 <- round(range_upper2, digits = 2)

  # output
  output[i,1] <- trait
  output[i,9] <- N2
  output[i,10] <- P2
  output[i,11] <- mean2
  output[i,12] <- sd2
  output[i,13] <- range_lower2
  output[i,14] <- range_upper2
  output[i,15] <- IQR2
}

write.csv(output, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/demo_suppl_table/immune_cell_estimates.csv", row.names = F)






