
###################################################################################

### Present results from comparisons in cox models 

###################################################################################

screen

R

library(ggplot2)
library(tidyverse)

# Read in full results files 

comb <- read.csv("/Local_Scratch/Danni/GDFBNP/00_Updated_results/Results_scores_cox/basic_v2.csv")

comb_full <- read.csv("/Local_Scratch/Danni/GDFBNP/00_Updated_results/Results_scores_cox/full_v2.csv")

# Do FDR correction
comb$FDR <- p.adjust(comb$P.Value, method = "BH")

# Keep any that meet significance threshold 
keep1 <- comb[which(comb$FDR < 0.05),] # 294 to 225 associations 
keep = paste(keep1$Predictor, keep1$Outcome, sep = "_")

# Keep only those passing the threshold from the basic model in the fully adjusted 
comb_full$retain <- paste(comb_full$Predictor, comb_full$Outcome, sep = "_")
comb_full$col <- ifelse(comb_full$retain %in% keep & comb_full$P.Value < 0.05, "red", "black")

write.csv(comb_full, "/Local_Scratch/Danni/GDFBNP/00_Updated_results/Results_scores_cox/comb_full_colours_EPIC_v2.csv", row.names = F)

write.csv(comb, "/Local_Scratch/Danni/GDFBNP/00_Updated_results/Results_scores_cox/comb_basic_FDR_v2.csv", row.names = F)

###############

## Make suppl table with all results

basic <- read.csv("/Local_Scratch/Danni/GDFBNP/00_Updated_results/Results_scores_cox/comb_basic_FDR_v2.csv")
full <- read.csv("/Local_Scratch/Danni/GDFBNP/00_Updated_results/Results_scores_cox/full_v2.csv")
WBC <- read.csv("/Local_Scratch/Danni/GDFBNP/00_Updated_results/Results_scores_cox/WBC_v2.csv")

basic <- basic[c(1,2,7:11,3:6, 12)]

basic <- basic[order(basic$P.Value),]
basic$retain <- paste(basic$Predictor, basic$Outcome, sep = "_")
full$retain <- paste(full$Predictor, full$Outcome, sep = "_")
WBC$retain <- paste(WBC$Predictor, WBC$Outcome, sep = "_")

join <- left_join(basic, full, by = 'retain')
join <- left_join(join, WBC, by = 'retain')

# Adjust names
join <- join[c(1:12, 20:24, 16:19, 27:30)]

#  [1] "Predictor.x"       "Outcome.x"         "No..of.Cases.x"
#  [4] "No..of.Controls.x" "cox.zph.x"         "tte_mean_sd.x"
#  [7] "max_tte"           "Hazard.Ratio.x"    "LCI.x"
# [10] "UCI.x"             "P.Value.x"         "FDR"
# [13] "No..of.Cases.y"    "No..of.Controls.y" "cox.zph.y"
# [16] "tte_mean_sd.y"     "NA..x"             "Hazard.Ratio.y"
# [19] "LCI.y"             "UCI.y"             "P.Value.y"
# [22] "Hazard.Ratio"      "LCI"               "UCI"
# [25] "P.Value"

names(join) <- c("Marker", "Outcome", "Basic N Cases", 
  "Basic N Controls", 'Basic zph', "Basic mean tte (sd)", 
  "Basic max tte", "Basic HR", "Basic LCI", 
  "Basic UCI", "Basic P", "Basic FDR P", 
  "Full N Cases", "Full N Controls", 'NA',
  "Full mean tte (sd)", "Full max tte","Full HR",
   "Full LCI", "Full UCI", "Full P", 
   "WBC HR", "WBC_LCI", "WBC UCI",
  "WBC P")
write.csv(join, "/Local_Scratch/Danni/GDFBNP/00_Updated_results/Results_scores_cox/Suppl_table_results.csv", row.names = F)


###############

res <- comb_full

# Separate to groups of measurements 
prot_GDF <- res[which(res$Predictor %in% c("gdf15_rnk")),]
prot_BNP <- res[which(res$Predictor %in% c("nt.probnp_rnk")),]
score_GDF <- res[which(res$Predictor %in% c("GDF_score")),]
score_BNP <- res[which(res$Predictor %in% c("BNP_score")),]


# Add identifiers for whether we have scores or proteins
prot_GDF$Biomarker <- "GDF15"
prot_BNP$Biomarker <- "Nt-proBNP"

score_GDF$Biomarker <- "EpiScore"
score_BNP$Biomarker <- "EpiScore"

# Add plot ready naming 
prot_GDF$Name <- "GDF15"
prot_BNP$Name <- "NT-proBNP"

score_GDF$Name <- "GDF15"
score_BNP$Name <- "NT-proBNP"

# Join together for each biomarker 
BNP <- rbind(prot_BNP, score_BNP)

GDF <- rbind(prot_GDF, score_GDF)

# # Make a unique category for each
# BNP$unique <- paste(BNP$Biomarker, BNP$Outcome, sep = "_")
# GDF$unique <- paste(GDF$Biomarker, GDF$Outcome, sep = "_")


library(stringr)

## Plot HR comparison for BNP

x <- BNP

x$Outcome <- str_replace(x$Outcome, "IHD", "Ischaemic Heart Disease")
x$Outcome <- str_replace(x$Outcome, "Diabetes", "Type 2 Diabetes")
x$Outcome <- str_replace(x$Outcome, "Stroke", "Ischaemic Stroke")

# x$Outcome2 <- x$Outcome
# x$TraitVar <- paste0(x$Name)
# x$TraitVar = factor(x$TraitVar, levels=unique(x$TraitVar[rev(order(x$Hazard.Ratio.x))]))

My_Theme = theme(
  panel.border = element_rect(colour="black",size=1, fill = NA),
  axis.title.x = element_text(size = 20), # controls HR label size 
  axis.text.x = element_text(size = 20),
  axis.text.y = element_text(size = 20),
  axis.title.y = element_text(size = 20),
  strip.text = element_text(size = 20, face = "bold"),
  legend.text=element_text(size=20),
  legend.title=element_text(size=20, face = "bold"), legend.position = "none",
  axis.title=element_text(size=20))

# collist1 <- c("blue", "black", "green", "purple", "pink", "grey", "orange", "yellow")
# collist2 <- c("blue", "black", "green", "purple", "pink", "grey", "orange", "yellow")

collist1 <- c("black", "red", "red", "black", "black", "black", "red", "black")
collist2 <- c("black", "black", "black", "black", "red", "black", "red", "red")

x$Biomarker <- str_replace(x$Biomarker, "Nt-proBNP", "NT-proBNP")

pdf("/Local_Scratch/Danni/GDFBNP/00_Updated_results/Results_scores_cox/BNP_plot_cox_EPIC_v3.pdf", width = 18, height = 5)
ggplot(x,aes(y=Hazard.Ratio, x=Biomarker)) + 
  geom_point(size = 4.5, colour = collist1)+
  geom_errorbar(aes(ymin = LCI, ymax = UCI),
                position = position_dodge(0.5), width = 0.1,
                colour = collist2)+
ylab("Hazard Ratio (95% CI)")+ xlab ("") + theme_classic() +
  geom_hline(yintercept = 1, linetype = "dotted")+
  theme(axis.text.x = element_text(size = 12, vjust = 0.5), axis.text.y = element_text(size = 12), legend.position = "right",
        plot.title = element_text(size = 12))+ theme(legend.title = element_text(hjust = 0.5)) +
  coord_flip() + facet_wrap(~Outcome, scales = "free_y") + My_Theme
  dev.off()

a <- ggplot(x,aes(y=Hazard.Ratio, x=Biomarker)) + 
  geom_point(size = 4.5, colour = collist1)+
  geom_errorbar(aes(ymin = LCI, ymax = UCI),
                position = position_dodge(0.5), width = 0.1,
                colour = collist2)+
ylab("Hazard Ratio (95% CI)")+ xlab ("") + theme_classic() +
  geom_hline(yintercept = 1, linetype = "dotted")+
  theme(axis.text.x = element_text(size = 12, vjust = 0.5), axis.text.y = element_text(size = 12), legend.position = "right",
        plot.title = element_text(size = 12))+ theme(legend.title = element_text(hjust = 0.5)) +
  coord_flip() + facet_wrap(~Outcome, scales = "free_y") + My_Theme

## Plot HR comparison for GDF

x <- GDF

x$Outcome <- str_replace(x$Outcome, "IHD", "Ischaemic Heart Disease")
x$Outcome <- str_replace(x$Outcome, "Diabetes", "Type 2 Diabetes")
x$Outcome <- str_replace(x$Outcome, "Stroke", "Ischaemic Stroke")


# x$Outcome2 <- x$Outcome
# x$TraitVar <- paste0(x$Name)
# x$TraitVar = factor(x$TraitVar, levels=unique(x$TraitVar[rev(order(x$Hazard.Ratio.x))]))

My_Theme = theme(
  panel.border = element_rect(colour="black",size=1, fill = NA),
  axis.title.x = element_text(size = 20), # controls HR label size 
  axis.text.x = element_text(size = 20),
  axis.text.y = element_text(size = 20),
  axis.title.y = element_text(size = 20),
  strip.text = element_text(size = 20, face = "bold"),
  legend.text=element_text(size=20),
  legend.title=element_text(size=20, face = "bold"), legend.position = "none",
  axis.title=element_text(size=20))

# collist1 <- c("blue", "black", "green", "purple", "pink", "grey", "orange", "yellow")
# collist2 <- c("blue", "black", "green", "purple", "pink", "grey", "orange", "yellow")

collist1 <- c("red", "red", "red", "black", "red", "red", "red", "black")
collist2 <- c("red", "red", "black", "black", "red", "red", "red", "red")

pdf("/Local_Scratch/Danni/GDFBNP/00_Updated_results/Results_scores_cox/GDF_plot_cox_EPIC_v3.pdf", width = 18, height = 5)
ggplot(x,aes(y=Hazard.Ratio, x=Biomarker)) + 
  geom_point(size = 4.5, colour = collist1)+
  geom_errorbar(aes(ymin = LCI, ymax = UCI),
                position = position_dodge(0.5), width = 0.1,
                colour = collist2)+
ylab("Hazard Ratio (95% CI)")+ xlab ("") + theme_classic() +
  geom_hline(yintercept = 1, linetype = "dotted")+
  theme(axis.text.x = element_text(size = 12, vjust = 0.5), axis.text.y = element_text(size = 12), legend.position = "right",
        plot.title = element_text(size = 12))+ theme(legend.title = element_text(hjust = 0.5)) +
  coord_flip() + facet_wrap(~Outcome, scales = "free_y") + My_Theme
  dev.off()

b <- ggplot(x,aes(y=Hazard.Ratio, x=Biomarker)) + 
  geom_point(size = 4.5, colour = collist1)+
  geom_errorbar(aes(ymin = LCI, ymax = UCI),
                position = position_dodge(0.5), width = 0.1,
                colour = collist2)+
ylab("Hazard Ratio (95% CI)")+ xlab ("") + theme_classic() +
  geom_hline(yintercept = 1, linetype = "dotted")+
  theme(axis.text.x = element_text(size = 12, vjust = 0.5), axis.text.y = element_text(size = 12), legend.position = "right",
        plot.title = element_text(size = 12))+ theme(legend.title = element_text(hjust = 0.5)) +
  coord_flip() + facet_wrap(~Outcome, scales = "free_y") + My_Theme

### Save patchwork plots
library(patchwork)
pdf("/Local_Scratch/Danni/GDFBNP/00_Updated_results/Results_scores_cox/JOINT_plot_cox_EPIC_v3.pdf", width = 18, height =12)
b / plot_spacer() / a + plot_layout(heights = c(4, 0.5 ,4))
dev.off()

# #############################

# ### PLOT originaly protein assocs in the 20k 


# screen

# R

# library(ggplot2)

# # read in results 20k proteins 
# comb <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_Cox_proteins_20k/Cox_proteins_210522/Cox_results/basic_210522.csv")
# comb_full <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_Cox_proteins_20k/Cox_proteins_210522/Cox_results/full_210522.csv")

# # Do FDR correction
# comb$FDR <- p.adjust(comb$P.Value, method = "BH")

# # Keep any that meet significance threshold 
# keep1 <- comb[which(comb$FDR < 0.05),] # 294 to 225 associations 
# keep = paste(keep1$Predictor, keep1$Outcome, sep = "_")

# # Keep only those passing the threshold from the basic model in the fully adjusted 
# comb_full$retain <- paste(comb_full$Predictor, comb_full$Outcome, sep = "_")
# comb_full$col <- ifelse(comb_full$retain %in% keep & comb_full$P.Value < 0.05, "red", "black")

# ## Plot HR comparison for GDF

# x <- comb_full

# library(stringr)
# x$Outcome <- str_replace(x$Outcome, "Alzheimer's Disease", "Alzheimer's Dementia")
# x$Outcome <- str_replace(x$Outcome, "IHD", "Ischaemic Heart Disease")
# x$Outcome <- str_replace(x$Outcome, "Diabetes", "Type 2 Diabetes")
# x$Outcome <- str_replace(x$Outcome, "Stroke", "Ischaemic Stroke")

# x$Predictor <- str_replace(x$Predictor, "gdf15_rnk", "GDF15")
# x$Predictor <- str_replace(x$Predictor, "nt.probnp_rnk", "Nt-proBNP")

# # x$Outcome2 <- x$Outcome
# # x$TraitVar <- paste0(x$Name)
# # x$TraitVar = factor(x$TraitVar, levels=unique(x$TraitVar[rev(order(x$Hazard.Ratio.x))]))

# My_Theme = theme(
#   panel.border = element_rect(colour="black",size=1, fill = NA),
#   axis.title.x = element_text(size = 20), # controls HR label size 
#   axis.text.x = element_text(size = 20),
#   axis.text.y = element_text(size = 20),
#   axis.title.y = element_text(size = 20),
#   strip.text = element_text(size = 20, face = "bold"),
#   legend.text=element_text(size=20),
#   legend.title=element_text(size=20, face = "bold"), legend.position = "none",
#   axis.title=element_text(size=20))


# # collist1 <- c("red", "blue", "black", "green", "purple", "pink", "grey", "orange", "yellow", "darkblue")
# collist1 <- c("black", "black", "red", "red", "red", "red", "red", "red", "black", "red")


# # collist2 <- c("red", "blue", "black", "green", "purple", "pink", "grey", "orange", "yellow", "darkblue")
# collist2 <- c("black", "black", "black", "red", "red", "red", "red", "red", "red", "red")



# pdf("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/Plots/Cox_ggplots/20k.pdf", width = 18, height = 5)
# ggplot(x,aes(y=Hazard.Ratio, x=Predictor)) + 
#   geom_point(size = 4.5, colour = collist1)+
#   geom_errorbar(aes(ymin = LCI, ymax = UCI),
#                 position = position_dodge(0.5), width = 0.1,
#                 colour = collist2)+
# ylab("Hazard Ratio (95% CI)")+ xlab ("") + theme_classic() +
#   geom_hline(yintercept = 1, linetype = "dotted")+
#   theme(axis.text.x = element_text(size = 12, vjust = 0.5), axis.text.y = element_text(size = 12), legend.position = "right",
#         plot.title = element_text(size = 12))+ theme(legend.title = element_text(hjust = 0.5)) +
#   coord_flip() + facet_wrap(~Outcome, scales = "free_y") + My_Theme
#   dev.off()

# b <- ggplot(x,aes(y=Hazard.Ratio, x=Biomarker)) + 
#   geom_point(size = 4.5, colour = collist1)+
#   geom_errorbar(aes(ymin = LCI, ymax = UCI),
#                 position = position_dodge(0.5), width = 0.1,
#                 colour = collist2)+
# ylab("Hazard Ratio (95% CI)")+ xlab ("") + theme_classic() +
#   geom_hline(yintercept = 1, linetype = "dotted")+
#   theme(axis.text.x = element_text(size = 12, vjust = 0.5), axis.text.y = element_text(size = 12), legend.position = "right",
#         plot.title = element_text(size = 12))+ theme(legend.title = element_text(hjust = 0.5)) +
#   coord_flip() + facet_wrap(~Outcome, scales = "free_y") + My_Theme


# # # Separate to groups of measurements 
# # prot_GDF <- res[which(res$Predictor %in% c("gdf15_rnk")),]
# # prot_BNP <- res[which(res$Predictor %in% c("nt.probnp_rnk")),]

# # # Add identifiers for whether we have scores or proteins
# # prot_GDF$Biomarker <- "GDF15"
# # prot_BNP$Biomarker <- "Nt-pro-BNP"

# # # Add plot ready naming 
# # prot_GDF$Name <- "GDF15"
# # prot_BNP$Name <- "Nt-pro-BNP"

# # BNP <- prot_BNP
# # GDF <- prot_GDF

# # # Make a unique category for each
# # BNP$unique <- paste(BNP$Biomarker, BNP$Outcome, sep = "_")
# # GDF$unique <- paste(GDF$Biomarker, GDF$Outcome, sep = "_")

# # x <- rbind(GDF, BNP)

# # ## Plot HR comparison for BNP

# # # x <- BNP

# # # x$Outcome2 <- x$Outcome
# # # x$TraitVar <- paste0(x$Name)
# # # x$TraitVar = factor(x$TraitVar, levels=unique(x$TraitVar[rev(order(x$Hazard.Ratio.x))]))

# # My_Theme = theme(
# #   panel.border = element_rect(colour="black",size=1, fill = NA),
# #   axis.title.x = element_text(size = 20), # controls HR label size 
# #   axis.text.x = element_text(size = 20),
# #   axis.text.y = element_text(size = 20),
# #   axis.title.y = element_text(size = 20),
# #   strip.text = element_text(size = 20, face = "bold"),
# #   legend.text=element_text(size=20),
# #   legend.title=element_text(size=20, face = "bold"), legend.position = "none",
# #   axis.title=element_text(size=20))

# # x <- x[-which(x$Outcome %in% "IHD"),]

# # pdf("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_Cox_proteins_20k/Results/Proteins_plot_diseases_cox_20k.pdf", width = 18, height = 5)
# # ggplot(x,aes(y=Hazard.Ratio, x=Biomarker)) + 
# #   geom_point(size = 4.5, colour = "dodgerblue4")+
# #   geom_errorbar(aes(ymin = LCI, ymax = UCI),
# #                 position = position_dodge(0.5), width = 0.1,
# #                 colour = "dodgerblue4")+
# # ylab("Hazard Ratio (95% CI) - FDR P<0.05")+ xlab ("") + theme_classic() +
# #   geom_hline(yintercept = 1, linetype = "dotted")+
# #   theme(axis.text.x = element_text(size = 12, vjust = 0.5), axis.text.y = element_text(size = 12), legend.position = "right",
# #         plot.title = element_text(size = 12))+ theme(legend.title = element_text(hjust = 0.5)) +
# #   coord_flip() + facet_wrap(~Outcome, scales = "free_y") + My_Theme
# #   dev.off()



