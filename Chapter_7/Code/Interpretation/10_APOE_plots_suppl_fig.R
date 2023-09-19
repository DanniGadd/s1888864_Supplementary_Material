
###################################################################################################

### Plot by haplotype for top proteins for APOE 

###################################################################################################

screen 

R

library(tidyverse)
library(readxl)

### APOE 
apoe <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS_run2/APOE/APOE_processed/APOE_JOINT_annotated_THR.csv")
apoe$type <- "APOE"
list1 <- apoe[which(apoe$Status == "pass"),]
list1 <- list1[c("SeqId", "phenotype", "type", "Gene.Name.Name", "beta", "SE", "Pcalc", "UniProt.Full.Name")] # 16 proteins for APOE 

data <- list1[c(1,4)] # get just seqID and gene names for proteins 

data <- data %>% unique() # 11 unique protein levels associated with either e2 or e4 allele status 

# Read in the protein data 
prot <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS/01_Phenotype_collation/prot_file_220221.csv", check.names = F)

# Plot for each of the proteins in relation to the main protein measurements and haplotypes 
library(ggpubr)
library(ggplot2)

pred <- data

plot_list <- list()
plot_list2 <- list()

### Individual plots, so that we can patch work them together as one figure 

# List APOE proteins 
# data 
#        SeqId Gene.Name.Name
# 1   11293-14          LRRN1
# 2   12501-10           TBCA
# 3  10082-251           NEFL
# 4   17671-58           ING4
# 5    7223-60        S100A13
# 6   10046-55          BIRC2
# 7    19158-1            PAF
# 8     6378-2        C5orf38
# 9    5744-12           MENT
# 10    8922-4          TMCC3
# 11   2797-56           APOB
# 12   4337-49            CRP
# 13   19272-9           PEF1
# 14   6433-57         FAM20A


# 2 TBCA
i <- 2
seqid_name <- as.character(data[i,1])
  protein_name <- as.character(data[i,2])

  d1 <- prot[,"apoe"] %>% as.data.frame()
  d2 <- prot[,seqid_name] %>% as.data.frame()
  d3 <- cbind(d1,d2)
  names(d3) <- c("apoe", "prot")
  d3 <- na.omit(d3)

  p = ggplot(d3, aes(x=apoe , y=prot)) +
  geom_violin(trim=FALSE, fill='#A4A4A4', color="darkred")+
  geom_boxplot(width=0.1) + theme_minimal() + ggtitle(protein_name)

TBCA <-  ggplot(d3, aes(x=apoe , y=prot)) +
  geom_violin(trim=FALSE, fill='lightskyblue3', color = 'lightskyblue3')+
  geom_boxplot(width=0.1) + theme_minimal() + xlab("APOE Haplotype") + ylab(protein_name) + theme_classic() + 
  theme(legend.title = element_blank(),
    legend.position = "none", 
    axis.title.x = element_text(size = 27), 
    axis.text.x = element_text(size = 27),
    axis.text.y = element_text(size = 27),
    axis.title.y = element_text(size = 27),
    plot.title = element_text(size = 27),
    axis.title.x.bottom = element_text(margin = margin(14, 0, 0, 0)))




# 3 NEFL
i <- 3
seqid_name <- as.character(data[i,1])
  protein_name <- as.character(data[i,2])

  d1 <- prot[,"apoe"] %>% as.data.frame()
  d2 <- prot[,seqid_name] %>% as.data.frame()
  d3 <- cbind(d1,d2)
  names(d3) <- c("apoe", "prot")
  d3 <- na.omit(d3)

  p = ggplot(d3, aes(x=apoe , y=prot)) +
  geom_violin(trim=FALSE, fill='#A4A4A4', color="darkred")+
  geom_boxplot(width=0.1) + theme_minimal() + ggtitle(protein_name)

NEFL <-  ggplot(d3, aes(x=apoe , y=prot)) +
  geom_violin(trim=FALSE, fill='lightskyblue3', color = 'lightskyblue3')+
  geom_boxplot(width=0.1) + theme_minimal() + xlab("APOE Haplotype") + ylab(protein_name) + theme_classic() + 
  theme(legend.title = element_blank(),
    legend.position = "none", 
    axis.title.x = element_text(size = 27), 
    axis.text.x = element_text(size = 27),
    axis.text.y = element_text(size = 27),
    axis.title.y = element_text(size = 27),
    plot.title = element_text(size = 27),
    axis.title.x.bottom = element_text(margin = margin(14, 0, 0, 0)))




# 4 ING4
i <- 4
seqid_name <- as.character(data[i,1])
  protein_name <- as.character(data[i,2])

  d1 <- prot[,"apoe"] %>% as.data.frame()
  d2 <- prot[,seqid_name] %>% as.data.frame()
  d3 <- cbind(d1,d2)
  names(d3) <- c("apoe", "prot")
  d3 <- na.omit(d3)

  p = ggplot(d3, aes(x=apoe , y=prot)) +
  geom_violin(trim=FALSE, fill='#A4A4A4', color="darkred")+
  geom_boxplot(width=0.1) + theme_minimal() + ggtitle(protein_name)

ING4 <-  ggplot(d3, aes(x=apoe , y=prot)) +
  geom_violin(trim=FALSE, fill='lightskyblue3', color = 'lightskyblue3')+
  geom_boxplot(width=0.1) + theme_minimal() + xlab("APOE Haplotype") + ylab(protein_name) + theme_classic() + 
  theme(legend.title = element_blank(),
    legend.position = "none", 
    axis.title.x = element_text(size = 27), 
    axis.text.x = element_text(size = 27),
    axis.text.y = element_text(size = 27),
    axis.title.y = element_text(size = 27),
    plot.title = element_text(size = 27),
    axis.title.x.bottom = element_text(margin = margin(14, 0, 0, 0)))




# 5 s100A13
i <- 5
seqid_name <- as.character(data[i,1])
  protein_name <- as.character(data[i,2])

  d1 <- prot[,"apoe"] %>% as.data.frame()
  d2 <- prot[,seqid_name] %>% as.data.frame()
  d3 <- cbind(d1,d2)
  names(d3) <- c("apoe", "prot")
  d3 <- na.omit(d3)

  p = ggplot(d3, aes(x=apoe , y=prot)) +
  geom_violin(trim=FALSE, fill='#A4A4A4', color="darkred")+
  geom_boxplot(width=0.1) + theme_minimal() + ggtitle(protein_name)

S100A13 <-  ggplot(d3, aes(x=apoe , y=prot)) +
  geom_violin(trim=FALSE, fill='lightskyblue3', color = 'lightskyblue3')+
  geom_boxplot(width=0.1) + theme_minimal() + xlab("APOE Haplotype") + ylab(protein_name) + theme_classic() + 
  theme(legend.title = element_blank(),
    legend.position = "none", 
    axis.title.x = element_text(size = 27), 
    axis.text.x = element_text(size = 27),
    axis.text.y = element_text(size = 27),
    axis.title.y = element_text(size = 27),
    plot.title = element_text(size = 27),
    axis.title.x.bottom = element_text(margin = margin(14, 0, 0, 0)))


# 6 BIRC2
i <- 6
seqid_name <- as.character(data[i,1])
  protein_name <- as.character(data[i,2])

  d1 <- prot[,"apoe"] %>% as.data.frame()
  d2 <- prot[,seqid_name] %>% as.data.frame()
  d3 <- cbind(d1,d2)
  names(d3) <- c("apoe", "prot")
  d3 <- na.omit(d3)

  p = ggplot(d3, aes(x=apoe , y=prot)) +
  geom_violin(trim=FALSE, fill='#A4A4A4', color="darkred")+
  geom_boxplot(width=0.1) + theme_minimal() + ggtitle(protein_name)

BIRC2 <-  ggplot(d3, aes(x=apoe , y=prot)) +
  geom_violin(trim=FALSE, fill='lightskyblue3', color = 'lightskyblue3')+
  geom_boxplot(width=0.1) + theme_minimal() + xlab("APOE Haplotype") + ylab(protein_name) + theme_classic() + 
  theme(legend.title = element_blank(),
    legend.position = "none", 
    axis.title.x = element_text(size = 27), 
    axis.text.x = element_text(size = 27),
    axis.text.y = element_text(size = 27),
    axis.title.y = element_text(size = 27),
    plot.title = element_text(size = 27),
    axis.title.x.bottom = element_text(margin = margin(14, 0, 0, 0)))


# 7 PAF
i <- 7
seqid_name <- as.character(data[i,1])
  protein_name <- as.character(data[i,2])

  d1 <- prot[,"apoe"] %>% as.data.frame()
  d2 <- prot[,seqid_name] %>% as.data.frame()
  d3 <- cbind(d1,d2)
  names(d3) <- c("apoe", "prot")
  d3 <- na.omit(d3)

  p = ggplot(d3, aes(x=apoe , y=prot)) +
  geom_violin(trim=FALSE, fill='#A4A4A4', color="darkred")+
  geom_boxplot(width=0.1) + theme_minimal() + ggtitle(protein_name)

PAF <-  ggplot(d3, aes(x=apoe , y=prot)) +
  geom_violin(trim=FALSE, fill='lightskyblue3', color = 'lightskyblue3')+
  geom_boxplot(width=0.1) + theme_minimal() + xlab("APOE Haplotype") + ylab(protein_name) + theme_classic() + 
  theme(legend.title = element_blank(),
    legend.position = "none", 
    axis.title.x = element_text(size = 27), 
    axis.text.x = element_text(size = 27),
    axis.text.y = element_text(size = 27),
    axis.title.y = element_text(size = 27),
    plot.title = element_text(size = 27),
    axis.title.x.bottom = element_text(margin = margin(14, 0, 0, 0)))




# 10  TMCC3
i <- 10
seqid_name <- as.character(data[i,1])
  protein_name <- as.character(data[i,2])

  d1 <- prot[,"apoe"] %>% as.data.frame()
  d2 <- prot[,seqid_name] %>% as.data.frame()
  d3 <- cbind(d1,d2)
  names(d3) <- c("apoe", "prot")
  d3 <- na.omit(d3)

  p = ggplot(d3, aes(x=apoe , y=prot)) +
  geom_violin(trim=FALSE, fill='#A4A4A4', color="darkred")+
  geom_boxplot(width=0.1) + theme_minimal() + ggtitle(protein_name)

TMCC3 <-  ggplot(d3, aes(x=apoe , y=prot)) +
  geom_violin(trim=FALSE, fill='lightskyblue3', color = 'lightskyblue3')+
  geom_boxplot(width=0.1) + theme_minimal() + xlab("APOE Haplotype") + ylab(protein_name) + theme_classic() + 
  theme(legend.title = element_blank(),
    legend.position = "none", 
    axis.title.x = element_text(size = 27), 
    axis.text.x = element_text(size = 27),
    axis.text.y = element_text(size = 27),
    axis.title.y = element_text(size = 27),
    plot.title = element_text(size = 27),
    axis.title.x.bottom = element_text(margin = margin(14, 0, 0, 0)))



# 12  CRP
i <- 12
seqid_name <- as.character(data[i,1])
  protein_name <- as.character(data[i,2])

  d1 <- prot[,"apoe"] %>% as.data.frame()
  d2 <- prot[,seqid_name] %>% as.data.frame()
  d3 <- cbind(d1,d2)
  names(d3) <- c("apoe", "prot")
  d3 <- na.omit(d3)

  p = ggplot(d3, aes(x=apoe , y=prot)) +
  geom_violin(trim=FALSE, fill='#A4A4A4', color="darkred")+
  geom_boxplot(width=0.1) + theme_minimal() + ggtitle(protein_name)

CRP <-  ggplot(d3, aes(x=apoe , y=prot)) +
  geom_violin(trim=FALSE, fill='lightskyblue3', color = 'lightskyblue3')+
  geom_boxplot(width=0.1) + theme_minimal() + xlab("APOE Haplotype") + ylab(protein_name) + theme_classic() + 
  theme(legend.title = element_blank(),
    legend.position = "none", 
    axis.title.x = element_text(size = 27), 
    axis.text.x = element_text(size = 27),
    axis.text.y = element_text(size = 27),
    axis.title.y = element_text(size = 27),
    plot.title = element_text(size = 27),
    axis.title.x.bottom = element_text(margin = margin(14, 0, 0, 0)))


# 1
i <- 1
seqid_name <- as.character(data[i,1])
  protein_name <- as.character(data[i,2])

  d1 <- prot[,"apoe"] %>% as.data.frame()
  d2 <- prot[,seqid_name] %>% as.data.frame()
  d3 <- cbind(d1,d2)
  names(d3) <- c("apoe", "prot")
  d3 <- na.omit(d3)

  p = ggplot(d3, aes(x=apoe , y=prot)) +
  geom_violin(trim=FALSE, fill='#A4A4A4', color="darkred")+
  geom_boxplot(width=0.1) + theme_minimal() + ggtitle(protein_name)


LRRN1 <-  ggplot(d3, aes(x=apoe , y=prot)) +
  geom_violin(trim=FALSE, fill='firebrick2', color = 'firebrick2')+
  geom_boxplot(width=0.1) + theme_minimal() + xlab("APOE Haplotype") + ylab(protein_name) + theme_classic() + 
  theme(legend.title = element_blank(),
    legend.position = "none", 
    axis.title.x = element_text(size = 27), 
    axis.text.x = element_text(size = 27),
    axis.text.y = element_text(size = 27),
    axis.title.y = element_text(size = 27),
    plot.title = element_text(size = 27),
    axis.title.x.bottom = element_text(margin = margin(14, 0, 0, 0)))


# 8  C5orf38
i <- 8
seqid_name <- as.character(data[i,1])
  protein_name <- as.character(data[i,2])

  d1 <- prot[,"apoe"] %>% as.data.frame()
  d2 <- prot[,seqid_name] %>% as.data.frame()
  d3 <- cbind(d1,d2)
  names(d3) <- c("apoe", "prot")
  d3 <- na.omit(d3)

  p = ggplot(d3, aes(x=apoe , y=prot)) +
  geom_violin(trim=FALSE, fill='#A4A4A4', color="darkred")+
  geom_boxplot(width=0.1) + theme_minimal() + ggtitle(protein_name)

C5orf38 <-  ggplot(d3, aes(x=apoe , y=prot)) +
  geom_violin(trim=FALSE, fill='firebrick2', color = 'firebrick2')+
  geom_boxplot(width=0.1) + theme_minimal() + xlab("APOE Haplotype") + ylab(protein_name) + theme_classic() + 
  theme(legend.title = element_blank(),
    legend.position = "none", 
    axis.title.x = element_text(size = 27), 
    axis.text.x = element_text(size = 27),
    axis.text.y = element_text(size = 27),
    axis.title.y = element_text(size = 27),
    plot.title = element_text(size = 27),
    axis.title.x.bottom = element_text(margin = margin(14, 0, 0, 0)))


# 9  MENT
i <- 9
seqid_name <- as.character(data[i,1])
  protein_name <- as.character(data[i,2])

  d1 <- prot[,"apoe"] %>% as.data.frame()
  d2 <- prot[,seqid_name] %>% as.data.frame()
  d3 <- cbind(d1,d2)
  names(d3) <- c("apoe", "prot")
  d3 <- na.omit(d3)

  p = ggplot(d3, aes(x=apoe , y=prot)) +
  geom_violin(trim=FALSE, fill='#A4A4A4', color="darkred")+
  geom_boxplot(width=0.1) + theme_minimal() + ggtitle(protein_name)

MENT <-  ggplot(d3, aes(x=apoe , y=prot)) +
  geom_violin(trim=FALSE, fill='firebrick2', color = 'firebrick2')+
  geom_boxplot(width=0.1) + theme_minimal() + xlab("APOE Haplotype") + ylab(protein_name) + theme_classic() + 
  theme(legend.title = element_blank(),
    legend.position = "none", 
    axis.title.x = element_text(size = 27), 
    axis.text.x = element_text(size = 27),
    axis.text.y = element_text(size = 27),
    axis.title.y = element_text(size = 27),
    plot.title = element_text(size = 27),
    axis.title.x.bottom = element_text(margin = margin(14, 0, 0, 0)))



# 11  APOB
i <- 11
seqid_name <- as.character(data[i,1])
  protein_name <- as.character(data[i,2])

  d1 <- prot[,"apoe"] %>% as.data.frame()
  d2 <- prot[,seqid_name] %>% as.data.frame()
  d3 <- cbind(d1,d2)
  names(d3) <- c("apoe", "prot")
  d3 <- na.omit(d3)

  p = ggplot(d3, aes(x=apoe , y=prot)) +
  geom_violin(trim=FALSE, fill='#A4A4A4', color="darkred")+
  geom_boxplot(width=0.1) + theme_minimal() + ggtitle(protein_name)

APOB <-  ggplot(d3, aes(x=apoe , y=prot)) +
  geom_violin(trim=FALSE, fill='firebrick2', color = 'firebrick2')+
  geom_boxplot(width=0.1) + theme_minimal() + xlab("APOE Haplotype") + ylab(protein_name) + theme_classic() + 
  theme(legend.title = element_blank(),
    legend.position = "none", 
    axis.title.x = element_text(size = 27), 
    axis.text.x = element_text(size = 27),
    axis.text.y = element_text(size = 27),
    axis.title.y = element_text(size = 27),
    plot.title = element_text(size = 27),
    axis.title.x.bottom = element_text(margin = margin(14, 0, 0, 0)))


# 14 FAM20A
i <- 14
seqid_name <- as.character(data[i,1])
  protein_name <- as.character(data[i,2])

  d1 <- prot[,"apoe"] %>% as.data.frame()
  d2 <- prot[,seqid_name] %>% as.data.frame()
  d3 <- cbind(d1,d2)
  names(d3) <- c("apoe", "prot")
  d3 <- na.omit(d3)

  p = ggplot(d3, aes(x=apoe , y=prot)) +
  geom_violin(trim=FALSE, fill='#A4A4A4', color="darkred")+
  geom_boxplot(width=0.1) + theme_minimal() + ggtitle(protein_name)

FAM20A <-  ggplot(d3, aes(x=apoe , y=prot)) +
  geom_violin(trim=FALSE, fill='lightskyblue3', color = 'lightskyblue3')+
  geom_boxplot(width=0.1) + theme_minimal() + xlab("APOE Haplotype") + ylab(protein_name) + theme_classic() + 
  theme(legend.title = element_blank(),
    legend.position = "none", 
    axis.title.x = element_text(size = 27), 
    axis.text.x = element_text(size = 27),
    axis.text.y = element_text(size = 27),
    axis.title.y = element_text(size = 27),
    plot.title = element_text(size = 27),
    axis.title.x.bottom = element_text(margin = margin(14, 0, 0, 0)))


# 15  PEF1
i <- 13
seqid_name <- as.character(data[i,1])
  protein_name <- as.character(data[i,2])

  d1 <- prot[,"apoe"] %>% as.data.frame()
  d2 <- prot[,seqid_name] %>% as.data.frame()
  d3 <- cbind(d1,d2)
  names(d3) <- c("apoe", "prot")
  d3 <- na.omit(d3)

  p = ggplot(d3, aes(x=apoe , y=prot)) +
  geom_violin(trim=FALSE, fill='#A4A4A4', color="darkred")+
  geom_boxplot(width=0.1) + theme_minimal() + ggtitle(protein_name)

PEF1 <-  ggplot(d3, aes(x=apoe , y=prot)) +
  geom_violin(trim=FALSE, fill='firebrick2', color = 'firebrick2')+
  geom_boxplot(width=0.1) + theme_minimal() + xlab("APOE Haplotype") + ylab(protein_name) + theme_classic() + 
  theme(legend.title = element_blank(),
    legend.position = "none", 
    axis.title.x = element_text(size = 27), 
    axis.text.x = element_text(size = 27),
    axis.text.y = element_text(size = 27),
    axis.title.y = element_text(size = 27),
    plot.title = element_text(size = 27),
    axis.title.x.bottom = element_text(margin = margin(14, 0, 0, 0)))


library(patchwork)

pdf("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/heatmap_phewas/APOE.pdf", width = 35, height = 22)
TBCA + NEFL + ING4 + S100A13 + BIRC2 + PAF + TMCC3 + CRP + FAM20A + LRRN1 + C5orf38 + MENT + APOB + PEF1
dev.off()









