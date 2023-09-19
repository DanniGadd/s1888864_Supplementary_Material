############################################################################################
############################### Plotting EWAS results ######################################
############################################################################################
############################################################################################

screen

R

all <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/cis_trans/slice_neuro_final_35_pQTMs.csv")

all <- all[which(all$Association.Type %in% "TRANS"),] # 15 associations
unique(all[,9]) # 9 seqids
unique(all[,10]) # 9 proteins
unique(all[,1]) # 11 cpgs

cpgs <- all

names(cpgs)[3] <- "Gene_of_Hit"
cpgs[which(cpgs$Gene_of_Hit %in% ""),"Gene_of_Hit"] <- "Unannotated"

## Make first bed file - positions of CpGss
bed1 <- matrix(nrow = nrow(cpgs), ncol = 4)
bed1 <- as.data.frame(bed1)
names(bed1) <- c("chr", "start", "end", "value1")

## Update file naming 
names(cpgs)[2] <- "Chromosome_of_Hit"
names(cpgs)[4] <- "Position_of_Hit"
names(cpgs)[1] <- "Probe"
names(cpgs)[6] <- "b"

## Set up loop to populate bed1 file 
list = cpgs$Probe
i <- 1
for(i in 1:length(list)){ 
## get chromosome 
chr = cpgs[i, "Chromosome_of_Hit"]
chr = paste0("chr", chr)
## get start/end 
start.end <- cpgs[i, "Position_of_Hit"]
## value = beta
value <- cpgs[i, "b"]
## populate the table
bed1[i, "chr"] <- chr
bed1[i, "start"] <- start.end
bed1[i, "end"] <- start.end+5e5
bed1[i, "value1"] <- value
}

## Create bed2 file - target of interest (genes of interest)
bed2 <- matrix(nrow = nrow(cpgs), ncol = 3)
bed2 <- as.data.frame(bed2)
names(bed2) <- c("chr", "start", "end")

# Update naming for somamers 
names(cpgs)[15] <- "Chromosome_of_Somamer"
names(cpgs)[13] <- "Position_of_Somamer"
names(cpgs)[10] <- "Gene_of_Somamer"

# Get the info for somamers 
for(i in 1:length(list)){ 
## get chromosome 
chr = cpgs[i, "Chromosome_of_Somamer"]
chr = paste0("chr", chr)
## get start/end 
start.end <- cpgs[i, "Position_of_Somamer"]
## populate the table
bed2[i, "chr"] <- chr
bed2[i, "start"] <- start.end
bed2[i, "end"] <- start.end+5e5
}

## Create annotation file 
names1 <- bed1
names2 <- bed2
## give value1 column in bed1 (now names1) gene names of hit + denote its a SNP 
names1$value1 <- cpgs$Gene_of_Hit
names1$value1 <- paste(names1$value1, "(CpG)")
## give value1 column in bed2 (now names2) gene names of Somamer + denote its a Somamer
names2$value1 <- 0
names2$value1 <- cpgs$Gene_of_Somamer
names2$value1 <- paste(names2$value1, "")
## create bed file - annotation 
bed = rbind(names1, names2)

# Write out file and edit so that each CpG is only annotated to a single key gene
write.csv(bed, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/circos/circus_g_dataset_cog_joint.csv", row.names = F)

# Read it in for plotting for the report 
library(readxl)
bed <- read_excel("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/circos/circus_g_dataset_cog_joint.xlsx")
bed <- as.data.frame(bed)

library(circlize)

# Write source data
write.csv(bed, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Source_data/Fig5_source_circos.csv", row.names = F)

pdf("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Interpretation/circos/circ4_V3.pdf", width = 12, height = 12)
## initalize plot 
circos.initializeWithIdeogram()
circos.genomicTrackPlotRegion(bed, ylim = c(0, 1), track.height = 0.1, bg.border = NA)
i_track = get.cell.meta.data("track.index") # remember this empty track, we'll come back
circos.genomicTrackPlotRegion(bed, ylim = c(0, 1),
                              panel.fun = function(region, value, ...) {
                                circos.genomicText(region, value, y = 1, labels.column = 1,
                                                   facing = "clockwise", adj = c(1, 0.4),
                                                   posTransform = posTransform.text, cex = 0.8, padding =0.4)
                              }, track.height = 0.19, bg.border = NA)

tr_track = get.cell.meta.data("track.index") # position transformation track
# because `circos.genomicPosTransformLines` is implemented by
# `circos.trackPlotRegion`, it accepts `track.index` argument.
circos.genomicPosTransformLines(bed,
                                posTransform = function(region, value)
                                  posTransform.text(region, y = 1, labels = value[[1]],
                                                    cex = 0.7, padding = 0.4, track.index = tr_track),
                                direction = "inside", track.index = i_track
)

## link the trans CpGss to their target Somamer/gene 
circos.genomicLink(bed1,bed2, col = rand_color(nrow(bed1)), border= NA)
dev.off()

