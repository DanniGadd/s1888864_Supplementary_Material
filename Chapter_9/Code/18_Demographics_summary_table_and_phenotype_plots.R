###############################################################################################

### SUMMARY TABLE

###############################################################################################

### SUMMARISE PHENOTYPES 

screen

R

library(tidyverse)
library(ggplot2)

###############################################################################################

# Train GDF
d1_20k <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_Cox_proteins_20k/d1.csv") 
trainGDF <- read.csv("/Local_Scratch/Danni/GDFBNP/01_Inputs/GDF_W1W3_combined_protein_pQTL.csv")

target <- readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/GS20k/GS20k_Targets.rds")
target <- target[c(1:2)]
trainGDF <- left_join(trainGDF, target, by = "Sample_Sentrix_ID")

d1_20k <- d1_20k[which(d1_20k$Sample_Name %in% trainGDF$Sample_Name),] # 8207

## Age 

mean(d1_20k$age.x)
sd(d1_20k$age.x)
min(d1_20k$age.x, na.rm = T)
max(d1_20k$age.x, na.rm = T)
IQR(d1_20k$age.x, na.rm = T)

# > mean(d1_20k$age.x)
# [1] 50.38778
# > sd(d1_20k$age.x)
# [1] 13.32607
# > min(d1_20k$age.x, na.rm = T)
# [1] 18
# > max(d1_20k$age.x, na.rm = T)
# [1] 94.5
# > IQR(d1_20k$age.x, na.rm = T)
# [1] 17.16667

## Sex
table(d1_20k$sex.x)
#   F    M
# 4821 3386

(4821 / 8207) * 100
(3386 / 8207) * 100


# > (4821 / 8207) * 100
# [1] 58.74254
# > (3386 / 8207) * 100
# [1] 41.25746


## Protein Ns

# Mean protein measurements (untransformed)

mean(d1_20k$gdf15, na.rm = T)
sd(d1_20k$gdf15, na.rm = T)

# [1] 1059.161
# > sd(d1_20k$gdf15, na.rm = T)
# [1] 914.3448

# Min max protein measurements and IQR

min(d1_20k$gdf15, na.rm = T)
max(d1_20k$gdf15, na.rm = T)
IQR(d1_20k$gdf15, na.rm = T)

# > min(d1_20k$gdf15, na.rm = T)
# [1] 400.7
# > max(d1_20k$gdf15, na.rm = T)
# [1] 19511
# > IQR(d1_20k$gdf15, na.rm = T)
# [1] 515.2


# BMI 

table(is.na(d1_20k$bmi)) 
mean(d1_20k$bmi, na.rm = T)
sd(d1_20k$bmi, na.rm = T)
IQR(d1_20k$bmi, na.rm = T)[1]
range(d1_20k$bmi, na.rm = T)[1]
range(d1_20k$bmi, na.rm = T)[2]

# > table(is.na(d1_20k$bmi))

# FALSE  TRUE
#  8117    90
# > mean(d1_20k$bmi, na.rm = T)
# [1] 26.78224
# > sd(d1_20k$bmi, na.rm = T)
# [1] 4.848638
# > IQR(d1_20k$bmi, na.rm = T)[1]
# [1] 6.15
# > range(d1_20k$bmi, na.rm = T)[1]
# [1] 15.96
# > range(d1_20k$bmi, na.rm = T)[2]
# [1] 45.02


## Alc

table(is.na(d1_20k$units)) 
mean(d1_20k$units, na.rm = T)
sd(d1_20k$units, na.rm = T)
IQR(d1_20k$units, na.rm = T)[1]
range(d1_20k$units, na.rm = T)[1]
range(d1_20k$units, na.rm = T)[2]

# > table(is.na(d1_20k$units))

# FALSE  TRUE
#  7410   797
# > mean(d1_20k$units, na.rm = T)
# [1] 10.11957
# > sd(d1_20k$units, na.rm = T)
# [1] 10.24536
# > IQR(d1_20k$units, na.rm = T)[1]
# [1] 13
# > range(d1_20k$units, na.rm = T)[1]
# [1] 0
# > range(d1_20k$units, na.rm = T)[2]
# [1] 54


## Smoking

table(is.na(d1_20k$smokingScore)) 
mean(d1_20k$smokingScore, na.rm = T)
sd(d1_20k$smokingScore, na.rm = T)
IQR(d1_20k$smokingScore, na.rm = T)[1]
range(d1_20k$smokingScore, na.rm = T)[1]
range(d1_20k$smokingScore, na.rm = T)[2]



# > table(is.na(d1_20k$smokingScore))

# FALSE
#  8207
# > mean(d1_20k$smokingScore, na.rm = T)
# [1] 18.03746
# > sd(d1_20k$smokingScore, na.rm = T)
# [1] 31.03583
# > IQR(d1_20k$smokingScore, na.rm = T)[1]
# [1] 31.24659
# > range(d1_20k$smokingScore, na.rm = T)[1]
# [1] -41.91793
# > range(d1_20k$smokingScore, na.rm = T)[2]
# [1] 155.8305


## EA

table(is.na(d1_20k$years)) 
mean(d1_20k$years, na.rm = T)
sd(d1_20k$years, na.rm = T)
IQR(d1_20k$years, na.rm = T)[1]
range(d1_20k$years, na.rm = T)[1]
range(d1_20k$years, na.rm = T)[2]


# > table(is.na(d1_20k$years))

# FALSE  TRUE
#  7841   366
# > mean(d1_20k$years, na.rm = T)
# [1] 4.552481
# > sd(d1_20k$years, na.rm = T)
# [1] 1.602081
# > IQR(d1_20k$years, na.rm = T)[1]
# [1] 3
# > range(d1_20k$years, na.rm = T)[1]
# [1] 0
# > range(d1_20k$years, na.rm = T)[2]
# [1] 10


## SIMD


table(is.na(d1_20k$rank)) 
mean(d1_20k$rank, na.rm = T)
sd(d1_20k$rank, na.rm = T)
IQR(d1_20k$rank, na.rm = T)[1]
range(d1_20k$rank, na.rm = T)[1]
range(d1_20k$rank, na.rm = T)[2]


# > table(is.na(d1_20k$rank))

# FALSE  TRUE
#  7796   411
# > mean(d1_20k$rank, na.rm = T)
# [1] 3891.804
# > sd(d1_20k$rank, na.rm = T)
# [1] 1879.429
# > IQR(d1_20k$rank, na.rm = T)[1]
# [1] 3216.5
# > range(d1_20k$rank, na.rm = T)[1]
# [1] 3
# > range(d1_20k$rank, na.rm = T)[2]
# [1] 6505



###############################################################################################

# Train Nt-proBNP
d1_20k <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_Cox_proteins_20k/d1.csv") 
trainBNP <- read.csv("/Local_Scratch/Danni/GDFBNP/01_Inputs/BNP_W1W3_combined_protein_pQTL.csv")

target <- readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/GS20k/GS20k_Targets.rds")
target <- target[c(1:2)]
trainBNP <- left_join(trainBNP, target, by = "Sample_Sentrix_ID")

d1_20k <- d1_20k[which(d1_20k$Sample_Name %in% trainBNP$Sample_Name),] # 8002


## Age 

mean(d1_20k$age.x)
sd(d1_20k$age.x)
min(d1_20k$age.x, na.rm = T)
max(d1_20k$age.x, na.rm = T)
IQR(d1_20k$age.x, na.rm = T)


# > mean(d1_20k$age.x)
# [1] 50.37968
# > sd(d1_20k$age.x)
# [1] 13.49243
# > min(d1_20k$age.x, na.rm = T)
# [1] 18
# > max(d1_20k$age.x, na.rm = T)
# [1] 94.5
# > IQR(d1_20k$age.x, na.rm = T)
# [1] 17.41667


## Sex
table(d1_20k$sex.x)

#    F    M
# 4910 3092

(4910 / 8002) * 100
(3092 / 8002) * 100


# > (4910 / 8002) * 100
# [1] 61.35966
# > (3092 / 8002) * 100
# [1] 38.64034



## Protein Ns

# Mean protein measurements (untransformed)

mean(d1_20k$nt.probnp, na.rm = T)
sd(d1_20k$nt.probnp, na.rm = T)

# > mean(d1_20k$nt.probnp, na.rm = T)
# [1] 93.91367
# > sd(d1_20k$nt.probnp, na.rm = T)
# [1] 195.2636


# Min max protein measurements and IQR

min(d1_20k$nt.probnp, na.rm = T)
max(d1_20k$nt.probnp, na.rm = T)
IQR(d1_20k$nt.probnp, na.rm = T)

# >
# > min(d1_20k$nt.probnp, na.rm = T)
# [1] 10
# > max(d1_20k$nt.probnp, na.rm = T)
# [1] 6204
# > IQR(d1_20k$nt.probnp, na.rm = T)
# [1] 69.8575


# BMI 

table(is.na(d1_20k$bmi)) 
mean(d1_20k$bmi, na.rm = T)
sd(d1_20k$bmi, na.rm = T)
IQR(d1_20k$bmi, na.rm = T)[1]
range(d1_20k$bmi, na.rm = T)[1]
range(d1_20k$bmi, na.rm = T)[2]


# > table(is.na(d1_20k$bmi))

# FALSE  TRUE
#  7899   103
# > mean(d1_20k$bmi, na.rm = T)
# [1] 26.69895
# > sd(d1_20k$bmi, na.rm = T)
# [1] 4.852885
# > IQR(d1_20k$bmi, na.rm = T)[1]
# [1] 6.15
# > range(d1_20k$bmi, na.rm = T)[1]
# [1] 15.96
# > range(d1_20k$bmi, na.rm = T)[2]
# [1] 45.02
# >



## Alc

table(is.na(d1_20k$units)) 
mean(d1_20k$units, na.rm = T)
sd(d1_20k$units, na.rm = T)
IQR(d1_20k$units, na.rm = T)[1]
range(d1_20k$units, na.rm = T)[1]
range(d1_20k$units, na.rm = T)[2]


# > table(is.na(d1_20k$units))

# FALSE  TRUE
#  7220   782
# > mean(d1_20k$units, na.rm = T)
# [1] 9.902216
# > sd(d1_20k$units, na.rm = T)
# [1] 10.05241
# > IQR(d1_20k$units, na.rm = T)[1]
# [1] 12
# > range(d1_20k$units, na.rm = T)[1]
# [1] 0
# > range(d1_20k$units, na.rm = T)[2]
# [1] 54



## Smoking

table(is.na(d1_20k$smokingScore)) 
mean(d1_20k$smokingScore, na.rm = T)
sd(d1_20k$smokingScore, na.rm = T)
IQR(d1_20k$smokingScore, na.rm = T)[1]
range(d1_20k$smokingScore, na.rm = T)[1]
range(d1_20k$smokingScore, na.rm = T)[2]

# FALSE
#  8002
# > mean(d1_20k$smokingScore, na.rm = T)
# [1] 17.83443
# > sd(d1_20k$smokingScore, na.rm = T)
# [1] 30.94602
# > IQR(d1_20k$smokingScore, na.rm = T)[1]
# [1] 31.00302
# > range(d1_20k$smokingScore, na.rm = T)[1]
# [1] -41.91793
# > range(d1_20k$smokingScore, na.rm = T)[2]
# [1] 155.8305


## EA

table(is.na(d1_20k$years)) 
mean(d1_20k$years, na.rm = T)
sd(d1_20k$years, na.rm = T)
IQR(d1_20k$years, na.rm = T)[1]
range(d1_20k$years, na.rm = T)[1]
range(d1_20k$years, na.rm = T)[2]

# > table(is.na(d1_20k$years))

# FALSE  TRUE
#  7639   363
# > mean(d1_20k$years, na.rm = T)
# [1] 4.559759
# > sd(d1_20k$years, na.rm = T)
# [1] 1.604298
# > IQR(d1_20k$years, na.rm = T)[1]
# [1] 3
# > range(d1_20k$years, na.rm = T)[1]
# [1] 0
# > range(d1_20k$years, na.rm = T)[2]
# [1] 10

## SIMD

table(is.na(d1_20k$rank)) 
mean(d1_20k$rank, na.rm = T)
sd(d1_20k$rank, na.rm = T)
IQR(d1_20k$rank, na.rm = T)[1]
range(d1_20k$rank, na.rm = T)[1]
range(d1_20k$rank, na.rm = T)[2]


# > table(is.na(d1_20k$rank))

# FALSE  TRUE
#  7604   398
# > mean(d1_20k$rank, na.rm = T)
# [1] 3899.805
# > sd(d1_20k$rank, na.rm = T)
# [1] 1873.698
# > IQR(d1_20k$rank, na.rm = T)[1]
# [1] 3187.25
# > range(d1_20k$rank, na.rm = T)[1]
# [1] 3
# > range(d1_20k$rank, na.rm = T)[2]
# [1] 6505
# >



###############################################################################################

# Test GDF
d1_w4 <- read.csv("/Local_Scratch/Danni/GDFBNP/03_Cox/d1_W4_sample.csv")
testGDF <- d1_w4[complete.cases(d1_w4$gdf15_rnk),] 

d1_20k <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_Cox_proteins_20k/d1.csv") 
d1_20k <- d1_20k[which(d1_20k$Sample_Name %in% testGDF$Sample_Name),] # 2954


## Age 

mean(d1_20k$age.x)
sd(d1_20k$age.x)
min(d1_20k$age.x, na.rm = T)
max(d1_20k$age.x, na.rm = T)
IQR(d1_20k$age.x, na.rm = T)

# > mean(d1_20k$age.x)
# [1] 48.94507
# > sd(d1_20k$age.x)
# [1] 16.0191
# > min(d1_20k$age.x, na.rm = T)
# [1] 17.08333
# > max(d1_20k$age.x, na.rm = T)
# [1] 91.33333
# > IQR(d1_20k$age.x, na.rm = T)
# [1] 24.75


## Sex
table(d1_20k$sex.x)


#    F    M
# 1719 1235


(1719 / 2954) * 100
(1235 / 2954) * 100


# > (1719 / 2954) * 100
# [1] 58.19228
# > (1235 / 2954) * 100
# [1] 41.80772


## Protein Ns

# Mean protein measurements (untransformed)

mean(d1_20k$gdf15, na.rm = T)
sd(d1_20k$gdf15, na.rm = T)


# > mean(d1_20k$gdf15, na.rm = T)
# [1] 1099.811
# > sd(d1_20k$gdf15, na.rm = T)
# [1] 923.1262


# Min max protein measurements and IQR

min(d1_20k$gdf15, na.rm = T)
max(d1_20k$gdf15, na.rm = T)
IQR(d1_20k$gdf15, na.rm = T)

# > min(d1_20k$gdf15, na.rm = T)
# [1] 400.7
# > max(d1_20k$gdf15, na.rm = T)
# [1] 19216
# > IQR(d1_20k$gdf15, na.rm = T)
# [1] 594.975


# BMI 

table(is.na(d1_20k$bmi)) 
mean(d1_20k$bmi, na.rm = T)
sd(d1_20k$bmi, na.rm = T)
IQR(d1_20k$bmi, na.rm = T)[1]
range(d1_20k$bmi, na.rm = T)[1]
range(d1_20k$bmi, na.rm = T)[2]


# > table(is.na(d1_20k$bmi))

# FALSE  TRUE
#  2920    34
# > mean(d1_20k$bmi, na.rm = T)
# [1] 26.97033
# > sd(d1_20k$bmi, na.rm = T)
# [1] 5.082325
# > IQR(d1_20k$bmi, na.rm = T)[1]
# [1] 6.535
# > range(d1_20k$bmi, na.rm = T)[1]
# [1] 13.16
# > range(d1_20k$bmi, na.rm = T)[2]
# [1] 44.75




## Alc

table(is.na(d1_20k$units)) 
mean(d1_20k$units, na.rm = T)
sd(d1_20k$units, na.rm = T)
IQR(d1_20k$units, na.rm = T)[1]
range(d1_20k$units, na.rm = T)[1]
range(d1_20k$units, na.rm = T)[2]


# > table(is.na(d1_20k$units))

# FALSE  TRUE
#  2604   350
# > mean(d1_20k$units, na.rm = T)
# [1] 9.722734
# > sd(d1_20k$units, na.rm = T)
# [1] 10.24269
# > IQR(d1_20k$units, na.rm = T)[1]
# [1] 13
# > range(d1_20k$units, na.rm = T)[1]
# [1] 0
# > range(d1_20k$units, na.rm = T)[2]
# [1] 54



## Smoking

table(is.na(d1_20k$smokingScore)) 
mean(d1_20k$smokingScore, na.rm = T)
sd(d1_20k$smokingScore, na.rm = T)
IQR(d1_20k$smokingScore, na.rm = T)[1]
range(d1_20k$smokingScore, na.rm = T)[1]
range(d1_20k$smokingScore, na.rm = T)[2]

# > table(is.na(d1_20k$smokingScore))

# FALSE
#  2954
# > mean(d1_20k$smokingScore, na.rm = T)
# [1] 13.37589
# > sd(d1_20k$smokingScore, na.rm = T)
# [1] 29.49829
# > IQR(d1_20k$smokingScore, na.rm = T)[1]
# [1] 31.04293
# > range(d1_20k$smokingScore, na.rm = T)[1]
# [1] -38.99736
# > range(d1_20k$smokingScore, na.rm = T)[2]
# [1] 136.8031


## EA

table(is.na(d1_20k$years)) 
mean(d1_20k$years, na.rm = T)
sd(d1_20k$years, na.rm = T)
IQR(d1_20k$years, na.rm = T)[1]
range(d1_20k$years, na.rm = T)[1]
range(d1_20k$years, na.rm = T)[2]


# > table(is.na(d1_20k$years))

# FALSE  TRUE
#  2716   238
# > mean(d1_20k$years, na.rm = T)
# [1] 4.447349
# > sd(d1_20k$years, na.rm = T)
# [1] 1.535278
# > IQR(d1_20k$years, na.rm = T)[1]
# [1] 3
# > range(d1_20k$years, na.rm = T)[1]
# [1] 0
# > range(d1_20k$years, na.rm = T)[2]
# [1] 10


## SIMD

table(is.na(d1_20k$rank)) 
mean(d1_20k$rank, na.rm = T)
sd(d1_20k$rank, na.rm = T)
IQR(d1_20k$rank, na.rm = T)[1]
range(d1_20k$rank, na.rm = T)[1]
range(d1_20k$rank, na.rm = T)[2]

#  table(is.na(d1_20k$rank))

# FALSE  TRUE
#  2750   204
# > mean(d1_20k$rank, na.rm = T)
# [1] 3763.106
# > sd(d1_20k$rank, na.rm = T)
# [1] 1834.621
# > IQR(d1_20k$rank, na.rm = T)[1]
# [1] 3066.75
# > range(d1_20k$rank, na.rm = T)[1]
# [1] 15
# > range(d1_20k$rank, na.rm = T)[2]
# [1] 6503



###############################################################################################

# Test Nt-proBNP
d1_w4 <- read.csv("/Local_Scratch/Danni/GDFBNP/03_Cox/d1_W4_sample.csv")
testBNP <- d1_w4[complete.cases(d1_w4$nt.probnp_rnk),]

d1_20k <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_Cox_proteins_20k/d1.csv") 
d1_20k <- d1_20k[which(d1_20k$Sample_Name %in% testBNP$Sample_Name),] # 2808


## Age 

mean(d1_20k$age.x)
sd(d1_20k$age.x)
min(d1_20k$age.x, na.rm = T)
max(d1_20k$age.x, na.rm = T)
IQR(d1_20k$age.x, na.rm = T)


# > mean(d1_20k$age.x)
# [1] 48.96026
# > sd(d1_20k$age.x)
# [1] 16.24782
# > min(d1_20k$age.x, na.rm = T)
# [1] 18
# > max(d1_20k$age.x, na.rm = T)
# [1] 91.33333
# > IQR(d1_20k$age.x, na.rm = T)
# [1] 25.10417

## Sex
table(d1_20k$sex.x)


#  F    M
# 1708 1100

(1708 / 2808) * 100
(1100 / 2808) * 100



# > (1708 / 2808) * 100
# [1] 60.82621
# > (1100 / 2808) * 100
# [1] 39.17379



## Protein Ns

# Mean protein measurements (untransformed)

mean(d1_20k$nt.probnp, na.rm = T)
sd(d1_20k$nt.probnp, na.rm = T)

# > mean(d1_20k$nt.probnp, na.rm = T)
# [1] 112.9348
# > sd(d1_20k$nt.probnp, na.rm = T)
# [1] 234.8663
# >

# Min max protein measurements and IQR

min(d1_20k$nt.probnp, na.rm = T)
max(d1_20k$nt.probnp, na.rm = T)
IQR(d1_20k$nt.probnp, na.rm = T)

# > min(d1_20k$nt.probnp, na.rm = T)
# [1] 10.03
# > max(d1_20k$nt.probnp, na.rm = T)
# [1] 4647
# > IQR(d1_20k$nt.probnp, na.rm = T)
# [1] 78.8575


# BMI 

table(is.na(d1_20k$bmi)) 
mean(d1_20k$bmi, na.rm = T)
sd(d1_20k$bmi, na.rm = T)
IQR(d1_20k$bmi, na.rm = T)[1]
range(d1_20k$bmi, na.rm = T)[1]
range(d1_20k$bmi, na.rm = T)[2]


# > table(is.na(d1_20k$bmi))

# FALSE  TRUE
#  2771    37
# > mean(d1_20k$bmi, na.rm = T)
# [1] 26.85827
# > sd(d1_20k$bmi, na.rm = T)
# [1] 5.101292
# > IQR(d1_20k$bmi, na.rm = T)[1]
# [1] 6.545
# > range(d1_20k$bmi, na.rm = T)[1]
# [1] 13.16
# > range(d1_20k$bmi, na.rm = T)[2]
# [1] 44.75




## Alc

table(is.na(d1_20k$units)) 
mean(d1_20k$units, na.rm = T)
sd(d1_20k$units, na.rm = T)
IQR(d1_20k$units, na.rm = T)[1]
range(d1_20k$units, na.rm = T)[1]
range(d1_20k$units, na.rm = T)[2]



# > table(is.na(d1_20k$units))

# FALSE  TRUE
#  2470   338
# > mean(d1_20k$units, na.rm = T)
# [1] 9.651012
# > sd(d1_20k$units, na.rm = T)
# [1] 10.10678
# > IQR(d1_20k$units, na.rm = T)[1]
# [1] 13
# > range(d1_20k$units, na.rm = T)[1]
# [1] 0
# > range(d1_20k$units, na.rm = T)[2]
# [1] 54
# >


## Smoking

table(is.na(d1_20k$smokingScore)) 
mean(d1_20k$smokingScore, na.rm = T)
sd(d1_20k$smokingScore, na.rm = T)
IQR(d1_20k$smokingScore, na.rm = T)[1]
range(d1_20k$smokingScore, na.rm = T)[1]
range(d1_20k$smokingScore, na.rm = T)[2]

# table(is.na(d1_20k$smokingScore))

# FALSE
#  2808
# > mean(d1_20k$smokingScore, na.rm = T)
# [1] 13.33283
# > sd(d1_20k$smokingScore, na.rm = T)
# [1] 29.42344
# > IQR(d1_20k$smokingScore, na.rm = T)[1]
# [1] 30.97186
# > range(d1_20k$smokingScore, na.rm = T)[1]
# [1] -38.99736
# > range(d1_20k$smokingScore, na.rm = T)[2]
# [1] 136.8031



## EA

table(is.na(d1_20k$years)) 
mean(d1_20k$years, na.rm = T)
sd(d1_20k$years, na.rm = T)
IQR(d1_20k$years, na.rm = T)[1]
range(d1_20k$years, na.rm = T)[1]
range(d1_20k$years, na.rm = T)[2]


# > table(is.na(d1_20k$years))

# FALSE  TRUE
#  2582   226
# > mean(d1_20k$years, na.rm = T)
# [1] 4.465143
# > sd(d1_20k$years, na.rm = T)
# [1] 1.541668
# > IQR(d1_20k$years, na.rm = T)[1]
# [1] 3
# > range(d1_20k$years, na.rm = T)[1]
# [1] 0
# > range(d1_20k$years, na.rm = T)[2]
# [1] 10


## SIMD

table(is.na(d1_20k$rank)) 
mean(d1_20k$rank, na.rm = T)
sd(d1_20k$rank, na.rm = T)
IQR(d1_20k$rank, na.rm = T)[1]
range(d1_20k$rank, na.rm = T)[1]
range(d1_20k$rank, na.rm = T)[2]

# > table(is.na(d1_20k$rank))

# FALSE  TRUE
#  2609   199
# > mean(d1_20k$rank, na.rm = T)
# [1] 3779.626
# > sd(d1_20k$rank, na.rm = T)
# [1] 1828.358
# > IQR(d1_20k$rank, na.rm = T)[1]
# [1] 3044
# > range(d1_20k$rank, na.rm = T)[1]
# [1] 15
# > range(d1_20k$rank, na.rm = T)[2]
# [1] 6503




###############################################################################################

### GROUP ONE - full sample with DNAm available and either BNP or GDF levels 

d1_20k <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_Cox_proteins_20k/d1.csv") 

table(is.na(d1_20k$smokingScore))
# 18413 with DNAm 

## Age 
mean(d1_20k$age.x)
# [1] 47.51012

sd(d1_20k$age.x)
# [1] 14.92897

min(d1_20k$age.x, na.rm = T)
max(d1_20k$age.x, na.rm = T)
IQR(d1_20k$age.x, na.rm = T)


# > min(d1_20k$age.x, na.rm = T)
# [1] 17.08333
# > max(d1_20k$age.x, na.rm = T)
# [1] 98.5
# > IQR(d1_20k$age.x, na.rm = T)
# [1] 22.41667


## Sex
table(d1_20k$sex.x)
#     F     M
# 10833  7580

(10833 / 18413) * 100

# [1] 58.83343

(7580 / 18413) * 100

# [1] 41.16657

## Protein Ns

table(is.na(d1$gdf15_rnk))
table(is.na(d1$nt.probnp_rnk))

# FALSE  TRUE
# 17489   924

# FALSE  TRUE
# 16963  1450

# Mean protein measurements (untransformed)

mean(d1_20k$gdf15, na.rm = T)
sd(d1_20k$gdf15, na.rm = T)

mean(d1_20k$nt.probnp, na.rm = T)
sd(d1_20k$nt.probnp, na.rm = T)

# > mean(d1_20k$gdf15, na.rm = T)
# [1] 1038.663
# > sd(d1_20k$gdf15, na.rm = T)
# [1] 928.0372
# >
# > mean(d1_20k$nt.probnp, na.rm = T)
# [1] 94.58287
# > sd(d1_20k$nt.probnp, na.rm = T)
# [1] 211.2274

# Min max protein measurements and IQR

min(d1_20k$gdf15, na.rm = T)
max(d1_20k$gdf15, na.rm = T)
IQR(d1_20k$gdf15, na.rm = T)


# > min(d1_20k$gdf15, na.rm = T)
# [1] 400.1
# > max(d1_20k$gdf15, na.rm = T)
# [1] 19511
# > IQR(d1_20k$gdf15, na.rm = T)
# [1] 516.3

min(d1_20k$nt.probnp, na.rm = T)
max(d1_20k$nt.probnp, na.rm = T)
IQR(d1_20k$nt.probnp, na.rm = T)

# > min(d1_20k$nt.probnp, na.rm = T)
# [1] 10
# > max(d1_20k$nt.probnp, na.rm = T)
# [1] 12967
# > IQR(d1_20k$nt.probnp, na.rm = T)
# [1] 69.2


# BMI 

table(is.na(d1_20k$bmi)) 
mean(d1_20k$bmi, na.rm = T)
sd(d1_20k$bmi, na.rm = T)
IQR(d1_20k$bmi, na.rm = T)[1]
range(d1_20k$bmi, na.rm = T)[1]
range(d1_20k$bmi, na.rm = T)[2]

# > table(is.na(d1_20k$bmi))

# FALSE  TRUE
# 18181   232
# > mean(d1_20k$bmi, na.rm = T)
# [1] 26.52022
# > sd(d1_20k$bmi, na.rm = T)
# [1] 4.853012
# > IQR(d1_20k$bmi, na.rm = T)[1]
# [1] 6.15
# > range(d1_20k$bmi, na.rm = T)[1]
# [1] 10.49
# > range(d1_20k$bmi, na.rm = T)[2]
# [1] 45.04

## Alc

table(is.na(d1_20k$units)) 
mean(d1_20k$units, na.rm = T)
sd(d1_20k$units, na.rm = T)
IQR(d1_20k$units, na.rm = T)[1]
range(d1_20k$units, na.rm = T)[1]
range(d1_20k$units, na.rm = T)[2]


# FALSE  TRUE
# 16522  1891
# > mean(d1_20k$units, na.rm = T)
# [1] 10.18339
# > sd(d1_20k$units, na.rm = T)
# [1] 10.36741
# > IQR(d1_20k$units, na.rm = T)[1]
# [1] 13
# > range(d1_20k$units, na.rm = T)[1]
# [1] 0
# > range(d1_20k$units, na.rm = T)[2]
# [1] 54
# >


## Smoking

table(is.na(d1_20k$smokingScore)) 
mean(d1_20k$smokingScore, na.rm = T)
sd(d1_20k$smokingScore, na.rm = T)
IQR(d1_20k$smokingScore, na.rm = T)[1]
range(d1_20k$smokingScore, na.rm = T)[1]
range(d1_20k$smokingScore, na.rm = T)[2]


# FALSE
# 18413
# > mean(d1_20k$smokingScore, na.rm = T)
# [1] 14.33258
# > sd(d1_20k$smokingScore, na.rm = T)
# [1] 29.69988
# > IQR(d1_20k$smokingScore, na.rm = T)[1]
# [1] 30.00076
# > range(d1_20k$smokingScore, na.rm = T)[1]
# [1] -43.44238
# > range(d1_20k$smokingScore, na.rm = T)[2]
# [1] 155.8305
# >

## EA

table(is.na(d1_20k$years)) 
mean(d1_20k$years, na.rm = T)
sd(d1_20k$years, na.rm = T)
IQR(d1_20k$years, na.rm = T)[1]
range(d1_20k$years, na.rm = T)[1]
range(d1_20k$years, na.rm = T)[2]


# FALSE  TRUE
# 17389  1024
# > mean(d1_20k$years, na.rm = T)
# [1] 4.625395
# > sd(d1_20k$years, na.rm = T)
# [1] 1.586582
# > IQR(d1_20k$years, na.rm = T)[1]
# [1] 3
# > range(d1_20k$years, na.rm = T)[1]
# [1] 0
# > range(d1_20k$years, na.rm = T)[2]
# [1] 10
# >

## SIMD


table(is.na(d1_20k$rank)) 
mean(d1_20k$rank, na.rm = T)
sd(d1_20k$rank, na.rm = T)
IQR(d1_20k$rank, na.rm = T)[1]
range(d1_20k$rank, na.rm = T)[1]
range(d1_20k$rank, na.rm = T)[2]


# FALSE  TRUE
# 17287  1126
# > mean(d1_20k$rank, na.rm = T)
# [1] 3890.199
# > sd(d1_20k$rank, na.rm = T)
# [1] 1848.004
# > IQR(d1_20k$rank, na.rm = T)[1]
# [1] 3114.5
# > range(d1_20k$rank, na.rm = T)[1]
# [1] 1
# > range(d1_20k$rank, na.rm = T)[2]
# [1] 6505
# >

# ####################################################################

# ## SUBSET TO testing sample and repeat 

# d1_w4 <- read.csv("/Local_Scratch/Danni/GDFBNP/03_Cox/d1_W4_sample.csv")

# table(is.na(d1_w4$gdf15_rnk))
# table(is.na(d1_w4$nt.probnp_rnk))


# # > table(is.na(d1_w4$gdf15_rnk))

# # FALSE  TRUE
# #  2954 21138
# # >
# # > table(is.na(d1_w4$nt.probnp_rnk))

# # FALSE  TRUE
# #  2808 21284

# d1_w4 <- d1_w4[which(d1_w4$Set %in% c("wave4")),] # 8876

# # Subset to those that have either GDF or BNP 


# d1_20k <- d1_w4

# ## Age 
# mean(d1_20k$age.x)
# sd(d1_20k$age.x)
# min(d1_20k$age.x, na.rm = T)
# max(d1_20k$age.x, na.rm = T)
# IQR(d1_20k$age.x, na.rm = T)


# # > mean(d1_20k$age.x)
# # [1] 44.99198
# # > sd(d1_20k$age.x)
# # [1] 15.76349
# # > min(d1_20k$age.x, na.rm = T)
# # [1] 17.08333
# # > max(d1_20k$age.x, na.rm = T)
# # [1] 98.5
# # > IQR(d1_20k$age.x, na.rm = T)
# # [1] 25.33333


# ## Sex
# table(d1_20k$sex.x)
# #    F    M
# # 5196 3680


# (5196 / 8876) * 100

# # [1] 58.53988

# (3680 / 8876) * 100

# # [1] 41.46012



# ## Protein Ns

# table(is.na(d1_20k$gdf15_rnk))
# table(is.na(d1_20k$nt.probnp_rnk))

# # > table(is.na(d1_20k$gdf15_rnk))

# # FALSE  TRUE
# #  8369   507
# # > table(is.na(d1_20k$nt.probnp_rnk))

# # FALSE  TRUE
# #  8123   753


# # Mean protein measurements (untransformed)

# mean(d1_20k$gdf15, na.rm = T)
# sd(d1_20k$gdf15, na.rm = T)

# mean(d1_20k$nt.probnp, na.rm = T)
# sd(d1_20k$nt.probnp, na.rm = T)

# # > mean(d1_20k$gdf15, na.rm = T)
# # [1] 1014.162
# # > sd(d1_20k$gdf15, na.rm = T)
# # [1] 956.4494
# # >
# # > mean(d1_20k$nt.probnp, na.rm = T)
# # [1] 93.93442
# # > sd(d1_20k$nt.probnp, na.rm = T)
# # [1] 223.5503

# # Min max protein measurements and IQR

# min(d1_20k$gdf15, na.rm = T)
# max(d1_20k$gdf15, na.rm = T)
# IQR(d1_20k$gdf15, na.rm = T)


# # > min(d1_20k$gdf15, na.rm = T)
# # [1] 400.1
# # > max(d1_20k$gdf15, na.rm = T)
# # [1] 19216
# # > IQR(d1_20k$gdf15, na.rm = T)
# # [1] 503.4


# min(d1_20k$nt.probnp, na.rm = T)
# max(d1_20k$nt.probnp, na.rm = T)
# IQR(d1_20k$nt.probnp, na.rm = T)

# # > min(d1_20k$nt.probnp, na.rm = T)
# # [1] 10
# # > max(d1_20k$nt.probnp, na.rm = T)
# # [1] 12967
# # > IQR(d1_20k$nt.probnp, na.rm = T)
# # [1] 67.51


# # BMI 

# table(is.na(d1_20k$bmi)) 
# mean(d1_20k$bmi, na.rm = T)
# sd(d1_20k$bmi, na.rm = T)
# IQR(d1_20k$bmi, na.rm = T)[1]
# range(d1_20k$bmi, na.rm = T)[1]
# range(d1_20k$bmi, na.rm = T)[2]

# # FALSE  TRUE
# #  8762   114
# # > mean(d1_20k$bmi, na.rm = T)
# # [1] 26.27654
# # > sd(d1_20k$bmi, na.rm = T)
# # [1] 4.841904
# # > IQR(d1_20k$bmi, na.rm = T)[1]
# # [1] 6.1475
# # > range(d1_20k$bmi, na.rm = T)[1]
# # [1] 10.49
# # > range(d1_20k$bmi, na.rm = T)[2]
# # [1] 45.04
# # >


# ## Alc

# table(is.na(d1_20k$units)) 
# mean(d1_20k$units, na.rm = T)
# sd(d1_20k$units, na.rm = T)
# IQR(d1_20k$units, na.rm = T)[1]
# range(d1_20k$units, na.rm = T)[1]
# range(d1_20k$units, na.rm = T)[2]


# # FALSE  TRUE
# #  7916   960
# # > mean(d1_20k$units, na.rm = T)
# # [1] 10.34058
# # > sd(d1_20k$units, na.rm = T)
# # [1] 10.51429
# # > IQR(d1_20k$units, na.rm = T)[1]
# # [1] 14
# # > range(d1_20k$units, na.rm = T)[1]
# # [1] 0
# # > range(d1_20k$units, na.rm = T)[2]
# # [1] 54
# # >



# ## Smoking

# table(is.na(d1_20k$smokingScore)) 
# mean(d1_20k$smokingScore, na.rm = T)
# sd(d1_20k$smokingScore, na.rm = T)
# IQR(d1_20k$smokingScore, na.rm = T)[1]
# range(d1_20k$smokingScore, na.rm = T)[1]
# range(d1_20k$smokingScore, na.rm = T)[2]


# # FALSE
# #  8876
# # > mean(d1_20k$smokingScore, na.rm = T)
# # [1] 10.70406
# # > sd(d1_20k$smokingScore, na.rm = T)
# # [1] 28.05125
# # > IQR(d1_20k$smokingScore, na.rm = T)[1]
# # [1] 26.76252
# # > range(d1_20k$smokingScore, na.rm = T)[1]
# # [1] -38.99736
# # > range(d1_20k$smokingScore, na.rm = T)[2]
# # [1] 138.307


# ## EA

# table(is.na(d1_20k$years)) 
# mean(d1_20k$years, na.rm = T)
# sd(d1_20k$years, na.rm = T)
# IQR(d1_20k$years, na.rm = T)[1]
# range(d1_20k$years, na.rm = T)[1]
# range(d1_20k$years, na.rm = T)[2]


# # FALSE  TRUE
# #  8282   594
# # > mean(d1_20k$years, na.rm = T)
# # [1] 4.704057
# # > sd(d1_20k$years, na.rm = T)
# # [1] 1.572606
# # > IQR(d1_20k$years, na.rm = T)[1]
# # [1] 3
# # > range(d1_20k$years, na.rm = T)[1]
# # [1] 0
# # > range(d1_20k$years, na.rm = T)[2]
# # [1] 10


# ## SIMD


# table(is.na(d1_20k$rank)) 
# mean(d1_20k$rank, na.rm = T)
# sd(d1_20k$rank, na.rm = T)
# IQR(d1_20k$rank, na.rm = T)[1]
# range(d1_20k$rank, na.rm = T)[1]
# range(d1_20k$rank, na.rm = T)[2]


# # FALSE  TRUE
# #  8244   632
# # > mean(d1_20k$rank, na.rm = T)
# # [1] 3895.499
# # > sd(d1_20k$rank, na.rm = T)
# # [1] 1822.693
# # > IQR(d1_20k$rank, na.rm = T)[1]
# # [1] 3019.5
# # > range(d1_20k$rank, na.rm = T)[1]
# # [1] 1
# # > range(d1_20k$rank, na.rm = T)[2]
# # [1] 6503
# # >



# ####################################################################

# ## SUBSET TO training sample and repeat 

# # Training samples after exclusions
# GDF_y <- read.csv("/Local_Scratch/Danni/GDFBNP/01_Inputs/GDF_W1W3_combined_protein_pQTL.csv")
# GDF_y <- GDF_y[1:2]
# BNP_y <- read.csv("/Local_Scratch/Danni/GDFBNP/01_Inputs/BNP_W1W3_combined_protein_pQTL.csv")
# BNP_y <- BNP_y[c(1,3)]

# d1_20k <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_Cox_proteins_20k/d1.csv") 
# d1_20k <- d1_20k[-c(139,141)]
# joint <- merge(GDF_y, BNP_y, all = T)

# target <- readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/GS20k/GS20k_Targets.rds")
# joint <- left_join(target, joint, by = "Sample_Sentrix_ID")

# d1_joint <- left_join(d1_20k, target, by = "Sample_Name")

# d1_joint <- left_join(joint, d1_joint, by = "Sample_Name")
# table(is.na(d1_joint$gdf15_rnk)) 
# table(is.na(d1_joint$nt.probnp_rnk)) 

# # > table(is.na(d1_joint$gdf15)) # 9120

# # FALSE  TRUE
# #  8207 10206
# # > table(is.na(d1_joint$nt.probnp)) # 8840

# # FALSE  TRUE
# #  8002 10411

# d1_joint <- d1_joint[which(d1_joint$Set.y %in% c("wave3", "wave1")),] # 9,537



# # Check to make sure phenotypic data is complete across the 9019, but proteins are subset 
# table(is.na(d1_joint$age.x))

# d1_20k <- d1_joint

# ## Age 
# mean(d1_20k$age.x, na.rm = T)
# sd(d1_20k$age.x, na.rm = T)
# min(d1_20k$age.x, na.rm = T)
# max(d1_20k$age.x, na.rm = T)
# IQR(d1_20k$age.x, na.rm = T)



# # > mean(d1_20k$age.x, na.rm = T)
# # [1] 49.85373
# # > sd(d1_20k$age.x, na.rm = T)
# # [1] 13.69907
# # > min(d1_20k$age.x, na.rm = T)
# # [1] 18
# # > max(d1_20k$age.x, na.rm = T)
# # [1] 94.5
# # > IQR(d1_20k$age.x, na.rm = T)
# # [1] 18.08333



# ## Sex
# table(d1_20k$sex.x)
# #    F    M
# # 5637 3900

# (5637 / 9537) * 100

# (3900 / 9537) * 100


# # > (5637 / 9537) * 100
# # [1] 59.10664
# # > (3900 / 9537) * 100
# # [1] 40.89336



# ## Protein Ns

# table(is.na(d1_20k$gdf15_rnk))
# table(is.na(d1_20k$nt.probnp_rnk))

# #  table(is.na(d1_20k$gdf15))

# # FALSE  TRUE
# #  8207  1330
# # > table(is.na(d1_20k$nt.probnp))

# # FALSE  TRUE
# #  8002  1535


# # Mean protein measurements 

# sub1 <- d1_20k[complete.cases(d1_20k$gdf15_rnk),]

# sub2 <- d1_20k[complete.cases(d1_20k$nt.probnp_rnk),]

# mean(d1_20k$gdf15, na.rm = T)
# sd(d1_20k$gdf15, na.rm = T)

# mean(d1_20k$nt.probnp, na.rm = T)
# sd(d1_20k$nt.probnp, na.rm = T)

# # > mean(d1_20k$gdf15, na.rm = T)
# # [1] 1061.145
# # > sd(d1_20k$gdf15, na.rm = T)
# # [1] 900.6431
# # >
# # > mean(d1_20k$nt.probnp, na.rm = T)
# # [1] 95.17872
# # > sd(d1_20k$nt.probnp, na.rm = T)
# # [1] 199.2442



# # Min max protein measurements and IQR

# min(d1_20k$gdf15, na.rm = T)
# max(d1_20k$gdf15, na.rm = T)
# IQR(d1_20k$gdf15, na.rm = T)



# # > min(d1_20k$gdf15, na.rm = T)
# # [1] 400.3
# # > max(d1_20k$gdf15, na.rm = T)
# # [1] 19511
# # > IQR(d1_20k$gdf15, na.rm = T)
# # [1] 524.15



# min(d1_20k$nt.probnp, na.rm = T)
# max(d1_20k$nt.probnp, na.rm = T)
# IQR(d1_20k$nt.probnp, na.rm = T)

# # > min(d1_20k$nt.probnp, na.rm = T)
# # [1] 10
# # > max(d1_20k$nt.probnp, na.rm = T)
# # [1] 6204
# # > IQR(d1_20k$nt.probnp, na.rm = T)
# # [1] 70.96



# # BMI 

# table(is.na(d1_20k$bmi)) 
# mean(d1_20k$bmi, na.rm = T)
# sd(d1_20k$bmi, na.rm = T)
# IQR(d1_20k$bmi, na.rm = T)[1]
# range(d1_20k$bmi, na.rm = T)[1]
# range(d1_20k$bmi, na.rm = T)[2]



# # FALSE  TRUE
# #  9419   118

# # > mean(d1_20k$bmi, na.rm = T)
# # [1] 26.7469
# # > sd(d1_20k$bmi, na.rm = T)
# # [1] 4.852604
# # > IQR(d1_20k$bmi, na.rm = T)[1]
# # [1] 6.155
# # > range(d1_20k$bmi, na.rm = T)[1]
# # [1] 14.78
# # > range(d1_20k$bmi, na.rm = T)[2]
# # [1] 45.02



# ## Alc

# table(is.na(d1_20k$units)) 
# mean(d1_20k$units, na.rm = T)
# sd(d1_20k$units, na.rm = T)
# IQR(d1_20k$units, na.rm = T)[1]
# range(d1_20k$units, na.rm = T)[1]
# range(d1_20k$units, na.rm = T)[2]




# # FALSE  TRUE
# #  8606   931

# # > mean(d1_20k$units, na.rm = T)
# # [1] 10.03881
# # > sd(d1_20k$units, na.rm = T)
# # [1] 10.22893
# # > IQR(d1_20k$units, na.rm = T)[1]
# # [1] 12
# # > range(d1_20k$units, na.rm = T)[1]
# # [1] 0
# # > range(d1_20k$units, na.rm = T)[2]
# # [1] 54




# ## Smoking

# table(is.na(d1_20k$smokingScore)) 
# mean(d1_20k$smokingScore, na.rm = T)
# sd(d1_20k$smokingScore, na.rm = T)
# IQR(d1_20k$smokingScore, na.rm = T)[1]
# range(d1_20k$smokingScore, na.rm = T)[1]
# range(d1_20k$smokingScore, na.rm = T)[2]





# # > table(is.na(d1_20k$smokingScore))

# # FALSE
# #  9537
# # > mean(d1_20k$smokingScore, na.rm = T)
# # [1] 17.70961
# # > sd(d1_20k$smokingScore, na.rm = T)
# # [1] 30.7754
# # > IQR(d1_20k$smokingScore, na.rm = T)[1]
# # [1] 30.85444
# # > range(d1_20k$smokingScore, na.rm = T)[1]
# # [1] -43.44238
# # > range(d1_20k$smokingScore, na.rm = T)[2]
# # [1] 155.8305


# # # >

# ## EA

# table(is.na(d1_20k$years)) 
# mean(d1_20k$years, na.rm = T)
# sd(d1_20k$years, na.rm = T)
# IQR(d1_20k$years, na.rm = T)[1]
# range(d1_20k$years, na.rm = T)[1]
# range(d1_20k$years, na.rm = T)[2]



# # > table(is.na(d1_20k$years))

# # FALSE  TRUE
# #  9107   430
# # > mean(d1_20k$years, na.rm = T)
# # [1] 4.55386
# # > sd(d1_20k$years, na.rm = T)
# # [1] 1.59591
# # > IQR(d1_20k$years, na.rm = T)[1]
# # [1] 3
# # > range(d1_20k$years, na.rm = T)[1]
# # [1] 0
# # > range(d1_20k$years, na.rm = T)[2]
# # [1] 10



# ## SIMD


# table(is.na(d1_20k$rank)) 
# mean(d1_20k$rank, na.rm = T)
# sd(d1_20k$rank, na.rm = T)
# IQR(d1_20k$rank, na.rm = T)[1]
# range(d1_20k$rank, na.rm = T)[1]
# range(d1_20k$rank, na.rm = T)[2]


# # > table(is.na(d1_20k$rank))

# # FALSE  TRUE
# #  9043   494
# # > mean(d1_20k$rank, na.rm = T)
# # [1] 3885.368
# # > sd(d1_20k$rank, na.rm = T)
# # [1] 1870.869
# # > IQR(d1_20k$rank, na.rm = T)[1]
# # [1] 3181.5
# # > range(d1_20k$rank, na.rm = T)[1]
# # [1] 3
# # > range(d1_20k$rank, na.rm = T)[2]
# # [1] 6505



# #########################################################################################

# # ### PLOT PHENOTYPES

# # screen

# # R

# # library(tidyverse)
# # library(ggplot2)

# # ## Cognitive scores 

# # # Read in the processed cognitive data from the script daniel shared with me - composite gf and g scores and other scores with outliers > 3.5 sd from mean removed 
# # comp <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS/01_Phenotype_collation/prot_file_220221.csv", check.names = F)

# # # Make plots for cognitive scores 

# # # Digit symbol 
# # d <- comp[,c("GS_id", "digit_symbol")]
# # d <- na.omit(d)
# # d$digit_symbol <- as.numeric(d$digit_symbol)
# # mean_value <- mean(d$digit_symbol)
# # bw <- 2 * IQR(d$digit_symbol) / length(d$digit_symbol)^(1/3)

# # plot1 <- ggplot(d, aes(x=digit_symbol)) + 
# #   geom_histogram(color="black", fill="azure3", binwidth = bw) + geom_vline(xintercept = mean_value,      
# #              col = "dodgerblue4",
# #              lwd = 2) + xlab("Processing Speed") + ylab("Count") +
# #  theme_classic() + theme(
# #     axis.text = element_text(size = 20), 
# #     axis.title = element_text(size = 20))

# # pdf("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS/01_Phenotype_collation/Plots/digit.pdf")
# # plot1
# # dev.off()

# # # Logical memory 
# # d <- comp[,c("GS_id", "LM")]
# # d <- na.omit(d)
# # d$LM <- as.numeric(d$LM)
# # mean_value <- mean(d$LM)
# # bw <- 2 * IQR(d$LM) / length(d$LM)^(1/3)

# # plot2 <- ggplot(d, aes(x=LM)) + 
# #   geom_histogram(color="black", fill="azure3", binwidth = bw) + geom_vline(xintercept = mean_value,      
# #              col = "dodgerblue4",
# #              lwd = 2) + xlab("Logical Memory") + ylab("Count") +
# #  theme_classic() + theme(
# #     axis.text = element_text(size = 20), 
# #     axis.title = element_text(size = 20))

# # pdf("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS/01_Phenotype_collation/Plots/LM.pdf")
# # plot2
# # dev.off()

# # # Verbal total (verbal reasoning)
# # d <- comp[,c("GS_id", "verbal_total")]
# # d <- na.omit(d)
# # d$verbal_total <- as.numeric(d$verbal_total)
# # mean_value <- mean(d$verbal_total)
# # bw <- 2 * IQR(d$verbal_total) / length(d$verbal_total)^(1/3)

# # plot3 <- ggplot(d, aes(x=verbal_total)) + 
# #   geom_histogram(color="black", fill="azure3", binwidth = bw) + geom_vline(xintercept = mean_value,      
# #              col = "dodgerblue4",
# #              lwd = 2) + xlab("Verbal Reasoning") + ylab("Count") +
# #  theme_classic() + theme(
# #     axis.text = element_text(size = 20), 
# #     axis.title = element_text(size = 20))

# # pdf("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS/01_Phenotype_collation/Plots/verbal_reasoning.pdf")
# # plot3
# # dev.off()


# # # Matrix reasoning (Non-Verbal reasoning)
# # d <- comp[,c("GS_id", "mr_correct")]
# # d <- na.omit(d)
# # d$verbal_total <- as.numeric(d$mr_correct)
# # mean_value <- mean(d$mr_correct)
# # bw <- 2 * IQR(d$mr_correct) / length(d$mr_correct)^(1/3)

# # plot4 <- ggplot(d, aes(x=mr_correct)) + 
# #   geom_histogram(color="black", fill="azure3", binwidth = bw) + geom_vline(xintercept = mean_value,      
# #              col = "dodgerblue4",
# #              lwd = 2) + xlab("Non-Verbal Reasoning") + ylab("Count") +
# #  theme_classic() + theme(
# #     axis.text = element_text(size = 20), 
# #     axis.title = element_text(size = 20))

# # pdf("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS/01_Phenotype_collation/Plots/nonverbal_reasoning.pdf")
# # plot4
# # dev.off()

# # # Vocbulary
# # d <- comp[,c("GS_id", "vocabulary")]
# # d <- na.omit(d)
# # d$verbal_total <- as.numeric(d$vocabulary)
# # mean_value <- mean(d$vocabulary)
# # bw <- 2 * IQR(d$vocabulary) / length(d$vocabulary)^(1/3)

# # plot5 <- ggplot(d, aes(x=vocabulary)) + 
# #   geom_histogram(color="black", fill="azure3", binwidth = bw) + geom_vline(xintercept = mean_value,      
# #              col = "dodgerblue4",
# #              lwd = 2) + xlab("Vocabulary") + ylab("Count") +
# #  theme_classic() + theme(
# #     axis.text = element_text(size = 20), 
# #     axis.title = element_text(size = 20))

# # pdf("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS/01_Phenotype_collation/Plots/vocabulary.pdf")
# # plot5
# # dev.off()

# # # g
# # d <- comp[,c("GS_id", "g")]
# # d <- na.omit(d)
# # d$g <- as.numeric(d$g)
# # mean_value <- mean(d$g)
# # bw <- 2 * IQR(d$g) / length(d$g)^(1/3)

# # plot6 <- ggplot(d, aes(x=g)) + 
# #   geom_histogram(color="black", fill="azure3", binwidth = bw) + geom_vline(xintercept = mean_value,      
# #              col = "dodgerblue4",
# #              lwd = 2) + xlab("General Cognitive Ability") + ylab("Count") +
# #  theme_classic() + theme(
# #     axis.text = element_text(size = 20), 
# #     axis.title = element_text(size = 20))

# # pdf("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS/01_Phenotype_collation/Plots/g.pdf")
# # plot6
# # dev.off()


# # # gf 
# # d <- comp[,c("GS_id", "gf")]
# # d <- na.omit(d)
# # d$gf <- as.numeric(d$gf)
# # mean_value <- mean(d$gf)
# # bw <- 2 * IQR(d$gf) / length(d$gf)^(1/3)

# # plot7 <- ggplot(d, aes(x=gf)) + 
# #   geom_histogram(color="black", fill="azure3", binwidth = bw) + geom_vline(xintercept = mean_value,      
# #              col = "dodgerblue4",
# #              lwd = 2) + xlab("General Fluid Cognitive Ability") + ylab("Count") +
# #  theme_classic() + theme(
# #     axis.text = element_text(size = 20), 
# #     axis.title = element_text(size = 20))

# # pdf("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS/01_Phenotype_collation/Plots/gf.pdf")
# # plot7
# # dev.off()

# # # join together as suppl plot 

# # library(patchwork)

# # pdf("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS/01_Phenotype_collation/Plots/joint_cognitive_V2.pdf", width = 15, height = 10)
# # plot1 + plot2 + plot3 + plot4 + plot5 + plot6 + plot7
# # dev.off()


# # ## Imaging 

# # # Global_GM_Volume
# # d <- comp[,c("GS_id", "Global_GM_Volume")]
# # d <- na.omit(d)
# # d$Global_GM_Volume <- as.numeric(d$Global_GM_Volume)
# # mean_value <- mean(d$Global_GM_Volume)
# # bw <- 2 * IQR(d$Global_GM_Volume) / length(d$Global_GM_Volume)^(1/3)

# # plot1 <- ggplot(d, aes(x=Global_GM_Volume)) + 
# #   geom_histogram(color="black", fill="azure3", binwidth = bw) + geom_vline(xintercept = mean_value,      
# #              col = "dodgerblue4",
# #              lwd = 2) + xlab("Global Grey Matter Volume") + ylab("Count") +
# #  theme_classic() + theme(
# #     axis.text = element_text(size = 20), 
# #     axis.title = element_text(size = 20))

# # pdf("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS/01_Phenotype_collation/Plots/global_gm_volume.pdf")
# # plot1
# # dev.off()


# # # WMHV
# # d <- comp[,c("GS_id", "WMH_Volume_Total")]
# # d <- na.omit(d)
# # d$WMH_Volume_Total <- as.numeric(d$WMH_Volume_Total)
# # mean_value <- mean(d$WMH_Volume_Total)
# # bw <- 2 * IQR(d$WMH_Volume_Total) / length(d$WMH_Volume_Total)^(1/3)

# # plot2 <- ggplot(d, aes(x=WMH_Volume_Total)) + 
# #   geom_histogram(color="black", fill="azure3", binwidth = bw) + geom_vline(xintercept = mean_value,      
# #              col = "dodgerblue4",
# #              lwd = 2) + xlab("White Matter Hyperintensity Volume") + ylab("Count") +
# #  theme_classic() + theme(
# #     axis.text = element_text(size = 20), 
# #     axis.title = element_text(size = 20))

# # pdf("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS/01_Phenotype_collation/Plots/WMHV.pdf")
# # plot2
# # dev.off()

# # # WBV
# # d <- comp[,c("GS_id", "WBV_No_Ventricles")]
# # d <- na.omit(d)
# # d$WBV_No_Ventricles <- as.numeric(d$WBV_No_Ventricles)
# # mean_value <- mean(d$WBV_No_Ventricles)
# # bw <- 2 * IQR(d$WBV_No_Ventricles) / length(d$WBV_No_Ventricles)^(1/3)

# # plot3 <- ggplot(d, aes(x=WBV_No_Ventricles)) + 
# #   geom_histogram(color="black", fill="azure3", binwidth = bw) + geom_vline(xintercept = mean_value,      
# #              col = "dodgerblue4",
# #              lwd = 2) + xlab("Whole Brain Volume") + ylab("Count") +
# #  theme_classic() + theme(
# #     axis.text = element_text(size = 20), 
# #     axis.title = element_text(size = 20))

# # pdf("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS/01_Phenotype_collation/Plots/WBV_no_vent.pdf")
# # plot3
# # dev.off()

# # # GFA
# # d <- comp[,c("GS_id", "gFA")]
# # d <- na.omit(d)
# # d$gFA <- as.numeric(d$gFA)
# # mean_value <- mean(d$gFA)
# # bw <- 2 * IQR(d$gFA) / length(d$gFA)^(1/3)

# # plot4 <- ggplot(d, aes(x=gFA)) + 
# #   geom_histogram(color="black", fill="azure3", binwidth = bw) + geom_vline(xintercept = mean_value,      
# #              col = "dodgerblue4",
# #              lwd = 2) + xlab("General Fractional Anisotropy") + ylab("Count") +
# #  theme_classic() + theme(
# #     axis.text = element_text(size = 20), 
# #     axis.title = element_text(size = 20))

# # pdf("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS/01_Phenotype_collation/Plots/gFA.pdf")
# # plot4
# # dev.off()


# # # GMD
# # d <- comp[,c("GS_id", "gMD")]
# # d <- na.omit(d)
# # d$gMD <- as.numeric(d$gMD)
# # mean_value <- mean(d$gMD)
# # bw <- 2 * IQR(d$gMD) / length(d$gMD)^(1/3)

# # plot5 <- ggplot(d, aes(x=gMD)) + 
# #   geom_histogram(color="black", fill="azure3", binwidth = bw) + geom_vline(xintercept = mean_value,      
# #              col = "dodgerblue4",
# #              lwd = 2) + xlab("General Mean Diffusivity") + ylab("Count") +
# #  theme_classic() + theme(
# #     axis.text = element_text(size = 20), 
# #     axis.title = element_text(size = 20))

# # pdf("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS/01_Phenotype_collation/Plots/gMD.pdf")
# # plot5
# # dev.off()

# # # Relative brain age 
# # d <- comp[,c("GS_id", "brain_accel")]
# # d <- na.omit(d)
# # d$brain_accel <- as.numeric(d$brain_accel)
# # mean_value <- mean(d$brain_accel)
# # bw <- 2 * IQR(d$brain_accel) / length(d$brain_accel)^(1/3)

# # plot6 <- ggplot(d, aes(x=brain_accel)) + 
# #   geom_histogram(color="black", fill="azure3", binwidth = bw) + geom_vline(xintercept = mean_value,      
# #              col = "dodgerblue4",
# #              lwd = 2) + xlab("Relative Brain Age") + ylab("Count") +
# #  theme_classic() + theme(
# #     axis.text = element_text(size = 20), 
# #     axis.title = element_text(size = 20))

# # pdf("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS/01_Phenotype_collation/Plots/Relative_brain_age.pdf")
# # plot6
# # dev.off()


# # # Fazekas WMH
# # d <- comp[,c("GS_id", "Fazekas_Score_Total")]
# # d <- na.omit(d)
# # d$Fazekas_Score_Total <- as.numeric(d$Fazekas_Score_Total)
# # mean_value <- mean(d$Fazekas_Score_Total)
# # bw <- 2 * IQR(d$Fazekas_Score_Total) / length(d$Fazekas_Score_Total)^(1/3)

# # plot7 <- ggplot(d, aes(x=Fazekas_Score_Total)) + 
# #   geom_histogram(color="black", fill="azure3", binwidth = bw) + geom_vline(xintercept = mean_value,      
# #              col = "dodgerblue4",
# #              lwd = 2) + xlab("Fazekas White Matter Hyperintensity Score") + ylab("Count") +
# #  theme_classic() + theme(
# #     axis.text = element_text(size = 20), 
# #     axis.title = element_text(size = 20))

# # pdf("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS/01_Phenotype_collation/Plots/Faz_WMH.pdf")
# # plot7
# # dev.off()

# # # join together as suppl plot 

# # library(patchwork)

# # pdf("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/PheWAS/01_Phenotype_collation/Plots/joint_imaging_V2.pdf", width = 16, height = 10)
# # plot1 + plot2 + plot3 + plot4 + plot5 + plot6 + plot7
# # dev.off()

# ###############################################################################################




