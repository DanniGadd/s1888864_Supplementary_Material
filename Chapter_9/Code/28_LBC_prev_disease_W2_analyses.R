#########################################################################################

### LBC disease binary 

#########################################################################################

# Load requisite libraries 
library(foreign)
library(data.table)


# Read in LBC data 
ph36=read.spss("/Cluster_Filespace/Marioni_Group/Elena/data/lbc_data/LBC1936_MultiOmicsOfAlcoholConsumption_RM_EB_13DEC2022.sav",to.data.frame=T)
ph36$lbc36no=gsub(" ", "", ph36$lbc36no)
# Extra from before 
b=read.spss("/Cluster_Filespace/Marioni_Group/Rob/GrimAge/LBC1936_GrimAge_RH_13MAY2019.sav",to.data.frame=T)
b=b[,c("lbc36no","bld_iron_w2","bld_hba1c_DCCT_w2","bld_creat_w2","fer_w2","pef_w2","fvc_w2")]
ph36=merge(ph36,b,by="lbc36no")

## need from other sources - HR,creatinine, g 

# variables for testing 
vars=c("lbc36no", "stroke_w2","cvdhist_w2","diab_w2") #bldcir_w2
ph36=ph36[,which(names(ph36)%in%vars),]
names(ph36)[1]="ID"

# Merge with w2 scores 
# scores=read.csv("/Cluster_Filespace/Marioni_Group/Rob/CRP/GS_LBC_Prediction/lbc_w2_output.csv")
scores <- read.csv('/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_Analysis_20k/04_LBC/original_data/W2_LBC_episcores_adjusted_technical.csv')#
names(scores)[4] <- 'ID'
df=merge(ph36,scores,by="ID")
# df$CRP=log(df$bld_hsCRP_w2+1)

# define lists to save output 
df$stroke_w2=ifelse(df$stroke_w2%in%"Yes",1,0)
df$cvdhist_w2=ifelse(df$cvdhist_w2%in%"Yes",1,0)
df$hibp_w2=ifelse(df$hibp_w2%in%"Yes",1,0)
df$diab_w2=ifelse(df$diab_w2%in%"Yes",1,0)


# > table(df$diab_w2)

#   0   1
# 706  90
# > table(df$stroke_w2)

#   0   1
# 746  50
# > table(df$cvdhist_w2)

#   0   1
# 565 231

library(tidyverse)
# Add age sex protein info in
w2 <- read.csv('/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_Analysis_20k/04_LBC/LBC_testing/GDF_w2.csv')
names(w2)[1] <- 'ID'
df <- left_join(df, w2, by = 'ID')
df$age_w2 <- df$agedays_w2.y / 365.25
df$sex <- as.factor(df$sex.y)

list3=list()
# Loop through predictors and input disease of interest as var variable 
for(j in c("GDF15_score", "Nt.proBNP_score")){ 
  for(var in vars[2:length(vars)]){  # starts at 3 due to hba1c included - didnt use here but did sensitivity with it 
    list3[[j]][[var]]<- summary(glm(df[,var] ~ scale(df[,j]) + df$age_w2 + df$sex,family="binomial"))$coefficients[2,]
  }
} 

# Summarise disease case/controls 
df2 <- df[complete.cases(df$GDF15_score),] # 796 individuals 
table(df2$stroke_w2) 
table(df2$cvdhist_w2)
table(df2$diab_w2)


# > table(df2$stroke_w2)

#   0   1
# 746  50
# > table(df2$cvdhist_w2)

#   0   1
# 565 231
# > table(df2$diab_w2)

#   0   1
# 706  90



# # Diabetes
# mod <- glm(df$diab_w2 ~ scale(df$GDF15_score) + df$age_w2 + df$sex,family="binomial")
# mod <- glm(df$diab_w2 ~ scale(df$Nt.proBNP_score) + df$age_w2 + df$sex,family="binomial")

# combine 
l4 <- do.call('c', list3)
# Tidy df 
l4=lapply(split(l4,sub('.*\\.', '', names(l4))),function(x) do.call(rbind, x))
l4=as.data.frame(do.call("rbind",l4))
l4$pred=gsub("\\..*", "", row.names(l4))
l4$trait=gsub(".*\\.", "", row.names(l4))
# Convert to OR and conf int 
l4$HCI=exp(l4[,1]+(1.96*l4[,2]))
l4$LCI=exp(l4[,1]-(1.96*l4[,2]))
l4$Odds=exp(l4[,1])
names(l4)=c("Beta","SE","t","P","Predictor","Trait","HCI","LCI","Odds")

l4 <- l4[-c(7:9),]

# Save out 
l4 <- l4[order(l4$P),]
l4$FDR <- p.adjust(l4$P, method = 'BH')
write.csv(l4, '/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_Analysis_20k/04_LBC/LBC_results_prev_basic.csv')

# Plot 

bind <- l4
bind <- bind[which(bind$FDR < 0.05),]

bind$Naming <- c('GDF15 EpiScore - Prevalent CVD',
  'GDF15 EpiScore - Prevalent type 2 Diabetes',
  'NT-proBNP - Prevalent stroke')

# Add colour column assignment
bind = bind %>% mutate(Col = case_when(
  bind$Odds < 1 ~ "royalblue",
  bind$Odds > 1 ~ "tomato2"))

# Set naming to match ordering 
bind$Naming = factor(bind$Naming, levels=unique(bind$Naming[order(bind$Odds)]))

# Set colours to match ordering 
bind$Col = factor(bind$Col, levels=unique(bind$Col[order(bind$Odds)]))

# # Sort variables for error bars
# bind$V7 <- bind$Beta - bind$SE
# bind$V8 <- bind$Beta + bind$SE

My_Theme = theme(
  axis.title.x = element_text(size = 20),
  axis.text.x = element_text(size = 20),
  axis.text.y = element_text(size = 20),
  axis.title.y = element_text(size = 24),
  strip.text = element_text(size = 20, face = "bold"),
  legend.text=element_text(size=24),
  legend.title=element_text(size=24, face = "bold"), legend.position = "none")

# bind2 <- bind[c('Odds', 'Naming', 'Col', 'LCI', 'HCI')]


pdf("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_Analysis_20k/04_LBC/PREV_JOINT_RESULTS_PLOT_dis.pdf", width = 20, height = 4)
ggplot(bind, aes(x = Odds, y = Naming, color = bind$Col)) + 
    geom_point(size = 4.5, color = bind$Col) +
        geom_errorbarh(aes(xmax = HCI, xmin = LCI), size = .9, height = 
                    .4, color = bind$Col) +
        geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") + 
    # coord_trans(x = scales:::exp_trans(10)) 
     theme_classic() + 
    My_Theme + 
    theme(panel.grid.minor = element_blank()) +
    ylab("") +
    xlab("Odds ratio") +
    ggtitle("") + theme(plot.title = element_text(hjust = 0.5, size = 27)) + xlim(0.95,2.2)
dev.off()
