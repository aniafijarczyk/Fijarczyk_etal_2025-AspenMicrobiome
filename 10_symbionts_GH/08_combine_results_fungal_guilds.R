rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(cowplot)


#########################
#--------DATA-----------#
#########################



### SOIL

dsA <- read.csv("guilds_03_clr_soil_inputs.tsv", sep="\t", header=T)
dsA <- dsA %>% filter(family == "Ectomycorrhizal") %>% dplyr::select(sample,Cluster,Population,Sample.type,CLR.abd)
dsA$metric <- "CLR.abd"
dsA <- dsA %>% rename(value = CLR.abd)
head(dsA)

dsG <- read.csv("guilds_04_alpha_soil.tsv", sep="\t", header=T)
dsG$log_gini <- log(dsG$diversity_gini_simpson + 0.01)
dsG <- dsG %>% dplyr::select(SampleID, Cluster,Population,Sample.type,log_gini,chao1) %>% rename(sample=SampleID)
dssG <- dsG %>% gather(key = "metric", value = "value",-sample,-Cluster,-Population,-Sample.type)
head(dssG)


### RHIZOSPHERE

drA <- read.csv("guilds_03_clr_rhizo_inputs.tsv", sep="\t", header=T)
drA <- drA %>% filter(family == "Ectomycorrhizal") %>% dplyr::select(sample,Cluster,Population,Sample.type,CLR.abd)
drA$metric <- "CLR.abd"
drA <- drA %>% rename(value = CLR.abd)
head(drA)


drG <- read.csv("guilds_04_alpha_rhizo.tsv", sep="\t", header=T)
drG$log_gini <- log(drG$diversity_gini_simpson + 0.01)
drG <- drG %>% dplyr::select(SampleID, Cluster,Population,Sample.type,log_gini,chao1) %>% rename(sample=SampleID)
drrG <- drG %>% gather(key = "metric", value = "value",-sample,-Cluster,-Population,-Sample.type)
head(drrG)

### ROOTS

doA <- read.csv("guilds_03_clr_roots_inputs.tsv", sep="\t", header=T)
doA <- doA %>% filter(family == "Ectomycorrhizal") %>% dplyr::select(sample,Cluster,Population,Sample.type,CLR.abd)
doA$metric <- "CLR.abd"
doA <- doA %>% rename(value = CLR.abd)
head(doA)


doG <- read.csv("guilds_04_alpha_roots.tsv", sep="\t", header=T)
doG$log_gini <- log(doG$diversity_gini_simpson + 0.01)
doG <- doG %>% dplyr::select(SampleID, Cluster,Population,Sample.type,log_gini,chao1) %>% rename(sample=SampleID)
dooG <- doG %>% gather(key = "metric", value = "value",-sample,-Cluster,-Population,-Sample.type)
head(dooG)

### Combining tables

dm <- rbind(dsA, dssG, drA, drrG, doA, dooG)
head(dm)






####### TEST RESULTS

# SOIL
dsAT <- read.csv("guilds_03_clr_soil_tests.tsv", sep="\t", header=T)
dsAT <- dsAT %>% filter(family == "Ectomycorrhizal" & EVAL == "OK") %>% dplyr::select(-family)
dsAT$metric <- "CLR.abd"
dsAT$Sample.type <- "bulk soil"
head(dsAT)

dsGT <- read.csv("guilds_04_alpha_soil_tests.tsv", sep="\t", header=T)
dsGT <- dsGT %>% filter(EVAL == "OK")
dsGT$metric <- "log_gini"
dsGT$Sample.type <- "bulk soil"
head(dsGT)

dsRT <- read.csv("guilds_04_chao1_soil_tests.tsv", sep="\t", header=T)
dsRT <- dsRT %>% filter(EVAL == "OK")
dsRT$metric <- "chao1"
dsRT$Sample.type <- "bulk soil"
head(dsRT)

# RHIZO

drAT <- read.csv("guilds_03_clr_rhizo_tests.tsv", sep="\t", header=T)
drAT <- drAT %>% filter(family == "Ectomycorrhizal" & EVAL == "OK") %>% dplyr::select(-family)
drAT$metric <- "CLR.abd"
drAT$Sample.type <- "rhizosphere"
head(drAT)

drGT <- read.csv("guilds_04_alpha_rhizo_tests.tsv", sep="\t", header=T)
drGT$EVAL <- "OK"
drGT <- drGT %>% filter(EVAL == "OK")
drGT$metric <- "log_gini"
drGT$Sample.type <- "rhizosphere"
head(drGT)

drRT <- read.csv("guilds_04_chao1_rhizo_tests.tsv", sep="\t", header=T)
drRT <- drRT %>% filter(EVAL == "OK")
drRT$metric <- "chao1"
drRT$Sample.type <- "rhizosphere"
head(drRT)


# ROOTS

doAT <- read.csv("guilds_03_clr_roots_tests.tsv", sep="\t", header=T)
doAT <- doAT %>% filter(family == "Ectomycorrhizal" & EVAL == "OK") %>% dplyr::select(-family)
doAT$metric <- "CLR.abd"
doAT$Sample.type <- "roots"
head(doAT)

doGT <- read.csv("guilds_04_alpha_roots_tests.tsv", sep="\t", header=T)
doGT$EVAL <- "OK"
doGT <- doGT %>% filter(EVAL == "OK")
doGT$metric <- "log_gini"
doGT$Sample.type <- "roots"
head(doGT)

doRT <- read.csv("guilds_04_chao1_roots_tests.tsv", sep="\t", header=T)
doRT <- drRT %>% filter(EVAL == "OK")
doRT$metric <- "chao1"
doRT$Sample.type <- "roots"
head(doRT)

### Combining tables

dt <- rbind(dsAT, dsGT, dsRT, drAT, drGT, drRT, doAT, doGT, doRT)
head(dt)





### Saving
write.table(dm, file="guilds_05_combine_inputs.tsv",sep="\t", row.names = F, col.names = T, append=F, quote=F)
write.table(dt, file="guilds_05_combine_tests.tsv",sep="\t", row.names = F, col.names = T, append=F, quote=F)











### CLR guilds only

### SOIL

dsA <- read.csv("guilds_03_clr_soil_inputs.tsv", sep="\t", header=T)
dsA <- dsA %>% dplyr::select(sample,Cluster,Population,Sample.type,CLR.abd, family)
dsA$metric <- "CLR.abd"
dsA <- dsA %>% rename(value = CLR.abd)

drA <- read.csv("guilds_03_clr_rhizo_inputs.tsv", sep="\t", header=T)
drA <- drA %>% dplyr::select(sample,Cluster,Population,Sample.type,CLR.abd, family)
drA$metric <- "CLR.abd"
drA <- drA %>% rename(value = CLR.abd)

doA <- read.csv("guilds_03_clr_roots_inputs.tsv", sep="\t", header=T)
doA <- doA %>% dplyr::select(sample,Cluster,Population,Sample.type,CLR.abd, family)
doA$metric <- "CLR.abd"
doA <- doA %>% rename(value = CLR.abd)
head(doA)


### Combining tables

dm2 <- rbind(dsA, drA, doA)
head(dm2)



####### TEST RESULTS

# SOIL
dsAT <- read.csv("guilds_03_clr_soil_tests.tsv", sep="\t", header=T)
dsAT <- dsAT %>% filter(EVAL == "OK")
dsAT$metric <- "CLR.abd"
dsAT$Sample.type <- "bulk soil"
head(dsAT)


# RHIZO

drAT <- read.csv("guilds_03_clr_rhizo_tests.tsv", sep="\t", header=T)
drAT <- drAT %>% filter(EVAL == "OK")
drAT$metric <- "CLR.abd"
drAT$Sample.type <- "rhizosphere"
head(drAT)

# ROOTS

doAT <- read.csv("guilds_03_clr_roots_tests.tsv", sep="\t", header=T)
doAT <- doAT %>% filter(EVAL == "OK")
doAT$metric <- "CLR.abd"
doAT$Sample.type <- "roots"
head(doAT)

### Combining tables

dt2 <- rbind(dsAT, drAT, doAT)
head(dt2)



### Saving
write.table(dm2, file="guilds_05_combine_guilds_inputs.tsv",sep="\t", row.names = F, col.names = T, append=F, quote=F)
write.table(dt2, file="guilds_05_combine_guilds_tests.tsv",sep="\t", row.names = F, col.names = T, append=F, quote=F)


