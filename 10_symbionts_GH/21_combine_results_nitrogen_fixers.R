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

dsA <- read.csv("nitro_04_clr_soil.tsv", sep="\t", header=T)
dsA <- dsA %>% filter(family == "NF") %>% dplyr::select(sample,Cluster,Population,Sample.type,CLR.abd)
dsA$metric <- "CLR.abd"
dsA <- dsA %>% rename(value = CLR.abd)
head(dsA)

dsG <- read.csv("nitro_05_alpha_soil.tsv", sep="\t", header=T)
dsG$log_gini <- log(dsG$diversity_gini_simpson + 0.01)
dsG <- dsG %>% dplyr::select(SampleID, Cluster,Population,Sample.type,log_gini,chao1) %>% rename(sample=SampleID)
dssG <- dsG %>% gather(key = "metric", value = "value",-sample,-Cluster,-Population,-Sample.type)
head(dssG)


### RHIZOSPHERE

drA <- read.csv("nitro_04_clr_rhizo.tsv", sep="\t", header=T)
drA <- drA %>% filter(family == "NF") %>% dplyr::select(sample,Cluster,Population,Sample.type,CLR.abd)
drA$metric <- "CLR.abd"
drA <- drA %>% rename(value = CLR.abd)
head(drA)


drG <- read.csv("nitro_05_alpha_rhizo.tsv", sep="\t", header=T)
drG$log_gini <- log(drG$diversity_gini_simpson + 0.01)
drG <- drG %>% dplyr::select(SampleID, Cluster,Population,Sample.type,log_gini,chao1) %>% rename(sample=SampleID)
drrG <- drG %>% gather(key = "metric", value = "value",-sample,-Cluster,-Population,-Sample.type)
head(drrG)

### ROOTS

doA <- read.csv("nitro_04_clr_roots.tsv", sep="\t", header=T)
doA <- doA %>% filter(family == "NF") %>% dplyr::select(sample,Cluster,Population,Sample.type,CLR.abd)
doA$metric <- "CLR.abd"
doA <- doA %>% rename(value = CLR.abd)
head(doA)


doG <- read.csv("nitro_05_alpha_roots.tsv", sep="\t", header=T)
doG$log_gini <- log(doG$diversity_gini_simpson + 0.01)
doG <- doG %>% dplyr::select(SampleID, Cluster,Population,Sample.type,log_gini,chao1) %>% rename(sample=SampleID)
dooG <- doG %>% gather(key = "metric", value = "value",-sample,-Cluster,-Population,-Sample.type)
head(dooG)

### Combining tables

dm <- rbind(dsA, dssG, drA, drrG, doA, dooG)
head(dm)






####### TEST RESULTS

# SOIL
dsAT <- read.csv("nitro_04_clr_soil_tests.tsv", sep="\t", header=T)
dsAT <- dsAT %>% filter(EVAL == "OK")
dsAT$metric <- "CLR.abd"
dsAT$Sample.type <- "bulk soil"
head(dsAT)

dsGT <- read.csv("nitro_05_alpha_soil_tests.tsv", sep="\t", header=T)
dsGT <- dsGT %>% filter(EVAL == "OK")
dsGT$metric <- "log_gini"
dsGT$Sample.type <- "bulk soil"
head(dsGT)

dsRT <- read.csv("nitro_06_chao_soil_tests.tsv", sep="\t", header=T)
dsRT <- dsRT %>% filter(EVAL == "OK")
dsRT$metric <- "chao1"
dsRT$Sample.type <- "bulk soil"
head(dsRT)

# RHIZO

drAT <- read.csv("nitro_04_clr_rhizo_tests.tsv", sep="\t", header=T)
drAT <- drAT %>% filter(EVAL == "OK")
drAT$metric <- "CLR.abd"
drAT$Sample.type <- "rhizosphere"
head(drAT)

drGT <- read.csv("nitro_05_alpha_rhizo_tests.tsv", sep="\t", header=T)
drGT <- drGT %>% filter(EVAL == "OK")
drGT$metric <- "log_gini"
drGT$Sample.type <- "rhizosphere"
head(drGT)

drRT <- read.csv("nitro_06_chao_rhizo_tests.tsv", sep="\t", header=T)
drRT <- drRT %>% filter(EVAL == "OK")
drRT$metric <- "chao1"
drRT$Sample.type <- "rhizosphere"
head(drRT)


# ROOTS

doAT <- read.csv("nitro_04_clr_roots_tests.tsv", sep="\t", header=T)
doAT <- doAT %>% filter(EVAL == "OK")
doAT$metric <- "CLR.abd"
doAT$Sample.type <- "roots"
head(doAT)

doGT <- read.csv("nitro_05_alpha_roots_tests.tsv", sep="\t", header=T)
doGT <- doGT %>% filter(EVAL == "OK")
doGT$metric <- "log_gini"
doGT$Sample.type <- "roots"
head(doGT)

doRT <- read.csv("nitro_06_chao_roots_tests.tsv", sep="\t", header=T)
doRT <- drRT %>% filter(EVAL == "OK")
doRT$metric <- "chao1"
doRT$Sample.type <- "roots"
head(doRT)

### Combining tables

dt <- rbind(dsAT, dsGT, dsRT, drAT, drGT, drRT, doAT, doGT, doRT)
head(dt)





### Saving
write.table(dm, file="nitro_07_combine_inputs.tsv",sep="\t", row.names = F, col.names = T, append=F, quote=F)
write.table(dt, file="nitro_07_combine_tests.tsv",sep="\t", row.names = F, col.names = T, append=F, quote=F)








