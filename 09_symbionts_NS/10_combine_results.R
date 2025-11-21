rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(cowplot)


#########################
#--------DATA-----------#
#########################



### EMF

dsA <- read.csv("../08_guilds/guilds_04_clr_EMF_inputs.tsv", sep="\t", header=T)
dsA <- dsA %>% filter(Phylum == "Ectomycorrhizal") %>% dplyr::select(sample,site,Country,region,Phylum,CLR.abd)
dsA$metric <- "CLR.abd"
dsA <- dsA %>% rename(value = CLR.abd)
head(dsA)

dsG <- read.csv("../08_guilds/guilds_05_alpha_EMF_inputs.tsv", sep="\t", header=T)
dsG$Phylum <- c("Ectomycorrhizal")
head(dsG)
dsG$log_gini <- log(dsG$diversity_gini_simpson + 0.01)
dsG <- dsG %>% dplyr::select(long.name, site,Country,region,Phylum,log_gini,chao1) %>% rename(sample=long.name)
dssG <- dsG %>% gather(key = "metric", value = "value",-sample,-site,-Country,-region,-Phylum)
head(dssG)


### AMF

drA <- read.csv("../09_guilds_amf/guilds_04_clr_AMF_inputs.tsv", sep="\t", header=T)
drA$Country <- ifelse(drA$cluster == "NENA","Canada","Mexico")
drA <- drA %>% filter(Phylum == "Arbuscular Mycorrhizal") %>% dplyr::select(sample,site,Country,region, Phylum,CLR.abd)
drA$metric <- "CLR.abd"
drA <- drA %>% rename(value = CLR.abd)
head(drA)

drG <- read.csv("../09_guilds_amf/guilds_05_alpha_AMF_inputs.tsv", sep="\t", header=T)
drG$log_gini <- log(drG$diversity_gini_simpson + 0.01)
drG$Country <- ifelse(drG$cluster == "NENA","Canada","Mexico")
drG$Phylum <- c("Arbuscular Mycorrhizal")
head(drG)
drG <- drG %>% dplyr::select(sampleid, site,Country,region,Phylum,log_gini,chao1) %>% rename(sample=sampleid)
drrG <- drG %>% gather(key = "metric", value = "value",-sample,-site,-Country,-region,-Phylum)
head(drrG)

### NITRO

doA <- read.csv("../10_nitro/nitro_04_clr_inputs.tsv", sep="\t", header=T)
doA$Phylum <- "NF"
doA <- doA %>% filter(Phylum == "NF") %>% dplyr::select(sample,site,Country,region,Phylum,CLR.abd)
doA$metric <- "CLR.abd"
doA <- doA %>% rename(value = CLR.abd)
head(doA)

doG <- read.csv("../10_nitro/nitro_05_alpha_inputs.tsv", sep="\t", header=T)
doG$log_gini <- log(doG$diversity_gini_simpson + 0.01)
doG$Phylum <- "NF"
doG <- doG %>% dplyr::select(long.name, site,Country,region,Phylum,log_gini,chao1) %>% rename(sample=long.name)
dooG <- doG %>% gather(key = "metric", value = "value",-sample,-site,-Country,-region,-Phylum)
head(dooG)


### Combining tables

dm <- rbind(dsA, dssG, drA, drrG, doA, dooG)
head(dm)






####### TEST RESULTS

# EMF
dsAT <- read.csv("../08_guilds/guilds_04_clr_EMF_country_tests.tsv", sep="\t", header=T)
dsAT <- dsAT %>% filter(family == "Ectomycorrhizal" & EVAL == "OK") %>% dplyr::select(-family)
dsAT$metric <- "CLR.abd"
dsAT$Phylum <- "Ectomycorrhizal"
dsAT$Group <- "Country"
head(dsAT)

dsAT2 <- read.csv("../08_guilds/guilds_04_clr_EMF_region_tests.tsv", sep="\t", header=T)
dsAT2 <- dsAT2 %>% filter(family == "Ectomycorrhizal" & EVAL == "OK") %>% dplyr::select(-family)
dsAT2$metric <- "CLR.abd"
dsAT2$Phylum <- "Ectomycorrhizal"
dsAT2$Group <- "region"
head(dsAT2)

dsGT <- read.csv("../08_guilds/guilds_05_alpha_country_tests.tsv", sep="\t", header=T)
dsGT <- dsGT %>% filter(EVAL == "OK")
dsGT$metric <- "log_gini"
dsGT$Phylum <- "Ectomycorrhizal"
dsGT$Group <- "Country"
head(dsGT)

dsGT2 <- read.csv("../08_guilds/guilds_05_alpha_region_tests.tsv", sep="\t", header=T)
dsGT2 <- dsGT2 %>% filter(EVAL == "OK")
dsGT2$metric <- "log_gini"
dsGT2$Phylum <- "Ectomycorrhizal"
dsGT2$Group <- "region"
head(dsGT2)

dsRT <- read.csv("../08_guilds/guilds_05_chao1_country_tests.tsv", sep="\t", header=T)
dsRT <- dsRT %>% filter(EVAL == "OK")
dsRT$metric <- "chao1"
dsRT$Phylum <- "Ectomycorrhizal"
dsRT$Group <- "Country"
head(dsRT)

dsRT2 <- read.csv("../08_guilds/guilds_05_chao1_region_tests.tsv", sep="\t", header=T)
dsRT2 <- dsRT2 %>% filter(EVAL == "OK")
dsRT2$metric <- "chao1"
dsRT2$Phylum <- "Ectomycorrhizal"
dsRT2$Group <- "region"
head(dsRT2)

ds.emf <- rbind(dsAT, dsAT2, dsGT, dsGT2, dsRT, dsRT2)



# AMF

dsAT <- read.csv("../09_guilds_amf/guilds_04_clr_AMF_country_tests.tsv", sep="\t", header=T)
dsAT <- dsAT %>% filter(EVAL == "OK")
dsAT$metric <- "CLR.abd"
dsAT$Phylum <- "Arbuscular Mycorrhizal"
dsAT$Group <- "Country"
head(dsAT)

dsAT2 <- read.csv("../09_guilds_amf/guilds_04_clr_AMF_region_tests.tsv", sep="\t", header=T)
dsAT2 <- dsAT2 %>% filter(EVAL == "OK")
dsAT2$metric <- "CLR.abd"
dsAT2$Phylum <- "Arbuscular Mycorrhizal"
dsAT2$Group <- "region"
head(dsAT2)

dsGT <- read.csv("../09_guilds_amf/guilds_05_alpha_AMF_country_tests.tsv", sep="\t", header=T)
dsGT <- dsGT %>% filter(EVAL == "OK")
dsGT$metric <- "log_gini"
dsGT$Phylum <- "Arbuscular Mycorrhizal"
dsGT$Group <- "Country"
head(dsGT)

dsGT2 <- read.csv("../09_guilds_amf/guilds_05_alpha_AMF_region_tests.tsv", sep="\t", header=T)
dsGT2 <- dsGT2 %>% filter(EVAL == "OK")
dsGT2$metric <- "log_gini"
dsGT2$Phylum <- "Arbuscular Mycorrhizal"
dsGT2$Group <- "region"
head(dsGT2)

dsRT <- read.csv("../09_guilds_amf/guilds_05_chao1_country_tests.tsv", sep="\t", header=T)
dsRT <- dsRT %>% filter(EVAL == "OK")
dsRT$metric <- "chao1"
dsRT$Phylum <- "Arbuscular Mycorrhizal"
dsRT$Group <- "Country"
head(dsRT)

dsRT2 <- read.csv("../09_guilds_amf/guilds_05_chao1_region_tests.tsv", sep="\t", header=T)
dsRT2 <- dsRT2 %>% filter(EVAL == "OK")
dsRT2$metric <- "chao1"
dsRT2$Phylum <- "Arbuscular Mycorrhizal"
dsRT2$Group <- "region"
head(dsRT2)

ds.amf <- rbind(dsAT, dsAT2, dsGT, dsGT2, dsRT, dsRT2)


# NITRO

dsAT <- read.csv("../10_nitro/nitro_04_clr_country_tests.tsv", sep="\t", header=T)
dsAT <- dsAT %>% filter(EVAL == "OK")
dsAT$metric <- "CLR.abd"
dsAT$Phylum <- "NF"
dsAT$Group <- "Country"
head(dsAT)

dsAT2 <- read.csv("../10_nitro/nitro_04_clr_region_tests.tsv", sep="\t", header=T)
dsAT2 <- dsAT2 %>% filter(EVAL == "OK")
dsAT2$metric <- "CLR.abd"
dsAT2$Phylum <- "NF"
dsAT2$Group <- "region"
head(dsAT2)

dsGT <- read.csv("../10_nitro/nitro_05_alpha_country_tests.tsv", sep="\t", header=T)
dsGT <- dsGT %>% filter(EVAL == "OK")
dsGT$metric <- "log_gini"
dsGT$Phylum <- "NF"
dsGT$Group <- "Country"
head(dsGT)

dsGT2 <- read.csv("../10_nitro/nitro_05_alpha_region_tests.tsv", sep="\t", header=T)
dsGT2 <- dsGT2 %>% filter(EVAL == "OK")
dsGT2$metric <- "log_gini"
dsGT2$Phylum <- "NF"
dsGT2$Group <- "region"
head(dsGT2)

dsRT <- read.csv("../10_nitro/nitro_05_chao1_country_tests.tsv", sep="\t", header=T)
dsRT <- dsRT %>% filter(EVAL == "OK")
dsRT$metric <- "chao1"
dsRT$Phylum <- "NF"
dsRT$Group <- "Country"
head(dsRT)

dsRT2 <- read.csv("../10_nitro/nitro_05_chao1_region_tests.tsv", sep="\t", header=T)
dsRT2 <- dsRT2 %>% filter(EVAL == "OK")
dsRT2$metric <- "chao1"
dsRT2$Phylum <- "NF"
dsRT2$Group <- "region"
head(dsRT2)

ds.nf <- rbind(dsAT, dsAT2, dsGT, dsGT2, dsRT, dsRT2)

### Combining tables

dt <- rbind(ds.emf, ds.amf, ds.nf)
head(dt)





### Saving
write.table(dm, file="guilds_05_combine_inputs.tsv",sep="\t", row.names = F, col.names = T, append=F, quote=F)
write.table(dt, file="guilds_05_combine_tests.tsv",sep="\t", row.names = F, col.names = T, append=F, quote=F)











### CLR guilds only

### EMF

dsA <- read.csv("../08_guilds/guilds_04_clr_EMF_inputs.tsv", sep="\t", header=T)
dsA <- dsA %>% dplyr::select(sample,site,Country,region,Phylum,CLR.abd)
dsA$metric <- "CLR.abd"
dsA <- dsA %>% dplyr::rename(value = CLR.abd)
head(dsA)

### AMF

drA <- read.csv("../09_guilds_amf/guilds_04_clr_AMF_inputs.tsv", sep="\t", header=T)
drA$Country <- ifelse(drA$cluster == "NENA","Canada","Mexico")
drA <- drA %>% dplyr::select(sample,site,Country,region, Phylum,CLR.abd)
drA$metric <- "CLR.abd"
drA <- drA %>% dplyr::rename(value = CLR.abd)
head(drA)

### NITRO

doA <- read.csv("../10_nitro/nitro_04_clr_inputs.tsv", sep="\t", header=T)
doA$Phylum <- "NF"
doA <- doA %>% dplyr::select(sample,site,Country,region,Phylum,CLR.abd)
doA$metric <- "CLR.abd"
doA <- doA %>% dplyr::rename(value = CLR.abd)
head(doA)

### Combining tables

dm2 <- rbind(dsA, drA, doA)
head(dm2)



####### TEST RESULTS


# EMF
dsAT <- read.csv("../08_guilds/guilds_04_clr_EMF_region_tests.tsv", sep="\t", header=T)
dsAT <- dsAT %>% filter(EVAL == "OK")
dsAT$metric <- "CLR.abd"
dsAT$Phylum <- dsAT$family
dsAT$Group <- "region"
dsAT <- dsAT %>% dplyr::select(-family)
head(dsAT)

# AMF

dsAT2 <- read.csv("../09_guilds_amf/guilds_04_clr_AMF_region_tests.tsv", sep="\t", header=T)
dsAT2 <- dsAT2 %>% filter(EVAL == "OK")
dsAT2$metric <- "CLR.abd"
dsAT2$Phylum <- "Arbuscular Mycorrhizal"
dsAT2$Group <- "region"
head(dsAT2)

# NITRO

dsAT3 <- read.csv("../10_nitro/nitro_04_clr_region_tests.tsv", sep="\t", header=T)
dsAT3 <- dsAT3 %>% filter(EVAL == "OK")
dsAT3$metric <- "CLR.abd"
dsAT3$Phylum <- "NF"
dsAT3$Group <- "region"
head(dsAT3)



### Combining tables

dt2 <- rbind(dsAT, dsAT2, dsAT3)
head(dt2)



### Saving
write.table(dm2, file="guilds_05_combine_guilds_inputs.tsv",sep="\t", row.names = F, col.names = T, append=F, quote=F)
write.table(dt2, file="guilds_05_combine_guilds_tests.tsv",sep="\t", row.names = F, col.names = T, append=F, quote=F)


