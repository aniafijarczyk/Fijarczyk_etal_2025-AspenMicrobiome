rm(list = ls())

library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(phyloseq)
library(microbiome)
library(readxl)
library(stringr)
library(rstatix) # Welsh anova
library(nlme) #gls()

source("run_NESTED_ANOVA_EXT.R")



#########################
#--------DATA-----------#
#########################


df <- read.csv("./nitro_02_read.csv", header=TRUE, sep=',')
head(df)
df.meta <- read.csv("./nitro_02_read_meta.csv", header=TRUE, sep=',')
head(df.meta)


rownames(df) <- df$NF
otu.table <- df %>% dplyr::select(-NF)
head(otu.table)
OTU = otu_table(otu.table, taxa_are_rows = TRUE)

df.tax <- df %>% dplyr::select(NF)
colnames(df.tax) <- c("Family")
df.tax.mat <- df.tax %>% as.matrix()
TAX = tax_table(df.tax.mat)
TAX



df.meta.format <- df.meta %>% dplyr::select(long.name, short.name, site, Cluster, Country, Lattitude, Longitude) %>% distinct()
rownames(df.meta.format) <- df.meta.format$long.name
head(df.meta.format)
samdata <- sample_data(df.meta.format)

# Importing to phyloseq
pseq.stands <- phyloseq(OTU, TAX, sam_data = samdata)





#####################
#---  TRANSFORM  ---#
#####################


pseq.stands.rel.abd <- microbiome::transform(pseq.stands, 'compositional')
nfs <- as.data.frame(t(pseq.stands.rel.abd@otu_table))
summary(nfs$NF)

# Convert to tab

# Transform
pseq.stands.rel <- microbiome::transform(pseq.stands, 'clr')
pseq.stands.rel@otu_table


# Convert to dataframe
dtab.stands <- as.data.frame(pseq.stands.rel@otu_table)
class(dtab.stands)
head(dtab.stands)

# Gather
dtab.stands$family <- rownames(dtab.stands)
dg.stands <- dtab.stands %>% gather(key = "sample", value = "CLR abd",-family)
dg.stands

# Add meta
dg.meta.stands <- merge(dg.stands, df.meta.format, by.x = "sample", by.y = "long.name", sort=FALSE)
head(dg.meta.stands)

# Add taxonomy
df.taxonomy <- as.data.frame(pseq.stands.rel@tax_table)
df.taxonomy$family <- rownames(df.taxonomy)

dm.stands <- merge(dg.meta.stands, df.taxonomy, by = 'family', sort=FALSE)
head(dm.stands)






#####################
#-----  PLOTS ------#
#####################

fams <- c('NF')

dsub.stands <- dm.stands %>% filter(family %in% fams) 
dsub.stands %>% head()

dsub.stands$Label <- factor(dsub.stands$Family,levels=c("NF"))

dsub.stands$region1 <- ifelse(dsub.stands$site %in% c("AMOS","STFE"),"Boreal","Rest")
dsub.stands$region2 <- ifelse(dsub.stands$site %in% c("STET","ESSI","FORE"),"Cold_temperate",dsub.stands$region1)
dsub.stands$region <- ifelse(dsub.stands$site %in% c("Santiago","FLOR1","FP2"),"Warm_temperate",dsub.stands$region2)
dsub.stands$region <- factor(dsub.stands$region, levels=c("Boreal","Cold_temperate","Warm_temperate"))
head(dsub.stands)


write.table(dsub.stands, file = "nitro_04_clr_inputs.tsv", sep="\t", row.names = F, col.names = T, append=F, quote=F)




#####################
#-----   SITE ------#
#####################

# Plotting average abundance for each site

head(dsub.stands)
dsub.stands$site <- factor(dsub.stands$site, levels = c("AMOS","STFE","STET","ESSI","FORE","Santiago","FLOR1","FP2"))

p1 <- ggplot(dsub.stands) +
  geom_boxplot(data=dsub.stands, aes(x = site, y = `CLR abd`, fill = Country), outlier.shape = NA) +
  geom_point(data=dsub.stands, aes(x = site, y = `CLR abd`), color="grey40", pch=21, size=1)

p1




#####################
#-----  TESTS  -----#
#####################



###############
### TESTS of difference between countries

head(dsub.stands)

# Transform CLR so it's positive
min_clr <- min(dsub.stands$`CLR abd`)
dsub.stands$CLR_abd <- dsub.stands$`CLR abd` + abs(min_clr) + 0.01
head(dsub.stands)
df.models <- run_nested_models(dsub.stands, "CLR_abd", "Country", "site")
df.models


write.table(df.models, file = "nitro_04_clr_country_tests.tsv", sep="\t", row.names = F, col.names = T, append=F, quote=F)


###############
### TESTS of difference between REGIONS

# Transform CLR so it's positive
min_clr <- min(dsub.stands$`CLR abd`)
dsub.stands$CLR_abd <- dsub.stands$`CLR abd` + abs(min_clr) + 0.01
head(dsub.stands)
df.models <- run_nested_models(dsub.stands, "CLR_abd", "region", "site")
df.models




write.table(df.models, file = "nitro_04_clr_region_tests.tsv", sep="\t", row.names = F, col.names = T, append=F, quote=F)


















