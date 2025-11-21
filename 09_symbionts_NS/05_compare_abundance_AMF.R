rm(list = ls())


library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(phyloseq)
library(microbiome)
library(readxl)
library(stringr)
library(rstatix)

source("run_NESTED_ANOVA_EXT.R")

#########################
#--------DATA-----------#
#########################


df <- read.csv("./guilds_01_read_stands.csv", header=TRUE, sep=',')
head(df)

df.meta <- read.csv("./guilds_01_read_stands_meta.csv", header=TRUE, sep=',')
df.meta$long.name <- paste0("X",gsub('-','\\.',df.meta$sample))
head(df.meta)
df.meta$long.name

df.stands <- df.meta %>% filter(location == "natural_stands") %>% filter(source == "Soil")
ids.stands <- df.stands %>% dplyr::select(long.name) %>% distinct() %>% pull(long.name)
ids.stands

df.greens <- df.meta %>% filter(location == "greenhouse") %>% filter(source == "Soil")
ids.greens <- df.greens %>% dplyr::select(long.name) %>% distinct() %>% pull(long.name)
ids.greens



# Getting otu_table - STANDS only
head(df)
colnames(df)

# Subset stands only
df_ <- df %>% dplyr::select(all_of(c("group",ids.stands)))
df_
df_$OTU <- c(1:length(df_$group))
rownames(df_) <- paste0('group',df_$OTU)
df_
otu.table <- df_ %>% dplyr::select(-OTU, -group)
head(otu.table)
OTU = otu_table(otu.table, taxa_are_rows = TRUE)

df.tax <- df_ %>% dplyr::select(group)
colnames(df.tax) <- c("Phylum")
df.tax.mat <- df.tax %>% as.matrix()
TAX = tax_table(df.tax.mat)
TAX

head(df.meta)
#df.meta.format <- df.meta %>% dplyr::rename(long.name = long.name.sequencing_R.,
#                                            Country = Country_x)
df.meta.format <- df.meta %>% dplyr::select(long.name, short_name, site, cluster) %>% distinct()
rownames(df.meta.format) <- df.meta.format$long.name
head(df.meta.format)
samdata <- sample_data(df.meta.format)

# Importing to phyloseq
pseq.stands <- phyloseq(OTU, TAX, sam_data = samdata)



# Getting otu_table - GREENHOUSE
# df.2 <- df %>% dplyr::select(all_of(c("group",ids.greens)))
# df.2$OTU <- c(1:length(df.2$group))
# rownames(df.2) <- paste0('group',df.2$OTU)
# otu.table <- df.2 %>% dplyr::select(-OTU, -group)
# head(otu.table)
# OTU = otu_table(otu.table, taxa_are_rows = TRUE)
# colnames(OTU)
# 
# df.tax <- df.2 %>% dplyr::select(group)
# colnames(df.tax) <- c("Phylum")
# df.tax.mat <- df.tax %>% as.matrix()
# TAX = tax_table(df.tax.mat)
# TAX
# 
# head(df.meta)
# #df.meta.format <- df.2.meta %>% dplyr::rename(Site = Cluster)
# df.meta.format <- df.meta %>% dplyr::select(sample, long.name, short_name, source, cluster, site) %>% distinct()
# #df.meta.format$Name <- chartr("-", ".", df.meta.format$SampleID)
# rownames(df.meta.format) <- df.meta.format$long.name
# head(df.meta.format)
# samdata <- sample_data(df.meta.format)
# 
# # Importing to phyloseq
# pseq.greenhouse <- phyloseq(OTU, TAX, sam_data = samdata)






#####################
#---  TRANSFORM  ---#
#####################


sample_data(pseq.stands)
#sample_data(pseq.greenhouse)


# Convert to tab

#### CLR Transform


# Transform
pseq.stands.clr <- microbiome::transform(pseq.stands, 'clr')
pseq.stands.clr@otu_table

# Convert to dataframe
dtab.stands <- as.data.frame(pseq.stands.clr@otu_table)
class(dtab.stands)
head(dtab.stands)

# Gather
dtab.stands$group <- rownames(dtab.stands)
dg.stands <- dtab.stands %>% gather(key = "sample", value = "CLR abd",-group)
head(dg.stands)


# Add meta
df.meta.format <- data.frame(sample_data(pseq.stands.clr))
head(df.meta.format)
dg.meta.stands <- merge(dg.stands, df.meta.format, by.x = "sample", by.y = "long.name", sort=FALSE)
head(dg.meta.stands)

# Add taxonomy
df.taxonomy <- as.data.frame(pseq.stands.clr@tax_table)
df.taxonomy$group <- rownames(df.taxonomy)

dm.stands <- merge(dg.meta.stands, df.taxonomy, by = 'group', sort=FALSE)
head(dm.stands)










#####################
#-----  PLOTS ------#
#####################

dm.stands$Phylum

guilds <- c("Arbuscular Mycorrhizal")
dsub.stands <- dm.stands %>% filter(Phylum %in% guilds) 
dsub.stands %>% head()

unique(dsub.stands$Phylum)
dsub.stands$Label <- factor(dsub.stands$Phylum, levels=c("Arbuscular Mycorrhizal"))

dsub.stands$region1 <- ifelse(dsub.stands$site %in% c("AMOS","STFE"),"Boreal","Rest")
dsub.stands$region2 <- ifelse(dsub.stands$site %in% c("STET","ESSI","FORE"),"Cold_temperate",dsub.stands$region1)
dsub.stands$region <- ifelse(dsub.stands$site %in% c("Santiago","FLOR1","FP2"),"Warm_temperate",dsub.stands$region2)
dsub.stands$region <- factor(dsub.stands$region, levels=c("Boreal","Cold_temperate","Warm_temperate"))
head(dsub.stands)


write.table(dsub.stands, file = "guilds_04_clr_AMF_inputs.tsv", sep="\t", row.names = F, col.names = T, append=F, quote=F)




#####################
#-----   SITE ------#
#####################
# Plotting average abundance for each site

head(dm.stands)
dm.stands$site <- factor(dm.stands$site, levels = c("AMOS","STFE","STET","ESSI","FORE","Santiago","FLOR1","FP2"))
dm.stands$Country <- factor(dm.stands$cluster, levels = c("NENA","MX"))

p1 <- ggplot(dm.stands) +
  geom_boxplot(data=dm.stands, aes(x = site, y = `CLR abd`, fill = cluster), outlier.shape = NA) +
  geom_point(data=dm.stands, aes(x = site, y = `CLR abd`), color="grey40", pch=21, size=1) +
  facet_wrap(~Phylum)
p1




###############
### TESTS of difference between countries

head(dsub.stands)

# Transform CLR so it's positive
min_clr <- min(dsub.stands$`CLR abd`)
dsub.stands$CLR_abd <- dsub.stands$`CLR abd` + abs(min_clr) + 0.01
head(dsub.stands)

### Raw model

df.models <- run_nested_models(dsub.stands, "CLR_abd", "cluster", "site")
df.models

write.table(df.models, file = "guilds_04_clr_AMF_country_tests.tsv", sep="\t", row.names = F, col.names = T, append=F, quote=F)






###############
### TESTS of difference between REGIONS

df.models <- run_nested_models(dsub.stands, "CLR_abd", "region", "site")
df.models

write.table(df.models, file = "guilds_04_clr_AMF_region_tests.tsv", sep="\t", row.names = F, col.names = T, append=F, quote=F)



