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


df <- read.csv("./guilds_01_read_stands.csv", header=TRUE, sep=',')
head(df)
df.meta <- read.csv("./guilds_01_read_stands_meta.csv", header=TRUE, sep=',')
head(df.meta)



# Getting otu_table - STANDS
df$OTU <- c(1:length(df$group))
rownames(df) <- paste0('group',df$OTU)
otu.table <- df %>% dplyr::select(-OTU, -group)
head(otu.table)
OTU = otu_table(otu.table, taxa_are_rows = TRUE)

df.tax <- df %>% dplyr::select(group)
colnames(df.tax) <- c("Phylum")
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


# Convert to tab

# Transform
pseq.stands.rel <- microbiome::transform(pseq.stands, 'clr')
pseq.stands.rel@otu_table


# Convert to dataframe
dtab.stands <- as.data.frame(pseq.stands.rel@otu_table)
class(dtab.stands)
head(dtab.stands)

# Gather
dtab.stands$group <- rownames(dtab.stands)
dg.stands <- dtab.stands %>% gather(key = "sample", value = "CLR abd",-group)
dg.stands

# Add meta
#df.meta.format <- df.meta %>% dplyr::rename(long.name = long.name.sequencing_R.,
#                                            Country = Country_x)
#df.meta.format <- df.meta.format %>% dplyr::select(long.name, short.name, site, Cluster, Country, Lattitude, Longitude) %>% distinct()
dg.meta.stands <- merge(dg.stands, df.meta.format, by.x = "sample", by.y = "long.name", sort=FALSE)
head(dg.meta.stands)

# Add taxonomy
df.taxonomy <- as.data.frame(pseq.stands.rel@tax_table)
df.taxonomy$group <- rownames(df.taxonomy)

dm.stands <- merge(dg.meta.stands, df.taxonomy, by = 'group', sort=FALSE)
head(dm.stands)






#####################
#-----  PLOTS ------#
#####################



guilds <- c("Ectomycorrhizal", "Undefined Saprotroph", "Ericoid Mycorrhizal",
            "Endophyte", "Plant Pathogen", "Wood Saprotroph","Animal Pathogen","Fungal Parasite")
dsub.stands <- dm.stands %>% filter(Phylum %in% guilds) 
dsub.stands %>% head()

unique(dsub.stands$Phylum)
dsub.stands$Label <- factor(dsub.stands$Phylum, levels=c("Ectomycorrhizal","Undefined Saprotroph","Ericoid Mycorrhizal","Endophyte",
                                                         "Plant Pathogen","Wood Saprotroph","Animal Pathogen","Fungal Parasite"))

dsub.stands$region1 <- ifelse(dsub.stands$site %in% c("AMOS","STFE"),"Boreal","Rest")
dsub.stands$region2 <- ifelse(dsub.stands$site %in% c("STET","ESSI","FORE"),"Cold_temperate",dsub.stands$region1)
dsub.stands$region <- ifelse(dsub.stands$site %in% c("Santiago","FLOR1","FP2"),"Warm_temperate",dsub.stands$region2)
dsub.stands$region <- factor(dsub.stands$region, levels=c("Boreal","Cold_temperate","Warm_temperate"))
head(dsub.stands)


write.table(dsub.stands, file = "guilds_04_clr_EMF_inputs.tsv", sep="\t", row.names = F, col.names = T, append=F, quote=F)



#####################
#-----   SITE ------#
#####################
# Plotting average abundance for each site

head(dsub.stands)
dsub.stands$site <- factor(dsub.stands$site, levels = c("AMOS","STFE","STET","ESSI","FORE","Santiago","FLOR1","FP2"))
dsub.stands$Country <- factor(dsub.stands$Country, levels = c("Canada","Mexico"))

p1 <- ggplot(dsub.stands) +
  geom_boxplot(data=dsub.stands, aes(x = site, y = `CLR abd`, fill = Country), outlier.shape = NA) +
  geom_point(data=dsub.stands, aes(x = site, y = `CLR abd`), color="grey40", pch=21, size=1) +
  facet_wrap(~Phylum)
p1




#####################
#-----  TESTS ------#
#####################



###############
### TESTS of difference between countries

head(dsub.stands)

# Transform CLR so it's positive

outputs <- list()
for (fam in unique(dsub.stands$Label)) {
  df_ <- dsub.stands[dsub.stands$Label == fam,]
  min_clr <- min(df_$`CLR abd`)
  df_$CLR_abd <- df_$`CLR abd` + abs(min_clr) + 0.01
  df.models <- run_nested_models(df_, "CLR_abd", "Country", "site")
  df.models$family <- fam
  outputs[[fam]] <- df.models
}

outputs[[2]]


dOUT <- bind_rows(outputs)
head(dOUT)
dOUT %>% filter(MOD == "GLS" & EVAL == "OK")

write.table(dOUT, file = "guilds_04_clr_EMF_country_tests.tsv", sep="\t", row.names = F, col.names = T, append=F, quote=F)


###############
### TESTS of difference between REGIONS

outputs <- list()
for (fam in unique(dsub.stands$Label)) {
  df_ <- dsub.stands[dsub.stands$Label == fam,]
  min_clr <- min(df_$`CLR abd`)
  df_$CLR_abd <- df_$`CLR abd` + abs(min_clr) + 0.01
  df.models <- run_nested_models(df_, "CLR_abd", "region", "site")
  df.models$family <- fam
  outputs[[fam]] <- df.models
}

outputs[[2]]


dOUT <- bind_rows(outputs)
head(dOUT)
dOUT %>% filter(MOD == "GLS" & EVAL == "OK")

write.table(dOUT, file = "guilds_04_clr_EMF_region_tests.tsv", sep="\t", row.names = F, col.names = T, append=F, quote=F)


