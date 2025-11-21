rm(list = ls())


library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(phyloseq)
library(microbiome)
library(readxl)
library(stringr)
#library(rstatix) # Welsh anova

library(MASS)
library(nlme)
library(emmeans) # to get Tukey test on gls model

source("run_NESTED_ANOVA.R")


#########################
#--------DATA-----------#
#########################


df <- read.csv("./guilds_02_read_greenhouse_roots.csv", header=TRUE, sep=',')
head(df)
df.meta <- read.csv("./guilds_02_read_greenhouse_meta_roots.csv", header=TRUE, sep=',')
head(df.meta)



rownames(df) <- df$group
otu.table <- df %>% dplyr::select(-group)
head(otu.table)
OTU = otu_table(otu.table, taxa_are_rows = TRUE)

df.tax <- df %>% dplyr::select(group)
colnames(df.tax) <- c("Family")
df.tax.mat <- df.tax %>% as.matrix()
TAX = tax_table(df.tax.mat)
TAX


head(df.meta)
df.meta.format <- df.meta %>% dplyr::select(sample, long.name, Cluster, Population, Sample.type) %>% distinct()
rownames(df.meta.format) <- df.meta.format$long.name
head(df.meta.format)
samdata <- sample_data(df.meta.format)

# Importing to phyloseq
pseq.stands <- phyloseq(OTU, TAX, sam_data = samdata)
pseq.stands




#####################
#---  TRANSFORM  ---#
#####################

df.meta.format$Sample.type
pseq.16s.soil <- subset_samples(pseq.stands, sample_data(pseq.stands)$Sample.type=="roots")
pseq.16s.soil

# Convert to tab

# Transform
pseq.stands.rel <- microbiome::transform(pseq.16s.soil, 'clr')
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
#df.meta.format <- df.meta %>% dplyr::rename(long.name = long.name.sequencing_R.,
#                                            Country = Country_x)
#df.meta.format <- df.meta.format %>% dplyr::select(long.name, short.name, site, Cluster, Country, Lattitude, Longitude) %>% distinct()
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


guilds <- c("Ectomycorrhizal", "Undefined Saprotroph", "Ericoid Mycorrhizal",
            "Endophyte", "Plant Pathogen", "Wood Saprotroph","Animal Pathogen","Fungal Parasite")
dsub.stands <- dm.stands %>% filter(family %in% guilds) 
dsub.stands %>% head()

unique(dsub.stands$Family)
dsub.stands$Label <- factor(dsub.stands$Family, levels=c("Ectomycorrhizal","Undefined Saprotroph","Ericoid Mycorrhizal","Endophyte",
                                                         "Plant Pathogen","Wood Saprotroph","Animal Pathogen","Fungal Parasite"))
#                            labels=c("Ecto-\nmycorrhizal","Undefined\nsaprotroph","Ericoid\nmycorrhizal","\nEndophyte",
#                                     "Plant\npathogen","Wood\nsaprotroph","Animal\npathogen","Fungal\nparasite"))

head(dsub.stands)
write.table(dsub.stands, file = "guilds_03_clr_roots_inputs.tsv", sep="\t", row.names = F, col.names = T, append=F, quote=F)



#####################
#-----   SITE ------#
#####################

# Plotting average abundance for each site

head(dsub.stands)
dsub.stands$Cluster
dsub.stands$site <- factor(dsub.stands$Cluster, levels = c("AB","FORE","FP","SP"))
dsub.stands$Population <- factor(dsub.stands$Population, levels = c("NENA","MX"))

p1 <- ggplot(dsub.stands) +
  geom_boxplot(data=dsub.stands, aes(x = Cluster, y = `CLR abd`, fill = Population), outlier.shape = NA) +
  geom_point(data=dsub.stands, aes(x = Cluster, y = `CLR abd`), color="grey40", pch=21, size=1) +
  facet_grid(~Label)

p1





#####################
#-----  TESTS  -----#
#####################


###############
### TESTS of difference between countries

head(dsub.stands)

# Transform CLR so it's positive

outputs <- list()
for (fam in dsub.stands$Label) {
  df_ <- dsub.stands[dsub.stands$Label == fam,]
  min_clr <- min(df_$`CLR abd`)
  df_$CLR_abd <- df_$`CLR abd` + abs(min_clr) + 0.01
  df.models <- run_nested_models(df_, "CLR_abd", "Population", "Cluster")
  df.models$family <- fam
  outputs[[fam]] <- df.models
}

outputs[[2]]


dOUT <- bind_rows(outputs)
head(dOUT)



write.table(dOUT, file = "guilds_03_clr_roots_tests.tsv", sep="\t", row.names = F, col.names = T, append=F, quote=F)


