rm(list=ls())



library(openxlsx)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(chemometrics)
library(MASS)
library(DAAG)
library(phyloseq)
library(microbiome)
library(cowplot)
library(nlme)


source("run_NESTED_ANOVA_EXT.R")

#########################
#--------DATA-----------#
#########################



# Raw rarefied abundances, with a few outlier samples removed
load("../../2022_MicrobiomeQuebecMexico/13_data_prep_18s/AVS_phyloseq_rarefied.RData")
sample_data(pseq.18s.AM) %>% head()
otu_table(pseq.18s.AM)
pseq.18s.AM

tab_blue <- "#1f77b4"
tab_orange <- "#ff7f0e"
ametrics <- c('diversity_shannon', 'chao1', 'diversity_gini_simpson', 'dominance_core_abundance','evenness_pielou')
site_order <- c("AMOS","STFE","STET","ESSI","FORE","Santiago","FLOR1","FP2")


#============================================#
#             Alpha diversity                #
#============================================#


# Metadata
meta <- data.frame(sample_data(pseq.18s.AM))
meta %>% head()
dim(meta)

# Alpha diversity
alpha.1 <- microbiome::alpha(pseq.18s.AM, index = c("diversity_shannon","chao1","diversity_gini_simpson","evenness_pielou"))
alpha.1

# The core_abundance function refers to the relative proportion of the core species
#alpha.1$core_abundance <- core_abundance(pseq.18s.AM, detection = .1/100, prevalence = 50/100)
# The coverage index gives the number of groups needed to have a given proportion of the ecosystem occupied (by default 0.5 ie 50%).
alpha.1$coverage <- coverage(pseq.18s.AM, threshold = 0.5)
head(alpha.1)
dim(alpha.1)
colnames(alpha.1)

# Merging
meta.alpha <- merge(meta, alpha.1, by="row.names", sort=FALSE)
meta.alpha %>% head()

# Chanching column names
meta.alpha <- meta.alpha %>% rename(Country = cluster,
                                    Site = site)
head(meta.alpha)


dalpha <- meta.alpha %>% filter(source == 'Soil') %>% filter(location == "natural_stands")
dalpha <- dalpha %>% rename(cluster = Country, site = Site)
dalpha$diversity_gini_simpson[is.na(dalpha$diversity_gini_simpson)] <- 0
dim(dalpha)
head(dalpha)


dalpha$region1 <- ifelse(dalpha$site %in% c("AMOS","STFE"),"Boreal","Rest")
dalpha$region2 <- ifelse(dalpha$site %in% c("STET","ESSI","FORE"),"Cold_temperate",dalpha$region1)
dalpha$region <- ifelse(dalpha$site %in% c("Santiago","FLOR1","FP2"),"Warm_temperate",dalpha$region2)
dalpha$region <- factor(dalpha$region, levels=c("Boreal","Cold_temperate","Warm_temperate"))
head(dalpha)


write.table(dalpha, file = "guilds_05_alpha_AMF_inputs.tsv", sep="\t", row.names = F, col.names = T, append=F, quote=F)



################
### BOXPLOTS ###
################

head(dalpha)
ametrics <- c('diversity_shannon', 'chao1', 'diversity_gini_simpson')
dg.alpha.soil <- dalpha %>% gather(key = 'measure', value = 'value', all_of(ametrics))
head(dg.alpha.soil)

# Basic plot
site_order <- c("AMOS","STFE","STET","ESSI","FORE","Santiago","FLOR1","FP2")
dg.alpha.soil$measure <- factor(dg.alpha.soil$measure, levels = ametrics)
dg.alpha.soil$site <- factor(dg.alpha.soil$site, levels = site_order)


p5 <- ggplot(dg.alpha.soil) + aes(y = value, x = site, fill = cluster) +
  geom_boxplot() +
  facet_wrap(~measure, scales = "free_y", ncol = 4)
p5








#####################
#-----  TESTS  -----#
#####################


###########################################################################
### TESTS of difference between countries

#dalpha <- dalpha %>% filter(diversity_gini_simpson>0)

head(dalpha)

dalpha$diversity_gini_simpson
dalpha$log_gini <- log(dalpha$diversity_gini_simpson + 0.01)
dalpha$log_gini
min_stat <- min(dalpha$log_gini)
dalpha$log_gini_trans <- dalpha$log_gini + abs(min_stat) + 0.01
dalpha$log_gini_trans
hist(dalpha$diversity_gini_simpson)
hist(dalpha$log_gini)
plot(x = dalpha$diversity_gini_simpson, y = dalpha$log_gini_trans)

# Testing raw model
df.models <- run_nested_models(dalpha, "log_gini_trans", "cluster", "site")
df.models

#write.table(df.models, file = "guilds_05_alpha_AMF_country_tests.tsv", sep="\t", row.names = F, col.names = T, append=F, quote=F)

####################
# Removing samples with gini index = 0
dalpha %>% filter(diversity_gini_simpson == 0)
dalpha_1 <- dalpha %>% filter(diversity_gini_simpson > 0)
dim(dalpha)
dim(dalpha_1)

df.models <- run_nested_models(dalpha_1, "log_gini_trans", "cluster", "site")
df.models


write.table(df.models, file = "guilds_05_alpha_AMF_country_tests.tsv", sep="\t", row.names = F, col.names = T, append=F, quote=F)


###########################################################################
### TESTS of difference between regions

df.models <- run_nested_models(dalpha_1, "log_gini_trans", "region", "site")
df.models


write.table(df.models, file = "guilds_05_alpha_AMF_region_tests.tsv", sep="\t", row.names = F, col.names = T, append=F, quote=F)







###########################################################################
###############
### TESTS of difference CHAO1

head(dalpha)
dalpha$chao1
dalpha$chao1_trans <- dalpha$chao1 + 1



df.models <- run_nested_models(dalpha, "chao1_trans", "cluster", "site")
df.models

write.table(df.models, file = "guilds_05_chao1_country_tests.tsv", sep="\t", row.names = F, col.names = T, append=F, quote=F)



## Region

df.models <- run_nested_models(dalpha, "chao1_trans", "region", "site")
df.models

write.table(df.models, file = "guilds_05_chao1_region_tests.tsv", sep="\t", row.names = F, col.names = T, append=F, quote=F)





