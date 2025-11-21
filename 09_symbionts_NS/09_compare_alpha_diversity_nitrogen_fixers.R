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
load("../01_data_prep/AVS_filtered_phyloseq_rarefied.RData")
sample_data(rseq.16s) %>% head()
otu_table(rseq.16s)
rseq.16s

tab_blue <- "#1f77b4"
tab_orange <- "#ff7f0e"
ametrics <- c('diversity_shannon', 'chao1', 'dominance_core_abundance','evenness_pielou')
site_order <- c("AMOS","STFE","STET","ESSI","FORE","Santiago","FLOR1","FP2")

### Subsetting samples
dmeta <- read.csv("../01_data_prep/dataPrep_01_overview_16s_META.tsv", sep="\t", header = T)
head(dmeta)
dmeta.soil <- dmeta %>% dplyr::filter(source == "Soil")
dim(dmeta.soil)
head(dmeta.soil)

head(sample_data(rseq.16s))
rseq <- subset_samples(rseq.16s, long.name %in% dmeta.soil$long.name)
rseq.16s
rseq


#####################################################
#--------     Subsetting  NF ASVs        -----------#
#####################################################

asv <- read.csv("nitro_02_read_NF_ASVs.txt", sep="\t", header=T)
asv$ASV
length(asv$ASV)
rseq.fams <- prune_taxa(asv$ASV, rseq)
rseq.fams



#============================================#
#             Alpha diversity                #
#============================================#


metrics <- c('diversity_shannon', 'chao1', 'dominance_core_abundance','evenness_pielou', 'diversity_gini_simpson')


meta <- data.frame(sample_data(rseq.fams))
alpha <- microbiome::alpha(rseq.fams, index = metrics)
dalpha <- merge(meta, alpha, by="row.names", sort=FALSE)
family <- c("NF")
dalpha$Family <- family
head(dalpha)

dalpha %>% head()
dalpha$diversity_gini_simpson[is.na(dalpha$diversity_gini_simpson)] <- 0

dalpha$region1 <- ifelse(dalpha$site %in% c("AMOS","STFE"),"Boreal","Rest")
dalpha$region2 <- ifelse(dalpha$site %in% c("STET","ESSI","FORE"),"Cold_temperate",dalpha$region1)
dalpha$region <- ifelse(dalpha$site %in% c("Santiago","FLOR1","FP2"),"Warm_temperate",dalpha$region2)
dalpha$region <- factor(dalpha$region, levels=c("Boreal","Cold_temperate","Warm_temperate"))
head(dalpha)


write.table(dalpha, file="nitro_05_alpha_inputs.tsv", sep='\t', row.names = FALSE, col.names = TRUE, quote = FALSE, append=FALSE)






################
### BOXPLOTS ###
################

metrics
dg.dalpha <- dalpha %>% gather(key = 'measure', value = 'value', all_of(metrics))
head(dg.dalpha)
unique(dg.dalpha$Family)

site_order <- c("AMOS","STFE","STET","ESSI","FORE","Santiago","FLOR1","FP2")
dg.dalpha$measure <- factor(dg.dalpha$measure, levels = metrics)
dg.dalpha$site <- factor(dg.dalpha$site, levels = site_order)


p5 <- ggplot(dg.dalpha) + aes(y = value, x = site, fill = Country) +
  geom_boxplot() +
  facet_wrap(Family~measure, scales = "free_y", ncol = 5)
p5






#####################
#-----  TESTS  -----#
#####################


###############
### TESTS of difference between countries

head(dalpha)
dalpha$diversity_gini_simpson

dalpha_1 <- dalpha %>% filter(diversity_gini_simpson > 0)
dalpha_1$log_gini <- log(dalpha_1$diversity_gini_simpson + 0.01)
dalpha_1$log_gini
min_stat <- min(dalpha_1$log_gini)
dalpha_1$log_gini_trans <- dalpha_1$log_gini + abs(min_stat) + 0.01
dalpha_1$log_gini_trans
hist(dalpha_1$diversity_gini_simpson)
hist(dalpha_1$log_gini)
plot(x = dalpha_1$diversity_gini_simpson, y = dalpha_1$log_gini_trans)


dalpha_1$log_gini_trans
df.models <- run_nested_models(dalpha_1, "log_gini_trans", "Country", "site")
df.models

write.table(df.models, file = "nitro_05_alpha_country_tests.tsv", sep="\t", row.names = F, col.names = T, append=F, quote=F)




################################################################################
### region

df.models <- run_nested_models(dalpha_1, "log_gini_trans", "region", "site")
df.models



write.table(df.models, file = "nitro_05_alpha_region_tests.tsv", sep="\t", row.names = F, col.names = T, append=F, quote=F)



##########################################################################################
### TESTS of difference CHAO1

head(dalpha)
dalpha$chao1
dalpha$chao1_trans <- dalpha$chao1 + 1


df.models <- run_nested_models(dalpha, "chao1_trans", "Country", "site")
df.models

write.table(df.models, file = "nitro_05_chao1_country_tests.tsv", sep="\t", row.names = F, col.names = T, append=F, quote=F)


#############################################
## Region

df.models <- run_nested_models(dalpha, "chao1_trans", "region", "site")
df.models

write.table(df.models, file = "nitro_05_chao1_region_tests.tsv", sep="\t", row.names = F, col.names = T, append=F, quote=F)





