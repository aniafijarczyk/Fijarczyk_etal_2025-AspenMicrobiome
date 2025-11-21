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
sample_data(rseq.its) %>% head()
otu_table(rseq.its)
rseq.its

tab_blue <- "#1f77b4"
tab_orange <- "#ff7f0e"
ametrics <- c('diversity_shannon', 'chao1', 'dominance_core_abundance','evenness_pielou')
site_order <- c("AMOS","STFE","STET","ESSI","FORE","Santiago","FLOR1","FP2")


#####################################################
#--------Subsetting only ecto-mycorrhiza -----------#
#####################################################

funfile <- "D:/NRCan/2022_MicrobiomeQuebecMexico/DATA/ITS/Full/Peuplier_REPRISE_ITS_oct2022.ASV_table_rarefied_11982_dnNA.guilds.txt"
df <- read.csv(funfile, sep="\t", header=TRUE, comment.char = "")
head(df)

# Select ectomycorrhizal
df.ecto <- df %>% filter(Confidence.Ranking %in% c("Probable","Highly Probable")) %>% filter(Guild == "Ectomycorrhizal")
otu.ecto <- df.ecto %>% pull(X.OTU.ID)

# Subset otu
ecto_subset <- subset(otu_table(rseq.its), rownames(otu_table(rseq.its)) %in% otu.ecto)
ecto.its <- merge_phyloseq(ecto_subset, tax_table(rseq.its), sample_data(rseq.its))
ecto.its


#============================================#
#             Alpha diversity                #
#============================================#


# Metadata
meta <- data.frame(sample_data(ecto.its))
meta %>% head()
dim(meta)

# Alpha diversity
alpha <- microbiome::alpha(ecto.its, index = "all")
# The core_abundance function refers to the relative proportion of the core species
alpha$core_abundance <- core_abundance(ecto.its, detection = .1/100, prevalence = 50/100)
# The coverage index gives the number of groups needed to have a given proportion of the ecosystem occupied (by default 0.5 ie 50%).
alpha$coverage <- coverage(ecto.its, threshold = 0.5)
head(alpha)
dim(alpha)
colnames(alpha)

# Merging
meta.alpha <- merge(meta, alpha, by="row.names", sort=FALSE)
meta.alpha %>% head()
dalpha <- meta.alpha %>% filter(source == 'Soil')

#write.table(dalpha, file="guilds_05_alpha_EMF.tsv", sep='\t', row.names = FALSE, col.names = TRUE, quote = FALSE, append=FALSE)



################
### BOXPLOTS ###
################

head(meta.alpha)

ametrics <- c('diversity_shannon', 'chao1', 'diversity_gini_simpson')
dg.alpha.soil <- dalpha %>% gather(key = 'measure', value = 'value', all_of(ametrics))
head(dg.alpha.soil)


site_order <- c("AMOS","STFE","STET","ESSI","FORE","Santiago","FLOR1","FP2")
dg.alpha.soil$measure <- factor(dg.alpha.soil$measure, levels = ametrics)
dg.alpha.soil$site <- factor(dg.alpha.soil$site, levels = site_order)


p5 <- ggplot(dg.alpha.soil) + aes(y = value, x = site, fill = Country) +
  geom_boxplot() +
  facet_wrap(~measure, scales = "free_y", ncol = 4)
p5








#####################
#-----  PLOTS ------#
#####################


head(dalpha)
dalpha$source


dalpha$region1 <- ifelse(dalpha$site %in% c("AMOS","STFE"),"Boreal","Rest")
dalpha$region2 <- ifelse(dalpha$site %in% c("STET","ESSI","FORE"),"Cold_temperate",dalpha$region1)
dalpha$region <- ifelse(dalpha$site %in% c("Santiago","FLOR1","FP2"),"Warm_temperate",dalpha$region2)
dalpha$region <- factor(dalpha$region, levels=c("Boreal","Cold_temperate","Warm_temperate"))
head(dalpha)


write.table(dalpha, file = "guilds_05_alpha_EMF_inputs.tsv", sep="\t", row.names = F, col.names = T, append=F, quote=F)








#####################
#-----  TESTS ------#
#####################


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


#############################################################
### Country

df.models <- run_nested_models(dalpha, "log_gini_trans", "Country", "site")
df.models

write.table(df.models, file = "guilds_05_alpha_country_tests.tsv", sep="\t", row.names = F, col.names = T, append=F, quote=F)


#############################################################
### Region
df.models <- run_nested_models(dalpha, "log_gini_trans", "region", "site")
df.models



write.table(df.models, file = "guilds_05_alpha_region_tests.tsv", sep="\t", row.names = F, col.names = T, append=F, quote=F)


############ CHAO1

head(dalpha)

dalpha$chao1


### Country
df.models <- run_nested_models(dalpha, "chao1", "Country", "site")
df.models

write.table(df.models, file = "guilds_05_chao1_country_tests.tsv", sep="\t", row.names = F, col.names = T, append=F, quote=F)


### Region
df.models <- run_nested_models(dalpha, "chao1", "region", "site")
df.models

write.table(df.models, file = "guilds_05_chao1_region_tests.tsv", sep="\t", row.names = F, col.names = T, append=F, quote=F)






