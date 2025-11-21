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


source("run_NESTED_ANOVA.R")

#########################
#--------DATA-----------#
#########################

# Raw rarefied abundances, with a few outlier samples removed
load("../../2022_MicrobiomeAspenGreenhouse/01_data_prep/AVS_filtered_phyloseq_rarefied.RData")
sample_data(pseq.its) %>% head()
otu_table(pseq.its)
pseq.its

tab_blue <- "#1f77b4"
tab_orange <- "#ff7f0e"
ametrics <- c('diversity_shannon', 'chao1', 'diversity_gini_simpson','dominance_core_abundance','evenness_pielou')
site_order <- c("AB","FORE","FP","SP")


#####################################################
#--------Subsetting only ecto-mycorrhiza -----------#
#####################################################

funfile <- "../../2022_MicrobiomeAspenGreenhouse/DATA/PeuplierGH_jan2024/ITS/all/PeuplierGH_ITS_jan2024.ASV_table_rarefied_13594_dnNA.guilds.txt"
df <- read.csv(funfile, sep="\t", header=TRUE, comment.char = "")
head(df)

# Select ectomycorrhizal
df.ecto <- df %>% filter(Confidence.Ranking %in% c("Probable","Highly Probable")) %>% filter(Guild == "Ectomycorrhizal")
otu.ecto <- df.ecto %>% pull(X.OTU.ID)

# Subset otu
ecto_subset <- subset(otu_table(pseq.its), rownames(otu_table(pseq.its)) %in% otu.ecto)
ecto.its <- merge_phyloseq(ecto_subset, tax_table(pseq.its), sample_data(pseq.its))
ecto.its

sample_data(ecto.its)
tax_table(ecto.its)

#============================================#
#             Alpha diversity                #
#============================================#


# Metadata
meta <- data.frame(sample_data(ecto.its))
meta %>% head()
dim(meta)

alpha.1 <- microbiome::alpha(ecto.its, index = c("observed","diversity_shannon","chao1",'dominance_core_abundance',"diversity_gini_simpson","evenness_pielou"))
alpha.1
# The coverage index gives the number of groups needed to have a given proportion of the ecosystem occupied (by default 0.5 ie 50%).
#alpha.1$coverage <- coverage(ecto.its, threshold = 0.5)
head(alpha.1)
dim(alpha.1)
colnames(alpha.1)

# Merging
meta.alpha <- merge(meta, alpha.1, by="row.names", sort=FALSE)
meta.alpha <- meta.alpha %>% filter(Sample.type == "bulk soil")
# Chanching column names
#meta.alpha <- meta.alpha %>% rename(Country = Population,
#                                    Site = Cluster)
head(meta.alpha)

meta.alpha$diversity_gini_simpson[is.na(meta.alpha$diversity_gini_simpson)] <- 0
write.table(meta.alpha, file="guilds_04_alpha_soil.tsv", sep='\t', row.names = FALSE, col.names = TRUE, quote = FALSE, append=FALSE)






################
### BOXPLOTS ###
################

head(meta.alpha)


ametrics <- c('diversity_shannon', 'chao1', 'diversity_gini_simpson','dominance_core_abundance','evenness_pielou')
dg.alpha.soil <- meta.alpha %>% gather(key = 'measure', value = 'value', all_of(ametrics))
head(dg.alpha.soil)


site_order <- c("AB","FORE","FP","SP")
dg.alpha.soil$measure <- factor(dg.alpha.soil$measure, levels = ametrics)
dg.alpha.soil$Cluster <- factor(dg.alpha.soil$Cluster, levels = site_order)


p5 <- ggplot(dg.alpha.soil) + aes(y = value, x = Cluster, fill = Population) +
  geom_boxplot() +
  facet_wrap(~measure, scales = "free_y", ncol = 3)
p5



#####################
#-----  TESTS  -----#
#####################


###############
### TESTS of difference between countries

head(meta.alpha)

meta.alpha$diversity_gini_simpson
meta.alpha$log_gini <- log(meta.alpha$diversity_gini_simpson + 0.01)
meta.alpha$log_gini
min_stat <- min(meta.alpha$log_gini)
meta.alpha$log_gini_trans <- meta.alpha$log_gini + abs(min_stat) + 0.01
meta.alpha$log_gini_trans
hist(meta.alpha$diversity_gini_simpson)
hist(meta.alpha$log_gini)
plot(x = meta.alpha$diversity_gini_simpson, y = meta.alpha$log_gini_trans)



df.models <- run_nested_models(meta.alpha, "log_gini_trans", "Population", "Cluster")
df.models

write.table(df.models, file = "guilds_04_alpha_soil_tests.tsv", sep="\t", row.names = F, col.names = T, append=F, quote=F)






###############
### TESTS of difference between countries

head(meta.alpha)
meta.alpha$chao1_trans <- meta.alpha$chao1 +1


df.models <- run_nested_models(meta.alpha, "chao1_trans", "Population", "Cluster")
df.models

write.table(df.models, file = "guilds_04_chao1_soil_tests.tsv", sep="\t", row.names = F, col.names = T, append=F, quote=F)
















