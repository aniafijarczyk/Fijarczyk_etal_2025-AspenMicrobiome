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

# Raw rarefied abundances
load("../01_data_prep/AVS_filtered_phyloseq_rarefied.RData")
sample_data(pseq.16s) %>% head()
otu_table(pseq.16s)
pseq.16s

tab_blue <- "#1f77b4"
tab_orange <- "#ff7f0e"
ametrics <- c('diversity_shannon', 'chao1', 'diversity_gini_simpson','dominance_core_abundance','evenness_pielou')
site_order <- c("AB","FORE","FP","SP")



#============================================#
#             Alpha diversity                #
#============================================#


# Metadata
meta <- data.frame(sample_data(pseq.16s))
meta %>% head()
dim(meta)

# Alpha diversity
alpha <- microbiome::alpha(pseq.16s, index = "all")
# The core_abundance function refers to the relative proportion of the core species
alpha$core_abundance <- core_abundance(pseq.16s, detection = .1/100, prevalence = 50/100)
# The coverage index gives the number of groups needed to have a given proportion of the ecosystem occupied (by default 0.5 ie 50%).
alpha$coverage <- coverage(pseq.16s, threshold = 0.5)
head(alpha)
dim(alpha)
colnames(alpha)

# Merging
meta.alpha <- merge(meta, alpha, by="row.names", sort=FALSE)
meta.alpha %>% head()

write.table(meta.alpha, file="01_alpha_16S.tsv", sep='\t', row.names = FALSE, col.names = TRUE, quote = FALSE, append=FALSE)






################
### Overview ###
################

head(meta.alpha)


meta.alpha_ <- meta.alpha %>% dplyr::rename(
  Country = Population,
  Site = Cluster
)
head(meta.alpha_)


meta.alpha.soil <- meta.alpha_ %>% filter(Sample.type == 'bulk soil')
meta.alpha.soil


ametrics <- c('diversity_shannon', 'chao1', 'diversity_gini_simpson','dominance_core_abundance','evenness_pielou')
dg.alpha.soil <- meta.alpha.soil %>% gather(key = 'measure', value = 'value', all_of(ametrics))
head(dg.alpha.soil)

site_order <- c("AB","FORE","FP","SP")
dg.alpha.soil$measure <- factor(dg.alpha.soil$measure, levels = ametrics)
dg.alpha.soil$Site <- factor(dg.alpha.soil$Site, levels = site_order)


p5 <- ggplot(dg.alpha.soil) + aes(y = value, x = Site, fill = Country) +
  geom_boxplot() +
  facet_wrap(~measure, scales = "free_y", ncol = 3)
p5







#####################
#-----  TESTS ------#
#####################
results <- list()


################# Soil

head(meta.alpha)
meta.alpha$Sample.type
dalpha <- meta.alpha %>% filter(Sample.type == "bulk soil") %>% dplyr::select(SampleID,Genotype,Sample.type,Population,Cluster,
                                                                              chao1,diversity_gini_simpson)
dim(dalpha)
head(dalpha)


### Gini simpson

dalpha$diversity_gini_simpson
dalpha$log_gini <- log(dalpha$diversity_gini_simpson + 0.01)
dalpha$log_gini
min_stat <- min(dalpha$log_gini)
dalpha$log_gini_trans <- dalpha$log_gini + abs(min_stat) + 0.01
dalpha$log_gini_trans
hist(dalpha$diversity_gini_simpson)
hist(dalpha$log_gini_trans)
plot(x = dalpha$log_gini, y = dalpha$log_gini_trans)


df.models <- run_nested_models(dalpha, "log_gini_trans", "Population", "Cluster")
df.models
df.models$marker <- "16S"
df.models$metric <- "log_gini"
df.models$source <- "soil"
results[[1]] <- df.models

#write.table(df.models, file = "01_alpha_16S_soil_tests.tsv", sep="\t", row.names = F, col.names = T, append=F, quote=F)


### chao1
dalpha$chao1
df.models <- run_nested_models(dalpha, "chao1", "Population", "Cluster")
df.models
df.models$marker <- "16S"
df.models$metric <- "chao1"
df.models$source <- "soil"
results[[2]] <- df.models

#write.table(df.models, file = "01_chao1_16S_soil_tests.tsv", sep="\t", row.names = F, col.names = T, append=F, quote=F)







############# Rhizo

dalpha <- meta.alpha %>% filter(Sample.type == "rhizosphere") %>% dplyr::select(SampleID,Genotype,Sample.type,Population,Cluster,
                                                                                chao1,diversity_gini_simpson)
dim(dalpha)
head(dalpha)

### Gini simpson
dalpha$diversity_gini_simpson
dalpha$log_gini <- log(dalpha$diversity_gini_simpson + 0.01)
dalpha$log_gini
min_stat <- min(dalpha$log_gini)
dalpha$log_gini_trans <- dalpha$log_gini + abs(min_stat) + 0.01
dalpha$log_gini_trans
hist(dalpha$diversity_gini_simpson)
hist(dalpha$log_gini_trans)
plot(x = dalpha$log_gini, y = dalpha$log_gini_trans)

df.models <- run_nested_models(dalpha, "log_gini_trans", "Population", "Cluster")
df.models
df.models$marker <- "16S"
df.models$metric <- "log_gini"
df.models$source <- "rhizo"
results[[3]] <- df.models

#write.table(df.models, file = "01_alpha_16S_rhizo_tests.tsv", sep="\t", row.names = F, col.names = T, append=F, quote=F)


### chao1

dalpha$chao1
df.models <- run_nested_models(dalpha, "chao1", "Population", "Cluster")
df.models
df.models$marker <- "16S"
df.models$metric <- "chao1"
df.models$source <- "rhizo"
results[[4]] <- df.models

#write.table(df.models, file = "01_chao1_16S_rhizo_tests.tsv", sep="\t", row.names = F, col.names = T, append=F, quote=F)



############# Roots

dalpha <- meta.alpha %>% filter(Sample.type == "roots" & Population %in% c("NENA","MX")) %>% dplyr::select(SampleID,Genotype,Sample.type,Population,Cluster,
                                                                                chao1,diversity_gini_simpson)
dim(dalpha)
head(dalpha)

### Gini simpson
dalpha$diversity_gini_simpson
dalpha$log_gini <- log(dalpha$diversity_gini_simpson + 0.01)
dalpha$log_gini
min_stat <- min(dalpha$log_gini)
dalpha$log_gini_trans <- dalpha$log_gini + abs(min_stat) + 0.01
dalpha$log_gini_trans
hist(dalpha$diversity_gini_simpson)
hist(dalpha$log_gini_trans)
plot(x = dalpha$log_gini, y = dalpha$log_gini_trans)

df.models <- run_nested_models(dalpha, "log_gini_trans", "Population", "Cluster")
df.models
df.models$marker <- "16S"
df.models$metric <- "log_gini"
df.models$source <- "roots"
results[[5]] <- df.models

#write.table(df.models, file = "01_alpha_16S_roots_tests.tsv", sep="\t", row.names = F, col.names = T, append=F, quote=F)


### chao1

dalpha$chao1
df.models <- run_nested_models(dalpha, "chao1", "Population", "Cluster")
df.models
df.models$marker <- "16S"
df.models$metric <- "chao1"
df.models$source <- "roots"
results[[6]] <- df.models

#write.table(df.models, file = "01_chao1_16S_roots_tests.tsv", sep="\t", row.names = F, col.names = T, append=F, quote=F)




###
dres <- bind_rows(results)
dres

write.table(dres, file = "01_alpha_16S_ALL_TEST_RESULTS.tsv", sep="\t", row.names = F, col.names = T, append=F, quote=F)



### END

