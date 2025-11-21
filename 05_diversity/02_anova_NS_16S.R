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



site_order <- c("AMOS","STFE","STET","ESSI","FORE","Santiago","FLOR1","FP2")

meta.alpha <- read.csv("03_alpha_stats.tsv", sep='\t', header=TRUE)



################
### Overview ###
################


ametrics <- c('diversity_shannon', 'diversity_gini_simpson', 'chao1', 'dominance_core_abundance','evenness_pielou')
meta.alpha.soil <- meta.alpha %>% filter(source == 'Soil' & marker == '16')
dg.alpha.soil <- meta.alpha.soil %>% gather(key = 'measure', value = 'value', all_of(ametrics))
site_order <- c("AMOS","STFE","STET","ESSI","FORE","Santiago","FLOR1","FP2")
meta.alpha.soil$site <- factor(meta.alpha.soil$site, levels = site_order)

#c(16:39)
dg.alpha.soil <- meta.alpha.soil %>% gather(key = 'measure', value = 'value', c(16:39))
dg.alpha.soil.1 <- meta.alpha.soil %>% gather(key = 'measure', value = 'value', c(16:21))
dg.alpha.soil.2 <- meta.alpha.soil %>% gather(key = 'measure', value = 'value', c(22:27))
dg.alpha.soil.3 <- meta.alpha.soil %>% gather(key = 'measure', value = 'value', c(28:33))
dg.alpha.soil.4 <- meta.alpha.soil %>% gather(key = 'measure', value = 'value', c(34:38))

dg.alpha.soil.1 %>% head()



p1 <- ggplot(dg.alpha.soil.1) + aes(x = value, y = site, fill = Country) +
  geom_boxplot() +
  facet_grid(~measure, scales = "free_x")
  
p2 <- ggplot(dg.alpha.soil.2) + aes(x = value, y = site, fill = Country) +
  geom_boxplot() +
  facet_grid(~measure, scales = "free_x")
  
p3 <- ggplot(dg.alpha.soil.3) + aes(x = value, y = site, fill = Country) +
  geom_boxplot() +
  facet_grid(~measure, scales = "free_x")

p4 <- ggplot(dg.alpha.soil.4) + aes(x = value, y = site, fill = Country) +
  geom_boxplot() +
  facet_grid(~measure, scales = "free_x")


plot_grid(p1,p2,p3,p4, nrow = 4)



colnames(dg.alpha.soil)
ametrics <- c('diversity_shannon', 'diversity_gini_simpson', 'chao1', 'dominance_core_abundance','evenness_pielou')
dg.alpha.soil %>% head()
dg.alpha.soil.sub <- dg.alpha.soil %>% filter(measure %in% ametrics)
dg.alpha.soil.sub %>% head()

dg.alpha.soil.sub$measure <- factor(dg.alpha.soil.sub$measure, levels = ametrics)
dg.alpha.soil.sub$site <- factor(dg.alpha.soil.sub$site, levels = site_order)
dg.alpha.soil.sub$value <- as.numeric(dg.alpha.soil.sub$value)

p5 <- ggplot(dg.alpha.soil.sub) + aes(y = value, x = site, fill = Country) +
  geom_boxplot() +
  facet_wrap(~measure, scales = "free_y", ncol = 5)
p5







#####################
#-----  TESTS ------#
#####################

head(meta.alpha)

dalpha <- meta.alpha %>% filter(marker == "16" & source == "Soil")
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


df.models <- run_nested_models(dalpha, "log_gini_trans", "Country", "site")
df.models

write.table(df.models, file = "03_alpha_16S_country_tests.tsv", sep="\t", row.names = F, col.names = T, append=F, quote=F)


### Region
df.models <- run_nested_models(dalpha, "log_gini_trans", "region", "site")
df.models

write.table(df.models, file = "03_alpha_16S_region_tests.tsv", sep="\t", row.names = F, col.names = T, append=F, quote=F)






### chao1

dalpha$chao1


df.models <- run_nested_models(dalpha, "chao1", "Country", "site")
df.models

write.table(df.models, file = "03_chao1_16S_country_tests.tsv", sep="\t", row.names = F, col.names = T, append=F, quote=F)


### Region
df.models <- run_nested_models(dalpha, "chao1", "region", "site")
df.models

write.table(df.models, file = "03_chao1_16S_region_tests.tsv", sep="\t", row.names = F, col.names = T, append=F, quote=F)










### END






