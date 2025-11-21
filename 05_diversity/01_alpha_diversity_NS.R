rm(list=ls())


library(openxlsx)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(phyloseq)
library(microbiome)
library(cowplot)


#########################
#--------DATA-----------#
#########################


# Raw rarefied abundances, with a few outlier samples removed
load("../01_data_prep/AVS_filtered_phyloseq_rarefied.RData")
sample_data(rseq.16s) %>% head()


site_order <- c("AMOS","STFE","STET","ESSI","FORE","Santiago","FLOR1","FP2")

#============================================#
#             Alpha diversity                #
#============================================#



###########
### 16S ###
###########

# Metadata
tab.16s <- data.frame(sample_data(rseq.16s))
tab.16s %>% head()
tab.16s %>% filter(source == "Soil") %>% dim()

dim(tab.16s)

# Alpha diversity
alpha.16s <- microbiome::alpha(rseq.16s, index = "all")
head(alpha.16s)
# The core_abundance function refers to the relative proportion of the core species
#alpha.16s$core_abundance <- core_abundance(rseq.16s, detection = .1/100, prevalence = 50/100)
# The coverage index gives the number of groups needed to have a given proportion of the ecosystem occupied (by default 0.5 ie 50%).
alpha.16s$coverage <- coverage(rseq.16s, threshold = 0.5)
head(alpha.16s)
dim(alpha.16s)
colnames(alpha.16s)

# Merging
meta.alpha.16s <- merge(tab.16s, alpha.16s, by="row.names", sort=FALSE)
meta.alpha.16s %>% head()
dim(meta.alpha.16s)

###########
### ITS ###
###########

# Metadata
tab.its <- data.frame(sample_data(rseq.its))
tab.its %>% head()
dim(tab.its)

# Alpha diversity
alpha.its <- microbiome::alpha(rseq.its, index = "all")
# The core_abundance function refers to the relative proportion of the core species
#alpha.its$core_abundance <- core_abundance(rseq.its, detection = .1/100, prevalence = 50/100)
# The coverage index gives the number of groups needed to have a given proportion of the ecosystem occupied (by default 0.5 ie 50%).
alpha.its$coverage <- coverage(rseq.its, threshold = 0.5)
head(alpha.its)
dim(alpha.its)
colnames(alpha.its)

# Merging
meta.alpha.its <- merge(tab.its, alpha.its, by="row.names", sort=FALSE)
meta.alpha.its %>% head()
dim(meta.alpha.its)

###############
### COMBINE ###
###############

meta.alpha <- rbind(meta.alpha.16s, meta.alpha.its)
dim(meta.alpha)
head(meta.alpha)

meta.alpha$region1 <- ifelse(meta.alpha$site %in% c("AMOS","STFE"),"Boreal","Rest")
meta.alpha$region2 <- ifelse(meta.alpha$site %in% c("STET","ESSI","FORE"),"Cold_temperate",meta.alpha$region1)
meta.alpha$region <- ifelse(meta.alpha$site %in% c("Santiago","FLOR1","FP2"),"Warm_temperate",meta.alpha$region2)
meta.alpha$region <- factor(meta.alpha$region, levels=c("Boreal","Cold_temperate","Warm_temperate"))
head(meta.alpha)

write.table(meta.alpha, file="03_alpha_stats.tsv", sep='\t', row.names = FALSE, col.names = TRUE, quote = FALSE, append=FALSE)


