rm(list = ls())

library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(phyloseq)
library(microbiome)
library(readxl)
library(stringr)
library(RColorBrewer)
library(plyr)
library(openxlsx)
library(ggrepel)
library(vegan)



#########################
#----- FUNCTIONS -------#
#########################


# Keep taxa present with min coverage of x or more in more than y fraction of the samples
filter_taxa_min_frac <- function(pseq.all, min.samp, frac.samp) {
  avs = genefilter_sample(pseq.all, filterfun_sample(function(x) x > min.samp), A=frac.samp*nsamples(pseq.all))
  pseq.new = prune_taxa(avs, pseq.all)
  return(pseq.new)
}




#########################
#--------DATA-----------#
#########################


# Abundances
load("../../2022_MicrobiomeAspenGreenhouse/01_data_prep/AVS_filtered_phyloseq.RData")
pseq.16s

# Removing controls
pseq.16s <- subset_samples(pseq.16s, sample_data(pseq.16s)$Population!="neg control extraction")
pseq.16s <- subset_samples(pseq.16s, sample_data(pseq.16s)$Population!="neg control PCR")
sample_names(pseq.16s)

pseq.its <- subset_samples(pseq.its, sample_data(pseq.its)$Population!="neg control extraction")
pseq.its <- subset_samples(pseq.its, sample_data(pseq.its)$Population!="neg control PCR")
sample_names(pseq.its)

pseq.p16s <- subset_samples(pseq.p16s, sample_data(pseq.p16s)$Population!="neg control extraction")
pseq.p16s <- subset_samples(pseq.p16s, sample_data(pseq.p16s)$Population!="neg control PCR")
pseq.p16s <- subset_samples(pseq.p16s, sample_data(pseq.p16s)$Population!="UNK")
sample_names(pseq.p16s)

pseq.pits <- subset_samples(pseq.pits, sample_data(pseq.pits)$Population!="neg control extraction")
pseq.pits <- subset_samples(pseq.pits, sample_data(pseq.pits)$Population!="neg control PCR")
pseq.pits <- subset_samples(pseq.pits, sample_data(pseq.pits)$Population!="UNK")
sample_names(pseq.pits)



dmeta.16s <- data.frame(sample_data(pseq.16s))
dmeta.16s %>% head()
dmeta.16s$long.name

dmeta.its <- data.frame(sample_data(pseq.its))
dmeta.its %>% head()
dmeta.its$long.name

dmeta.p16s <- data.frame(sample_data(pseq.p16s))
dmeta.p16s %>% head()

dmeta.pits <- data.frame(sample_data(pseq.pits))
dmeta.pits %>% head()



### 16S

pseq.16s.2 = filter_taxa_min_frac(pseq.16s, 1, 0.05)
pseq.16s.2.rel = microbiome::transform(pseq.16s.2, 'compositional')
tab_16s <- as.data.frame(t(pseq.16s.2.rel@otu_table))
tab_16s[1:5, 1:5]


### ITS

pseq.its.2 = filter_taxa_min_frac(pseq.its, 1, 0.05)
pseq.its.2.rel = microbiome::transform(pseq.its.2, 'compositional')
tab_its <- as.data.frame(t(pseq.its.2.rel@otu_table))
tab_its[1:5, 1:5]

### 16S - phyllo

pseq.p16s.2 = filter_taxa_min_frac(pseq.p16s, 1, 0.05)
pseq.p16s.2.rel = microbiome::transform(pseq.p16s.2, 'compositional')
tab_p16s <- as.data.frame(t(pseq.p16s.2.rel@otu_table))
tab_p16s[1:5, 1:5]


### ITS

pseq.pits.2 = filter_taxa_min_frac(pseq.pits, 1, 0.05)
pseq.pits.2.rel = microbiome::transform(pseq.16s.2, 'compositional')
tab_pits <- as.data.frame(t(pseq.pits.2.rel@otu_table))
tab_pits[1:5, 1:5]


### Colors
source_cols <- c("#0072B2","#D55E00","#F0E442")
country_cols <- c("#E84A5F","#DDCC77")
site_cols <- brewer.pal(9,"BrBG")
phylum_cols <- c("#332288","#117733","#44AA99","#88CCEE","#DDCC77","#CC6677","#AA4499","#882255")
phylum_cols_wong <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

tab_blue <- "#1f77b4"
tab_orange <- "#ff7f0e"






##########################
#---    PERMANOVA      --#
##########################

# Checking of rownames of two datasets are the same
rownames(tab_16s)
dmeta.16s %>% head()
dmeta.16s$long.name

### Distances
bc_dist <- vegan::vegdist(tab_16s, method = "bray")
bc_dist


get_distance <- function(table) {
  bc_dist <- vegan::vegdist(table, method = "bray")
  return(bc_dist)
}

run_permanova_source <- function(distances, metatable) {
  perm <- adonis2(distances ~ Sample.type, data=metatable, permutations=999)
  df <- as.data.frame(perm)
  df$Set <- "Source"
  df$part <- c("Model","Residual","Total")
  return(df)
}

run_permanova_group <- function(distances, metatable) {
  perm <- adonis2(distances ~ Population, data=metatable, permutations=999)
  df <- as.data.frame(perm)
  df$Set <- "Group"
  df$part <- c("Model","Residual","Total")
  return(df)
}

run_permanova_group_strata <- function(distances, metatable) {
  perm <- adonis2(distances ~ Population, data=metatable, permutations=999, strata = metatable$Sample.type)
  df <- as.data.frame(perm)
  df$Set <- "Group"
  df$part <- c("Model","Residual","Total")
  return(df)
}

run_permanova_site_strata <- function(distances, metatable) {
  perm <- adonis2(distances ~ Cluster, data=metatable, permutations=999, strata = metatable$Population)
  df <- as.data.frame(perm)
  df$Set <- "Site"
  df$part <- c("Model","Residual","Total")
  return(df)
}

run_permanova_site_groupsource <- function(distances, metatable) {
  grouping <- paste0(metatable$Population,metatable$Sample.type)
  perm <- adonis2(distances ~ Cluster, data=metatable, permutations=999, strata = grouping)
  df <- as.data.frame(perm)
  df$Set <- "Site"
  df$part <- c("Model","Residual","Total")
  return(df)
}





################
#     16S      #
################


################################################################################################################################################################
#     ~ Sample.type / Population / Cluster
################################################################################################################################################################

tab_16s[1:5,1:5]
head(dmeta.16s)
relative_abundances <- tab_16s
metadata <- dmeta.16s

df.dist <- get_distance(relative_abundances)
df.source <- run_permanova_source(df.dist, metadata)
df.group <- run_permanova_group_strata(df.dist, metadata)
df.site <- run_permanova_site_groupsource(df.dist, metadata)
df.out <- rbind(df.source, df.group, df.site)
df.out
#write.table(df.out, "ordination_03_permanova_16s_groups.tab", sep="\t", row.names=TRUE, col.names=TRUE,quote=FALSE,append=FALSE)



#########################################
#     soil :   ~ Population / Cluster
#########################################

dmeta_sub <- dmeta.16s %>% filter(Sample.type == "bulk soil")
samples <- dmeta_sub$long.name
samples
tab_sub <- tab_16s %>% filter(row.names(tab_16s) %in% samples)
relative_abundances <- tab_sub
metadata <- dmeta_sub

df.dist <- get_distance(relative_abundances)
df.group <- run_permanova_group(df.dist, metadata)
df.site <- run_permanova_site_strata(df.dist, metadata)
df.out1 <- rbind(df.group, df.site)
df.out1
#write.table(df.out1, "ordination_03_permanova_16s_soil_only.tab", sep="\t", row.names=TRUE, col.names=TRUE,quote=FALSE,append=FALSE)


#########################################
#     roots :   ~ Population / Cluster
#########################################

dmeta_sub <- dmeta.16s %>% filter(Sample.type == "roots")
samples <- dmeta_sub$long.name
samples
tab_sub <- tab_16s %>% filter(row.names(tab_16s) %in% samples)
relative_abundances <- tab_sub
metadata <- dmeta_sub

df.dist <- get_distance(relative_abundances)
df.group <- run_permanova_group(df.dist, metadata)
df.site <- run_permanova_site_strata(df.dist, metadata)
df.out2 <- rbind(df.group, df.site)
df.out2
#write.table(df.out2, "ordination_03_permanova_16s_roots_only.tab", sep="\t", row.names=TRUE, col.names=TRUE,quote=FALSE,append=FALSE)


#########################################
#     rhizo :   ~ Population / Cluster
#########################################

dmeta.16s$Sample.type
dmeta_sub <- dmeta.16s %>% filter(Sample.type == "rhizosphere")
samples <- dmeta_sub$long.name
samples
tab_sub <- tab_16s %>% filter(row.names(tab_16s) %in% samples)
relative_abundances <- tab_sub
metadata <- dmeta_sub

df.dist <- get_distance(relative_abundances)
df.group <- run_permanova_group(df.dist, metadata)
df.site <- run_permanova_site_strata(df.dist, metadata)
df.out3 <- rbind(df.group, df.site)
df.out3
#write.table(df.out3, "ordination_03_permanova_16s_rhizo_only.tab", sep="\t", row.names=TRUE, col.names=TRUE,quote=FALSE,append=FALSE)


########################

df.out$type <- "all"
df.out1$type <- "soil"
df.out2$type <- "roots"
df.out3$type <- "rhizo"

res <- rbind(df.out, df.out1, df.out2, df.out3)
res <- res %>% filter(part == "Model")
res$p.adj <- p.adjust(res$`Pr(>F)`, "BH")
res

write.table(res, "ordination_03_permanova_16s_ALL.tab", sep="\t", row.names=TRUE, col.names=TRUE,quote=FALSE,append=FALSE)








################
#     ITS      #
################


################################################################################################################################################################
#     ~ Sample.type / Population / Cluster
################################################################################################################################################################

tab_its[1:5,1:5]
head(dmeta.its)
relative_abundances <- tab_its
metadata <- dmeta.its

df.dist <- get_distance(relative_abundances)
df.source <- run_permanova_source(df.dist, metadata)
df.group <- run_permanova_group_strata(df.dist, metadata)
df.site <- run_permanova_site_groupsource(df.dist, metadata)
df.out <- rbind(df.source, df.group, df.site)
df.out
#write.table(df.out, "ordination_03_permanova_ITS_groups.tab", sep="\t", row.names=TRUE, col.names=TRUE,quote=FALSE,append=FALSE)



#########################################
#     soil :   ~ Population / Cluster
#########################################

dmeta_sub <- dmeta.its %>% filter(Sample.type == "bulk soil")
samples <- dmeta_sub$long.name
samples
tab_sub <- tab_its %>% filter(row.names(tab_its) %in% samples)
relative_abundances <- tab_sub
metadata <- dmeta_sub

df.dist <- get_distance(relative_abundances)
df.group <- run_permanova_group(df.dist, metadata)
df.site <- run_permanova_site_strata(df.dist, metadata)
df.out1 <- rbind(df.group, df.site)
df.out1
#write.table(df.out, "ordination_03_permanova_ITS_soil_only.tab", sep="\t", row.names=TRUE, col.names=TRUE,quote=FALSE,append=FALSE)


#########################################
#     roots :   ~ Population / Cluster
#########################################

dmeta_sub <- dmeta.its %>% filter(Sample.type == "roots")
samples <- dmeta_sub$long.name
samples
tab_sub <- tab_its %>% filter(row.names(tab_its) %in% samples)
relative_abundances <- tab_sub
metadata <- dmeta_sub

df.dist <- get_distance(relative_abundances)
df.group <- run_permanova_group(df.dist, metadata)
df.site <- run_permanova_site_strata(df.dist, metadata)
df.out2 <- rbind(df.group, df.site)
df.out2
#write.table(df.out, "ordination_03_permanova_ITS_roots_only.tab", sep="\t", row.names=TRUE, col.names=TRUE,quote=FALSE,append=FALSE)


#########################################
#     rhizo :   ~ Population / Cluster
#########################################

dmeta.its$Sample.type
dmeta_sub <- dmeta.its %>% filter(Sample.type == "rhizosphere")
samples <- dmeta_sub$long.name
samples
tab_sub <- tab_its %>% filter(row.names(tab_its) %in% samples)
relative_abundances <- tab_sub
metadata <- dmeta_sub

df.dist <- get_distance(relative_abundances)
df.group <- run_permanova_group(df.dist, metadata)
df.site <- run_permanova_site_strata(df.dist, metadata)
df.out3 <- rbind(df.group, df.site)
df.out3
#write.table(df.out, "ordination_03_permanova_ITS_rhizo_only.tab", sep="\t", row.names=TRUE, col.names=TRUE,quote=FALSE,append=FALSE)



########################

df.out$type <- "all"
df.out1$type <- "soil"
df.out2$type <- "roots"
df.out3$type <- "rhizo"

res <- rbind(df.out, df.out1, df.out2, df.out3)
res <- res %>% filter(part == "Model")
res$p.adj <- p.adjust(res$`Pr(>F)`, "BH")
res


write.table(res, "ordination_03_permanova_ITS_ALL.tab", sep="\t", row.names=TRUE, col.names=TRUE,quote=FALSE,append=FALSE)



