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
library(ANCOMBC)




#########################
#--------DATA-----------#
#########################



df_asv <- read.csv("ancombc_01_input_16S_ASV.tsv", sep="\t", header=T)
df_asv[1:10, 1:10]
df_asv$genus_taxon

df_tax <- read.csv("ancombc_01_input_16S_TAX.tsv", sep="\t", header=T)
df_tax %>% head()

df_meta <-  read.csv("../01_data_prep/dataPrep_01_overview_16S_META.tsv", sep="\t", header=T)
df_meta$Region0 <- ifelse(df_meta$site %in% c("AMOS","STFE"), "Boreal", df_meta$site)
df_meta$Region1 <- ifelse(df_meta$Region0 %in% c("STET","ESSI","FORE"), "Cold_temperate", df_meta$Region0)
df_meta$Region <- ifelse(df_meta$Region1 %in% c("FLOR1","FP2","Santiago"), "Warm_temperate", df_meta$Region1)
head(df_meta)

## Create phyloseq
row.names(df_asv) <- df_asv$genus_taxon
df_asv <- df_asv %>% dplyr::select(-genus_taxon)
row.names(df_asv)
head(df_asv)
OTU = otu_table(df_asv, taxa_are_rows = TRUE)
OTU

row.names(df_tax) <- df_tax$genus_taxon
df_tax <- df_tax %>% dplyr::select(-genus_taxon) %>% as.matrix()
head(df_tax)
TAX = tax_table(df_tax)
TAX

ds <- data.frame("long.name"=colnames(OTU))
dm <- merge(ds, df_meta, by = "long.name", sort=F)
rownames(dm) <- dm$long.name
head(dm)
samdata <- sample_data(dm)


pseq <- phyloseq(OTU, TAX, sam_data = samdata)
pseq

pseq.soil <- subset_samples(pseq, sample_data(pseq)$source=="Soil")
pseq.soil


row.names(OTU)



### Filtering

# Remove scarce species.
# A minimum treshold would be 0.05% for the less abundant species.
# Remove all species that occur in less than 5% of all samples

# Keep taxa present with min coverage of x or more in more than y fraction of the samples
filter_taxa_min_frac <- function(pseq.all, min.samp, frac.samp) {
  avs = genefilter_sample(pseq.all, filterfun_sample(function(x) x > min.samp), A=frac.samp*nsamples(pseq.all))
  pseq.new = prune_taxa(avs, pseq.all)
  return(pseq.new)
}

filter_taxa_min_n <- function(pseq.all, min.samp, n.samp) {
  avs = genefilter_sample(pseq.all, filterfun_sample(function(x) x > min.samp), A=n.samp)
  pseq.new = prune_taxa(avs, pseq.all)
  return(pseq.new)
}



### Dataset filtering

### Dataset filtering
pseq.soil
pseq.soil.2 = filter_taxa_min_frac(pseq.soil, 1, 0.35)
#pseq.soil.2 = filter_taxa_min_frac(pseq.soil, 1, 0.01)
pseq.soil.2

# Convert to tse
#tse.soil = mia::makeTreeSummarizedExperimentFromPhyloseq(pseq.soil.2)
tse.soil = mia::convertFromPhyloseq(pseq.soil.2)
tse.soil

rownames(tse.soil)










#########################
#-----ANCOMB-BC2--------#
#########################
sample_data(pseq.soil)



# https://bioconductor.org/packages/release/bioc/manuals/ANCOMBC/man/ANCOMBC.pdf

output.genus <- ancombc2(data = tse.soil, assay_name = "counts", tax_level = "Genus",
                         fix_formula = "Region", rand_formula = NULL,
                         p_adj_method = "holm", pseudo_sens = TRUE,
                         prv_cut = 0.10, lib_cut = 0, s0_perc = 0.05,
                         group = NULL, struc_zero = FALSE, neg_lb = FALSE,
                         alpha = 0.05, n_cl = 2, verbose = TRUE,
                         global = FALSE, pairwise = FALSE, 
                         dunnet = FALSE, trend = FALSE,
                         iter_control = list(tol = 1e-5, max_iter = 20, 
                                             verbose = FALSE),
                         em_control = list(tol = 1e-5, max_iter = 100),
                         lme_control = NULL,
                         mdfdr_control = list(fwer_ctrl_method = "holm", B = 100), 
                         trend_control = NULL)


#output.genus
res = output.genus$res
head(res)
res$taxon
rownames(tse.soil)
res$taxon_org <- rownames(tse.soil)
res %>% filter(diff_RegionCold_temperate == TRUE)
res %>% filter(diff_RegionWarm_temperate == TRUE)
#output.genus$feature_table

# Save results & feature table
write.table(res, "ancombc_03_run_16s_results.tsv", sep="\t",row.names = FALSE, col.names = TRUE, quote=FALSE, append=FALSE)
write.table(output.genus$feature_table, "ancombc_03_run_16s_feature_table.tsv", sep="\t",row.names = TRUE, col.names = TRUE, quote=FALSE, append=FALSE)





#################################################################################
### With site

output.genus.site <- ancombc2(data = tse.soil, assay_name = "counts", tax_level = "Genus",
                         fix_formula = "Region", rand_formula = "(1|site)",
                         p_adj_method = "holm", pseudo_sens = TRUE,
                         prv_cut = 0.10, lib_cut = 0, s0_perc = 0.05,
                         group = "Region", struc_zero = FALSE, neg_lb = FALSE,
                         alpha = 0.05,
                         n_cl = 2, verbose = TRUE,
                         global = TRUE, pairwise = TRUE, 
                         dunnet = FALSE, trend = FALSE,
                         iter_control = list(tol = 1e-5, max_iter = 20, 
                                             verbose = FALSE),
                         em_control = list(tol = 1e-5, max_iter = 100),
                         mdfdr_control = list(fwer_ctrl_method = "holm", B = 100), 
                         trend_control = NULL)


#output.genus.site
res = output.genus.site$res
res_pair = output.genus.site$res_pair
head(res)
head(res_pair)
dim(res_pair)
res_pair$taxon
# The assignment of unid taxa is WRONG - they do not correspond to ones in my table
#res_pair$taxon_org <- rownames(tse.soil)
res_pair %>% filter(diff_RegionWarm_temperate_RegionCold_temperate == TRUE)
res_pair %>% filter(diff_RegionCold_temperate == TRUE)
res_pair %>% filter(diff_RegionWarm_temperate == TRUE)

# Checking names
rownames(tse.soil)
rownames(tse.soil)[order(rownames(tse.soil))]
res_pair$taxon_org <- rownames(tse.soil)[order(rownames(tse.soil))]
name_check <- res_pair %>% dplyr::select(taxon, taxon_org)
for (i in 1:length(name_check$taxon)) {
  if (name_check$taxon[i] != name_check$taxon_org[i]) {
    print(name_check[i,])
  }
}

head(res_pair)
# Save results & feature table
write.table(res_pair, "ancombc_03_run_16s_results_random_effect.tsv", sep="\t",row.names = FALSE, col.names = TRUE, quote=FALSE, append=FALSE)
write.table(output.genus.site$feature_table, "ancombc_03_run_16s_feature_table_random_effect.tsv", sep="\t",row.names = TRUE, col.names = TRUE, quote=FALSE, append=FALSE)


