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


df_asv <- read.csv("ancombc_02_input_ITS_ASV.tsv", sep="\t", header=T)
df_asv[1:10, 1:10]
df_asv$genus_taxon

df_tax <- read.csv("ancombc_02_input_ITS_TAX.tsv", sep="\t", header=T)
df_tax %>% head()

df_meta <-  read.csv("../../2022_MicrobiomeAspenGreenhouse/01_data_prep/dataPrep_01_overview_ITS_META.tsv", sep="\t", header=T)
df_meta <- df_meta %>% filter(Population %in% c("NENA","MX"))
samples <- df_meta$long.name
samples
head(df_meta)


## Create phyloseq
row.names(df_asv) <- df_asv$genus_taxon
df_asv <- df_asv %>% dplyr::select(-genus_taxon) %>% dplyr::select(all_of(samples))
row.names(df_asv)
colnames(df_asv)
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
samdata

pseq <- phyloseq(OTU, TAX, sam_data = samdata)
pseq

pseq.soil <- subset_samples(pseq, sample_data(pseq)$Sample.type=="bulk soil")
pseq.soil


row.names(OTU)


#########################
#      Filtering        #
#########################

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
pseq
pseq.soil.2 = filter_taxa_min_frac(pseq.soil, 1, 0.35)
#pseq.soil.2 = filter_taxa_min_frac(pseq.soil, 1, 0.01)
pseq.soil.2

# Convert to tse
tse.soil = mia::convertFromPhyloseq(pseq.soil.2)
tse.soil

rownames(tse.soil)






#########################
#-----ANCOMB-BC2--------#
#########################
sample_data(pseq.soil)



# https://bioconductor.org/packages/release/bioc/manuals/ANCOMBC/man/ANCOMBC.pdf


### With random effect


output.genus.site <- ancombc2(data = tse.soil, assay_name = "counts", tax_level = "Genus",
                              fix_formula = "Population", rand_formula = "(1|Cluster)",
                              p_adj_method = "holm", pseudo_sens = TRUE,
                              prv_cut = 0.10, lib_cut = 0, s0_perc = 0.05,
                              group = "Population", struc_zero = FALSE, neg_lb = FALSE,
                              alpha = 0.05,
                              n_cl = 2, verbose = TRUE,
                              global = TRUE, pairwise = FALSE, 
                              dunnet = FALSE, trend = FALSE,
                              iter_control = list(tol = 1e-5, max_iter = 20, 
                                                  verbose = FALSE),
                              em_control = list(tol = 1e-5, max_iter = 100),
                              mdfdr_control = list(fwer_ctrl_method = "holm", B = 100), 
                              trend_control = NULL)


#output.genus.site
res = output.genus.site$res
head(res)
res$taxon
# The assignment of unid taxa is WRONG - they do not correspond to ones in my table
#res_pair$taxon_org <- rownames(tse.soil)
res %>% filter(diff_PopulationNENA == TRUE)

# Checking names
rownames(tse.soil)
rownames(tse.soil)[order(rownames(tse.soil))]
res$taxon_org <- rownames(tse.soil)[order(rownames(tse.soil))]
name_check <- res %>% dplyr::select(taxon, taxon_org)
for (i in 1:length(name_check$taxon)) {
  if (name_check$taxon[i] != name_check$taxon_org[i]) {
    print(name_check[i,])
  }
}

head(res)
# Save results & feature table
write.table(res, "ancombc_04_run_ITS_results_random_effect.tsv", sep="\t",row.names = FALSE, col.names = TRUE, quote=FALSE, append=FALSE)
write.table(output.genus.site$feature_table, "ancombc_04_run_ITS_feature_table_random_effect.tsv", sep="\t",row.names = TRUE, col.names = TRUE, quote=FALSE, append=FALSE)



