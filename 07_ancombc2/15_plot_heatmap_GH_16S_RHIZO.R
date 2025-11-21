rm(list = ls())
#setwd("/media/anna/ADATA/NRCan/2022_MicrobiomeQuebecMexico/02_data_prep")
setwd("D:/NRCan/2024_Aspen/20_ancombc_gh_rhizo")

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
library(vegan)
library(ggnewscale)

#########################
#--------DATA-----------#
#########################

# log fold changes
res <- read.csv("ancombc_03_run_16S_results_random_effect.tsv", sep="\t", header=T)
res <- res %>% dplyr::select(-taxon)
dim(res)
head(res)
# Combine with short taxon name 
tax <- read.csv("../19_ancombc_gh_soil/ancombc_01_input_16S_TAX.tsv", sep="\t", header=T)
tax <- tax %>% dplyr::select(Genus, taxon)
head(tax)
mres <- merge(res, tax, by.x = "taxon_org", by.y = "Genus", sort=F)
head(mres)


# metadata
dmeta <- read.csv("../../2022_MicrobiomeAspenGreenhouse/01_data_prep/dataPrep_01_overview_16S_META.tsv", sep="\t", header=T)
#dmeta$Region0 <- ifelse(dmeta$site %in% c("AMOS","STFE"), "Boreal", dmeta$site)
#dmeta$Region1 <- ifelse(dmeta$Region0 %in% c("STET","ESSI","FORE"), "Cold_temperate", dmeta$Region0)
#dmeta$Region <- ifelse(dmeta$Region1 %in% c("FLOR1","FP2","Santiago"), "Warm_temperate", dmeta$Region1)
head(dmeta)
samples <- dmeta %>% filter(Sample.type == "rhizosphere") %>% filter(Population %in% c("NENA","MX")) %>% pull(long.name)
length(samples)
samples

# abundances
abd <- read.csv("../19_ancombc_gh_soil/ancombc_01_input_16S_ASV.tsv", sep="\t", header=T)
#abd <- abd %>% filter(taxon %in% res$taxon_org)
row.names(abd) <- abd$genus_taxon
head(abd)
abd <- abd %>% dplyr::select(-genus_taxon) %>% dplyr::select(any_of(samples))
dim(abd)
abd[1:10, 1:10]
ASV <- data.frame(t(abd))
ASV[1:10, 1:10]


# Remove unidentified_1 (non fungi)



#########################
#------- FILTER --------#
#########################

# First, verify that all samples contain the same amount of OTU: -> NO
apply(ASV, 1, sum)

# Caluclate relative abundance (in %)
sample_summed <- apply(ASV, 1, sum)
ASV_rel <- sweep(ASV, 1, sample_summed, "/")
ASV_rel[1:10, 1:10]
ASV_rel <- ASV_rel*100
ASV_rel[1:10, 1:10]
apply(ASV_rel,1,sum)


# Remove scarce species.
# A minimum treshold would be 0.05% for the less abundant species.
ASV_mean <- apply(ASV_rel, 2, mean)
summary(ASV_mean)
ASV_scarce <- ASV_rel[ , ! ASV_mean < 0.05]
dim(ASV_scarce)

# Remove rare species.
ASV_pa <- decostand(ASV_scarce, "pa") # "pa" stands for presence-absence
ASV_sum <- apply(ASV_pa, 2, sum) # sums the recurrence of species
# Remove all species that occur in less than 5% of all samples
ASV_rare <- ASV_scarce[ , ! ASV_sum < nrow(ASV_scarce)*0.05] # treshold of 5%
dim(ASV_rare)
head(ASV_rare)
apply(ASV_rare,1,sum)
dim(ASV_rare)
ASV_rare[1:10, 1:10]




#########################
#------ TOP TAXA -------#
#########################

dab <- as.data.frame(apply(ASV_rare,2,sum))
colnames(dab) <- "sum_abd"
dab$taxon <- rownames(dab)
dab
dab <- dab %>% arrange(-sum_abd)
dab$rank <- c(1:length(dab$taxon))
head(dab, n=30)




#########################
#------  FORMAT  -------#
#########################

### 
ASV_rare %>% head()
df_asv <- ASV_rare
df_asv$SampleID <- rownames(df_asv)
dg_abd <- df_asv %>% gather(-SampleID, key = "taxon", value = "rel_abd")
head(dg_abd)

# Adding rank
dg_abd <- merge(dg_abd, dab, by = "taxon", sort=F)
head(dg_abd)

# Adding metadata
dg_abd <- merge(dg_abd, dmeta, by.x = "SampleID", by.y = "long.name", sort=F)
head(dg_abd)
dg_abd <- dg_abd %>% dplyr::rename(taxon_org = taxon)

mg_abd <- merge(dg_abd, tax, by.x = "taxon_org", by.y = "Genus", sort=F)
head(mg_abd)

write.table(mg_abd, "ancombc_06_plot_16S_abundances_with_meta.tsv", sep="\t", row.names=F, col.names=T, quote=F, append=F)


### Region where each taxon is most abundant
mg_abd %>% head()
dh <- mg_abd %>% group_by(taxon_org, Population) %>% dplyr::summarise(mean_reg_abd = mean(rel_abd)) %>% arrange(taxon_org, -mean_reg_abd)
#%>% spread(Region, mean_reg_abd)
dh2 <- dh %>% group_by(taxon_org) %>% dplyr::mutate(top = first(mean_reg_abd))
dh3 <- dh2 %>% dplyr::filter(mean_reg_abd == top)
dh3 <- dh3 %>% dplyr::select(-top)
head(dh3)
dim(dh3)

#########################
#--- LOG FOLD CHANGE ---#
#########################

# lfc_RegionCold_temperate : lfc Cold temperate - boreal (positive are more abundant in cold temperate) 
# lfc_RegionWarm_temperate : lfc Warm temperate - boreal (positive are more abundant in warm temperate) 
# lfc_RegionWarm_temperate_RegionCold_temperate : lfc Warm temperate - Cold temperate (positive are more abundant in warm temperate) 

# lfc_PopulationNENA : lfc NENA - MX (positive are more abundant in Canada)

### Adding ranks
head(mres)
dab <- dab %>% dplyr::rename(taxon_org = taxon)
head(dab)

wres <- merge(mres, dab, by = "taxon_org", sort=F)
head(wres)
dim(wres)

tres <- merge(wres, dh3, by = "taxon_org", sort=F)
head(tres)

## Merging with ranks
#mmres <- merge(mres, dab, by.x = "taxon.y", by.y = "taxon", sort=F)
#head(mmres)




##########################
#   SELECTING TOP TAXA   #
##########################

head(tres)
tres %>% filter(endsWith(taxon_org,"_Bacteria"))
tres %>% filter(endsWith(taxon_org,"_Archaea"))

top <- tres %>% filter(!endsWith(taxon_org,"_Bacteria")) %>% filter(!endsWith(taxon_org,"_Archaea")) %>% arrange(rank) %>% head(20)
top %>% head()


top <- top %>% arrange(Population, -mean_reg_abd)
head(top)
top_ordered <- top %>% dplyr::select(taxon_org, lfc_PopulationNENA,
                                     p_PopulationNENA,  
                                     q_PopulationNENA,
                                     Population, mean_reg_abd) %>% distinct() %>% arrange(lfc_PopulationNENA)
top_ordered %>% head()
#top$taxon_org <- factor(top$taxon_org, levels = top_ordered$taxon_org)
#top$Country <- ifelse(top$lfc_CountryMexico>0,"Mexico","Canada")
#top$sign <- paste0(top$Country,"-",top$diff_CountryMexico)

head(top)
#top$pval <- ifelse(top$p_CountryMexico<0.05,0,1)
dim(top)
head(top)




#################################
#----      Z-SCORES       ------#
#################################

head(mg_abd)
dg.site <- mg_abd %>% group_by(taxon_org, taxon, Cluster, Population) %>% dplyr::summarise(mean_abd = mean(rel_abd), rank = first(rank),sum=sum(rel_abd))
head(dg.site)


#zscores <- ddply(dg.site, .(taxon), summarize, z_score=scale(mean_abd))
zscores <- dg.site %>% group_by(taxon_org) %>% dplyr::reframe(z_score=(sum-mean(sum))/sd(sum))
head(zscores)

zscoresm <- dg.site %>% group_by(taxon_org) %>% dplyr::reframe(z_score=(mean_abd-mean(mean_abd))/sd(mean_abd))
head(zscoresm)

#head(zscores)
dg.site$zscore <- zscores$z_score
dg.site$zscore_m <- zscoresm$z_score
head(dg.site)

top_ordered$taxon_org
top_avs <- dg.site %>% filter(taxon_org %in% top_ordered$taxon_org)
dim(top_avs)
head(top_avs)
length(top_ordered$taxon_org)
length(unique(top_avs$taxon))
length(top_avs$taxon_org)

head(top_avs)
dim(top_avs)

spread_asv <- top_avs %>% dplyr::select(taxon, Cluster, zscore_m) %>% spread(Cluster, zscore_m)
head(spread_asv)

#################################
#----   HEATMAP BY SITE   ------#
#################################


head(top_avs)
dord <- data.frame("taxon_org" = top_ordered$taxon_org)
mord <- merge(dord, top_avs, by = "taxon_org", sort=F)
tax_order <- mord %>% pull(taxon) %>% unique()
tax_order

top_avs$taxon <- factor(top_avs$taxon, levels = tax_order)
#taxon_order <- top %>% dplyr::select(taxon, rank) %>% distinct() %>% arrange(rank) %>% pull(taxon)
#top_avs$taxon <- factor(top_avs$taxon, levels = top_ordered$taxon_org)
top_avs$Site <- factor(top_avs$Cluster, levels = c("AB","FORE","FP","SP"), labels = c("AB","FORE","FP2","SP"))
#top_avs$Country <- factor(top_avs$Country, levels = c("Canada","Mexico"))
top_avs$Pop <- factor(top_avs$Population, levels = c("NENA","MX"), labels = c("Canada","Mexico"))

head(top_avs)
top_avs$group <- "16S"
head(top_avs)

p6 <- ggplot(top_avs) + aes(x = Site, y = taxon) +
  geom_tile(aes(fill = zscore_m), color='white') +
  #coord_fixed() +
  #scale_fill_distiller(palette = "Spectral", name="Relative\nabundance") +
  scale_fill_distiller(palette = "Greys", name="Abundance\nz-score",direction=1) +
  ylab('') +
  xlab('Site') +
  #facet_grid(~ Country) +
  facet_grid(~ Pop, scales = "free_x", space = "free",
             as.table = TRUE,drop = TRUE, switch="x") +
  theme(panel.background = element_blank(),
        #strip.placement = "outside",
        strip.text = element_blank(),
        legend.position = "top",
        legend.key.size = unit(0.4, 'cm'),
        strip.background = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(angle=90,hjust=1),
  )
p6




#############################
#----   HEATMAP LCF   ------#
#############################

spread_asv <- top_avs %>% dplyr::select(taxon, Cluster, zscore_m) %>% spread(Cluster, zscore_m)
head(spread_asv)
dim(spread_asv)

spread_asv_lfc <- merge(spread_asv, top, by = c("taxon_org","taxon"), sort=F)
head(spread_asv_lfc)

# Change taxon name (remove unind_*_ tag)

vtax <- c() 
for (i in 1:length(tax_order)) {
  taxid <- as.character(tax_order[i])
  print(taxid)
  if (str_sub(taxid,1,5) == "unid_") {
    ntax <- paste0("uc_",unlist(str_split(taxid,"_"))[3])
  }
  else {
    ntax <- taxid
  }
  print(ntax)
  vtax[i] <- ntax
}


tab_blue <- "#1f77b4"
tab_orange <- "#ff7f0e"

spread_asv_lfc$taxon <- factor(spread_asv_lfc$taxon, levels = tax_order, labels = vtax)
head(spread_asv_lfc)



p1 <- ggplot(spread_asv_lfc) + 
  
  geom_tile(data=spread_asv_lfc, aes(x = 1, y = taxon, fill = AB)) +
  geom_tile(data=spread_asv_lfc, aes(x = 2, y = taxon, fill = FORE)) +
  
  geom_tile(data=spread_asv_lfc, aes(x = 4, y = taxon, fill = FP)) +
  geom_tile(data=spread_asv_lfc, aes(x = 5, y = taxon, fill = SP)) +
  
  scale_fill_distiller(palette = "Greys", name="Relative\nabundance\nz-score",direction=1,
                       guide = guide_legend(direction="horizontal", ncol=10,order = 1)) +
  
  new_scale_fill() +
  geom_tile(data=spread_asv_lfc, aes(x = 7, y = taxon, fill = lfc_PopulationNENA)) +
  scale_fill_gradient2(high="#1f77b4", mid="white", low="#ff7f0e", name="LFC\nCanada-\nMexico",
                       guide = guide_legend(direction="horizontal", order = 2, nrow=1),limits=c(-2,2)) +
  coord_fixed() +
  
  labs(y = "Genus") +
  
  theme(panel.background = element_blank(),
        axis.text.x = element_text(size=18,angle=90,hjust=1,vjust=0.5),
        axis.text.y = element_text(size=18),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=20),
        legend.title = element_text(size=18),
        legend.text = element_text(size=16),
        legend.spacing.y = unit(0.001,"cm"),
        legend.key.spacing = unit(0.0,"cm"),
        legend.text.position = "bottom",
        legend.title.position = "top") +
  
  scale_x_discrete(breaks = c(1,2,4,5,7), limits = c(1,2,4,5,7),
                   labels=c("AB","FORE","FP2","Santiago","CAN\nvs. MX"))
  
p1
# gives warning but works!

png('ancombc_06_plot_16S_heatmap.png', w=2800,h=1600, res=300)
p1
dev.off()


