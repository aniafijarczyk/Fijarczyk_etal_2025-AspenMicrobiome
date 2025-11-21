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
tax <- read.csv("ancombc_01_input_16S_TAX.tsv", sep="\t", header=T)
tax <- tax %>% dplyr::select(Genus, taxon)
head(tax)
mres <- merge(res, tax, by.x = "taxon_org", by.y = "Genus", sort=F)
head(mres)


# metadata
dmeta <- read.csv("../01_data_prep/dataPrep_01_overview_16S_META.tsv", sep="\t", header=T)
dmeta$Region0 <- ifelse(dmeta$site %in% c("AMOS","STFE"), "Boreal", dmeta$site)
dmeta$Region1 <- ifelse(dmeta$Region0 %in% c("STET","ESSI","FORE"), "Cold_temperate", dmeta$Region0)
dmeta$Region <- ifelse(dmeta$Region1 %in% c("FLOR1","FP2","Santiago"), "Warm_temperate", dmeta$Region1)
head(dmeta)
samples <- dmeta %>% filter(source == "Soil") %>% pull(long.name)
length(samples)


# abundances
abd <- read.csv("ancombc_01_input_16S_ASV.tsv", sep="\t", header=T)
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
dh <- mg_abd %>% group_by(taxon_org, Region) %>% dplyr::summarise(mean_reg_abd = mean(rel_abd)) %>% arrange(taxon_org, -mean_reg_abd)
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

top <- tres %>% filter(!endsWith(taxon_org,"_Bacteria")) %>% filter(!endsWith(taxon_org,"_Archaea")) %>% arrange(rank) %>% head(30)
top %>% head()


top <- top %>% arrange(Region, -mean_reg_abd)
head(top)
top_ordered <- top %>% dplyr::select(taxon_org, lfc_RegionCold_temperate,lfc_RegionWarm_temperate, lfc_RegionWarm_temperate_RegionCold_temperate,
                                     p_RegionCold_temperate, p_RegionWarm_temperate, p_RegionWarm_temperate_RegionCold_temperate,  
                                     q_RegionCold_temperate, q_RegionWarm_temperate, q_RegionWarm_temperate_RegionCold_temperate,
                                     Region, mean_reg_abd) %>% distinct() %>% arrange(Region, -mean_reg_abd)
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
dg.site <- mg_abd %>% group_by(taxon_org, taxon, site, Region) %>% dplyr::summarise(mean_abd = mean(rel_abd), rank = first(rank),sum=sum(rel_abd))
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

spread_asv <- top_avs %>% dplyr::select(taxon, site, zscore_m) %>% spread(site, zscore_m)
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
top_avs$site <- factor(top_avs$site, levels = c("AMOS","STFE","STET","ESSI","FORE","Santiago","FLOR1","FP2"))
#top_avs$Country <- factor(top_avs$Country, levels = c("Canada","Mexico"))
top_avs$Region <- factor(top_avs$Region, levels = c("Boreal","Cold_temperate","Warm_temperate"))

head(top_avs)
top_avs$group <- "16S"
head(top_avs)

p6 <- ggplot(top_avs) + aes(x = site, y = taxon) +
  geom_tile(aes(fill = zscore_m), color='white') +
  #coord_fixed() +
  #scale_fill_distiller(palette = "Spectral", name="Relative\nabundance") +
  scale_fill_distiller(palette = "Greys", name="Abundance\nz-score",direction=1) +
  ylab('') +
  xlab('Site') +
  #facet_grid(~ Country) +
  facet_grid(~ Region, scales = "free_x", space = "free",
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

spread_asv <- top_avs %>% dplyr::select(taxon, site, zscore_m) %>% spread(site, zscore_m)
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
head(spread_asv)

#write.table(spread_asv_lfc, file="ancombc_06_plot_16S_table.tsv", sep="\t", col.names = T,row.names = F, append=F, quote=F)



p1 <- ggplot(spread_asv_lfc) + 
  
  geom_tile(data=spread_asv_lfc, aes(x = 1, y = taxon, fill = AMOS)) +
  geom_tile(data=spread_asv_lfc, aes(x = 2, y = taxon, fill = STFE)) +
  
  geom_tile(data=spread_asv_lfc, aes(x = 4, y = taxon, fill = STET)) +
  geom_tile(data=spread_asv_lfc, aes(x = 5, y = taxon, fill = ESSI)) +
  geom_tile(data=spread_asv_lfc, aes(x = 6, y = taxon, fill = FORE)) +
  
  geom_tile(data=spread_asv_lfc, aes(x = 8, y = taxon, fill = Santiago)) +
  geom_tile(data=spread_asv_lfc, aes(x = 9, y = taxon, fill = FLOR1)) +
  geom_tile(data=spread_asv_lfc, aes(x = 10, y = taxon, fill = FP2)) +
  scale_fill_distiller(palette = "Greys", name="Relative\nabundance\nz-score",direction=1,
                       guide = guide_legend(direction="horizontal", ncol=10,order = 1)) +
  
  new_scale_fill() +
  geom_tile(data=spread_asv_lfc, aes(x = 12, y = taxon, fill = lfc_RegionCold_temperate)) +
  scale_fill_gradient2(high="orchid4", mid="white", low="#1f77b4", name="LFC\nCold temp-\nBoreal",
                       guide = guide_legend(direction="horizontal", order = 2, nrow=1), limits=c(-2,2)) +
  
  new_scale_fill() +
  geom_tile(data=spread_asv_lfc, aes(x = 13, y = taxon, fill = lfc_RegionWarm_temperate)) +
  scale_fill_gradient2(high="#ff7f0e", mid="white", low="#1f77b4", name="LFC\nWarm temp-\nBoreal",
                       guide = guide_legend(direction="horizontal", order = 3,nrow=1),limits=c(-1,3)) +
  
  new_scale_fill() +
  geom_tile(data=spread_asv_lfc, aes(x = 14, y = taxon, fill = lfc_RegionWarm_temperate_RegionCold_temperate)) +
  scale_fill_gradient2(high="#ff7f0e", mid="white", low="orchid4", name="LFC\nWarm temp-\nCold temp",
                       guide = guide_legend(direction="horizontal", order = 4,nrow=1),limits=c(-1,3.5)) +
  coord_fixed() +
  
  labs(y = "Genus") +
  
  theme(panel.background = element_blank(),
        axis.text.x = element_text(size=16,angle=90,hjust=1,vjust=0.5),
        axis.text.y = element_text(size=16),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=18),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        legend.spacing.y = unit(0.001,"cm"),
        legend.key.spacing = unit(0.0,"cm"),
        legend.text.position = "bottom",
        legend.title.position = "top") +
  
  scale_x_discrete(breaks = c(1,2,4,5,6,8,9,10,12,13,14), limits = c(1,2,4,5,6,8,9,10,12,13,14),
                   labels=c("AMOS","STFE","STET","ESSI","FORE","Santiago","FLOR1","FP2",
                            "CT vs B", "WT vs B", "WT vs CT"))
  
p1

png('ancombc_06_plot_16S_heatmap.png', w=2800,h=2400, res=300)
p1
dev.off()


