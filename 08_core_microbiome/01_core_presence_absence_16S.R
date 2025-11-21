rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(phyloseq)
library(microbiome)
library(cowplot)




#########################
#--------DATA-----------#
#########################

# Non rarefied taxa - as suggested in this review: https://www.pnas.org/doi/full/10.1073/pnas.2104429118
load("../01_data_prep/AVS_filtered_phyloseq.RData")
sample_data(pseq.16s) %>% head()
otu_table(pseq.16s)
pseq.16s

tab_blue <- "#1f77b4"
tab_orange <- "#ff7f0e"
site_order <- c("AMOS","STFE","STET","ESSI","FORE","Santiago","FLOR1","FP2")

# Subsetting soil samples
pseq.sub <- subset_samples(pseq.16s, sample_data(pseq.16s)$source=="Soil")
pseq.sub

# Metadata
dmeta <- data.frame(sample_data(pseq.sub))
dmeta %>% head()
dim(dmeta)
dmeta$long.name

# Transform to compositional abundance 
pseq.rel <- microbiome::transform(pseq.sub, "compositional")
pseq.rel



#============================================#
#             ASV frequencies                #
#============================================#


### Count at min >1/1000 detection 

# Min abundance 
prevs <- prevalence(pseq.rel, detection = 0, sort = TRUE, count=T) %>% as.data.frame()
colnames(prevs) <- "count"
prevs$ASV <- row.names(prevs)
head(prevs)


prevs.core <- prevs %>% filter(count >= 25)
dim(prevs.core)

## Subset otu table by ASV names
abd <- as.data.frame(otu_table(pseq.rel))
abd[1:10, 1:10]
dim(abd)
apply(abd, 2, sum)
asv.abd <- as.data.frame(apply(abd, 1, sum))
colnames(asv.abd) <- "asv_abd"
head(asv.abd)

#abd.sub <- abd %>% filter(row.names(abd) %in% row.names(prevs.core))
abd.sub <- abd
dim(abd.sub)
head(abd.sub)

abd.pres <- as.data.frame(lapply(abd.sub, function(x) ifelse(x <=0, 0, 1)))
row.names(abd.pres) <- row.names(abd.sub)
head(abd.pres)
dim(abd.pres)

dtax <- as.data.frame(tax_table(pseq.rel))
dm <- merge(abd.pres, dtax, by=0, sort=F)
dm2 <- merge(dm, prevs, by.x = "Row.names", by.y = "ASV", sort=F)
dm3 <- merge(dm2, asv.abd, by.x = "Row.names", by.y=0, sort=F)
dim(dm3)
dm3 %>% head()

write.table(dm3, file = "02_core_16S_sample_ASV.tsv", sep="\t", row.names=F, col.names=T, quote=F, append=F)


dtax %>% head()
dtax$Genus
dm$Genus

###############
#   By site
###############

head(abd.pres)

dt <- t(abd.pres) %>% as.data.frame()
head(dt)
dt$long.name <- row.names(dt)
dmeta.sub <- dmeta %>% dplyr::select(long.name, site)
dt2 <- merge(dt, dmeta.sub, by = "long.name", sort=F)
head(dt2)

gt2 <- dt2 %>% dplyr::select(-long.name) %>% group_by(site) %>% summarise_all(sum) %>% as.data.frame()
rownames(gt2) <- gt2$site
gt2 <- gt2 %>% dplyr::select(-site)
head(gt2)
gt3 <- as.data.frame(t(gt2))
head(gt3)
rownames(gt3) <- dm2$Row.names
head(gt3)

# Occurence in sites
prevs <- rowSums(gt3 > 0) %>% as.data.frame()
colnames(prevs) <- "count"
prevs$ASV <- row.names(prevs)
head(prevs)
#prevs %>% filter(count == 8)

gt4 <- as.data.frame(lapply(gt3, function(x) ifelse(x < 1, 0, 1)))
row.names(gt4) <- row.names(gt3)
head(gt4)
dim(gt4)

out <- merge(gt4, dtax, by=0, sort=F)
out2 <- merge(out, prevs, by.x = "Row.names", by.y = "ASV", sort=F)
out3 <- merge(out2, asv.abd, by.x = "Row.names", by.y = 0, sort=F)
out3 %>% head()

write.table(out3, file = "02_core_16S_site_ASV.tsv", sep="\t", row.names=F, col.names=T, quote=F, append=F)




#============================================#
#            Genus frequencies               #
#============================================#


# Get abundances 
abd <- read.csv("../15_ancombc_region/ancombc_01_input_16S_ASV.tsv", sep="\t", header=T)
#abd <- as.data.frame(otu_table(pseq.rel))
abd[1:10, 1:10]
dim(abd)


# Combine abundances
rownames(abd) <- abd$genus_taxon
head(abd)
abd <- abd %>% dplyr::select(-genus_taxon)
abd[1:10, 1:10]

genus.abd <- as.data.frame(apply(abd, 1, sum))
colnames(genus.abd) <- "genus_abd"
head(genus.abd)


# Prevelance of samples
#head(abd.rel)
prevs <- rowSums(abd > 0.00) %>% as.data.frame()
colnames(prevs) <- "count"
prevs$Genus <- row.names(prevs)
head(prevs)
prevs.core <- prevs %>% filter(count >= 25)
dim(prevs.core)


# Convert rel abd to presence/absence
abd.pres <- as.data.frame(lapply(abd, function(x) ifelse(x <=0, 0, 1)))
row.names(abd.pres) <- row.names(abd)
head(abd.pres)
dim(abd.pres)

dm2 <- merge(abd.pres, prevs, by.x = 0, by.y = "Genus", sort=F)
dm3 <- merge(dm2, genus.abd, by.x = "Row.names", by.y = 0, sort=F)
dm3 %>% head()
dim(dm3)

write.table(dm3, file = "02_core_16S_sample_GENUS.tsv", sep="\t", row.names=F, col.names=T, quote=F, append=F)



###############
#   By site
###############

head(dm2)

dm4 <- dm2 %>% dplyr::select(-Row.names, -count)
head(dm4)
dt <- t(dm4) %>% as.data.frame()
head(dt)
dt$long.name <- row.names(dt)
dmeta.sub <- dmeta %>% dplyr::select(long.name, site)
dt2 <- merge(dt, dmeta.sub, by = "long.name", sort=F)
head(dt2)

gt2 <- dt2 %>% dplyr::select(-long.name) %>% group_by(site) %>% summarise_all(sum) %>% as.data.frame()
rownames(gt2) <- gt2$site
gt2 <- gt2 %>% dplyr::select(-site)
head(gt2)
gt3 <- as.data.frame(t(gt2))
head(gt3)
rownames(gt3) <- dm2$Row.names
head(gt3)

gt4 <- as.data.frame(lapply(gt3, function(x) ifelse(x < 1, 0, 1)))
row.names(gt4) <- row.names(gt3)
head(gt4)
dim(gt4)

# Occurence in sites
prevs <- rowSums(gt3 > 0) %>% as.data.frame()
colnames(prevs) <- "count"
head(prevs)
prevs %>% filter(count == 8)

out <- merge(gt4, prevs, by = 0, sort=F)
out2 <- merge(out, genus.abd, by.x = "Row.names", by.y=0, sort=F)
out2 %>% head()

write.table(out2, file = "02_core_16S_site_GENUS.tsv", sep="\t", row.names=F, col.names=T, quote=F, append=F)




