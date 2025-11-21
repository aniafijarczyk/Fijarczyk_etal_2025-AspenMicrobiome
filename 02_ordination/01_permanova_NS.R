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
load("../01_data_prep/AVS_filtered_phyloseq.RData")
pseq.16s


### 16S

# Subsetting soil samples
pseq.16s.soil <- subset_samples(pseq.16s, sample_data(pseq.16s)$source=="Soil")
pseq.16s.soil

dmeta.16s <- data.frame(sample_data(pseq.16s.soil))
dmeta.16s %>% head()

# Filtering
pseq.16s.2 = filter_taxa_min_frac(pseq.16s.soil, 1, 0.05)
pseq.16s.2.rel = microbiome::transform(pseq.16s.2, 'compositional') # transform to relative abundance

tab_16s <- as.data.frame(t(pseq.16s.2.rel@otu_table))
tab_16s[1:10, 1:10]




### ITS

pseq.its.soil <- subset_samples(pseq.its, sample_data(pseq.its)$source=="Soil")
pseq.its.soil

dmeta.its <- data.frame(sample_data(pseq.its.soil))
dmeta.its %>% head()

# Filtering
pseq.its.2 = filter_taxa_min_frac(pseq.its.soil, 1, 0.05)
pseq.its.2.rel = microbiome::transform(pseq.its.2, 'compositional') # transform to relative abundance

tab_its <- as.data.frame(t(pseq.its.2.rel@otu_table))
tab_its[1:10, 1:10]



########
### SOIL

data_soil <- read.xlsx(xlsxFile = "../../2022_MicrobiomeQuebecMexico/DATA/results/soil.xlsx")
data_soil %>% head()
data_soil[8:23] <- data_soil[8:23] %>% mutate_if(is.character, as.numeric)
data_soil %>% head()
data_soil$Region0 <- ifelse(data_soil$Site %in% c("AMOS","STFE"), "Boreal", data_soil$Site)
data_soil$Region1 <- ifelse(data_soil$Region0 %in% c("STET","ESSI","FORE"), "Cold_temperate", data_soil$Region0)
data_soil$Region <- ifelse(data_soil$Region1 %in% c("FLOR1","FP2","Santiago"), "Warm_temperate", data_soil$Region1)
data_soil %>% dplyr::select(Site, Region, Country) %>% distinct()
data_soil <- data_soil %>% dplyr::select(-Country, -long.name, -Region0, -Region1)
dmeta.16s.soil <- merge(dmeta.16s, data_soil, by = "short.name", sort=FALSE)
dmeta.16s.soil %>% head()
dim(dmeta.16s.soil)

dmeta.its.soil <- merge(dmeta.its, data_soil, by = "short.name",sort=FALSE)
dmeta.its.soil %>% head()
dim(dmeta.its.soil)





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


###########
### 16S ###
###########

# Checking of rownames of two datasets are the same
rownames(tab_16s)
dmeta.16s.soil %>% head()
dmeta.16s.soil$long.name

### Distances
bc_dist <- vegan::vegdist(tab_16s, method = "bray")
bc_dist


### COUNTRY ###
perm.country <- adonis2(bc_dist ~ Country, data=dmeta.16s.soil, permutations=999)
perm.country
df.country <- as.data.frame(perm.country)
df.country$Set <- "Country"
df.country$part <- c("Model","Residual","Total")

### REGION ###
perm.region <- adonis2(bc_dist ~ Region, data=dmeta.16s.soil, permutations=999, strata = dmeta.16s.soil$Country)
perm.region
df.region <- as.data.frame(perm.region)
df.region$Set <- "Region"
df.region$part <- c("Model","Residual","Total")
df.region

### SITE ###
perm.site <- adonis2(bc_dist ~ Site, data=dmeta.16s.soil, permutations=999, strata = dmeta.16s.soil$Region)
perm.site
df.site <- as.data.frame(perm.site)
df.site$Set <- "Site"
df.site$part <- c("Model","Residual","Total")
df.site

### combining outputs
df.output <- rbind(df.country, df.region, df.site)
df.output

write.table(df.output, "03_permanova_16s_3models.tab", sep="\t", row.names=TRUE, col.names=TRUE,quote=FALSE,append=FALSE)



bc_dist <- vegan::vegdist(tab_16s, method = "bray")
perm.3 <- adonis2(bc_dist ~ Country / Region/ Site, data=dmeta.16s.soil, permutations=999, strata=dmeta.16s.soil$Country, by="margin")
perm.3 <- adonis2(bc_dist ~ Country:Region:Site, data=dmeta.16s.soil, permutations=999, by="margin")
perm.3
perm.3 <- adonis2(bc_dist ~ Country:Region, data=dmeta.16s.soil, permutations=999, by="margin")
perm.3



###########
### ITS ###
###########

# Checking of rownames of two datasets are the same
rownames(tab_its)
dmeta.its.soil %>% head()
dmeta.its.soil$long.name

### Distances
bc_dist <- vegan::vegdist(tab_its, method = "bray")
bc_dist


### COUNTRY ###
perm.country <- adonis2(bc_dist ~ Country, data=dmeta.its.soil, permutations=999)
perm.country
df.country <- as.data.frame(perm.country)
df.country$Set <- "Country"
df.country$part <- c("Model","Residual","Total")

### REGION ###
perm.region <- adonis2(bc_dist ~ Region, data=dmeta.its.soil, permutations=999, strata = dmeta.its.soil$Country)
perm.region
df.region <- as.data.frame(perm.region)
df.region$Set <- "Region"
df.region$part <- c("Model","Residual","Total")
df.region

### SITE ###
perm.site <- adonis2(bc_dist ~ Site, data=dmeta.its.soil, permutations=999, strata = dmeta.its.soil$Region)
perm.site
df.site <- as.data.frame(perm.site)
df.site$Set <- "Site"
df.site$part <- c("Model","Residual","Total")
df.site

### combining outputs
df.output <- rbind(df.country, df.region, df.site)
df.output

write.table(df.output, "03_permanova_ITS_3models.tab", sep="\t", row.names=TRUE, col.names=TRUE,quote=FALSE,append=FALSE)

