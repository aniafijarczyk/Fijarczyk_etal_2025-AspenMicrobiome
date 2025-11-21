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
library(ggnewscale)




#########################
#--------DATA-----------#
#########################


# Abundances
load("../01_data_prep/AVS_filtered_phyloseq.RData")
pseq.16s


# Subsetting soil samples
pseq.16s.soil <- subset_samples(pseq.16s, sample_data(pseq.16s)$source=="Soil")
pseq.16s.soil

dmeta.16s <- data.frame(sample_data(pseq.16s.soil))
dmeta.16s %>% head()

pseq.its.soil <- subset_samples(pseq.its, sample_data(pseq.its)$source=="Soil")
pseq.its.soil

dmeta.its <- data.frame(sample_data(pseq.its.soil))
dmeta.its %>% head()


########
### SOIL

data_soil <- read.xlsx(xlsxFile = "../../2022_MicrobiomeQuebecMexico/DATA/results/soil.xlsx")
data_soil %>% head()
data_soil[8:23] <- data_soil[8:23] %>% mutate_if(is.character, as.numeric)
data_soil <- data_soil %>% dplyr::select(-Country, -long.name)
data_soil$C_N <- data_soil$C_total/data_soil$N_total
dim(data_soil)
data_soil %>% head()

dmeta.16s.soil <- merge(dmeta.16s, data_soil, by = "short.name", sort=FALSE)
dmeta.16s.soil %>% head()
dim(dmeta.16s.soil)

dmeta.its.soil <- merge(dmeta.its, data_soil, by = "short.name",sort=FALSE)
dmeta.its.soil %>% head()
dim(dmeta.its.soil)



site_order <- c("AMOS","STFE","STET","ESSI","FORE","Santiago","FLOR1","FP2")

### Colors
source_cols <- c("#0072B2","#D55E00","#F0E442")
country_cols <- c("#E84A5F","#DDCC77")
site_cols <- brewer.pal(9,"BrBG")
phylum_cols <- c("#332288","#117733","#44AA99","#88CCEE","#DDCC77","#CC6677","#AA4499","#882255")
phylum_cols_wong <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

tab_blue <- "#1f77b4"
tab_orange <- "#ff7f0e"



#################
#---FILTERING---#
#################

# Keep taxa present with min coverage of x or more in more than y fraction of the samples
filter_taxa_min_frac <- function(pseq.all, min.samp, frac.samp) {
  avs = genefilter_sample(pseq.all, filterfun_sample(function(x) x > min.samp), A=frac.samp*nsamples(pseq.all))
  pseq.new = prune_taxa(avs, pseq.all)
  return(pseq.new)
}



######################
#---     16s      ---#
######################


##########################
#---  TRANSFORMATION  ---#
##########################

pseq.16s.2 = filter_taxa_min_frac(pseq.16s.soil, 1, 0.05)
pseq.16s.2.rel = microbiome::transform(pseq.16s.2, 'compositional') # transform to relative abundance



####################################
### Soil

# Table with soil properties
soil_metrics <- c("C_total","N_total","pH_H2O","P",
                  "K","Ca","Mg","Mn","Al","Fe",
                  "Na","CEC","Sable","Limon","Argile","C_N")

dim(dmeta.16s.soil)
dmeta.16s.soil$long.name
dmeta.16s.soil %>% head()
tab_soil <- dmeta.16s.soil %>% dplyr::select(all_of(soil_metrics))
tab_soil %>% head()
rownames(tab_soil) <- dmeta.16s.soil$long.name

# Renaming
tab_soil <- tab_soil %>% dplyr::rename(
  Sand = Sable,
  Silt = Limon,
  Clay = Argile,
  `C:N ratio` = C_N,
  `pH` = pH_H2O,
  `total C` = C_total,
  `total N` = N_total,
)
head(tab_soil)




##################### PCoA
pcoa <- ordinate(pseq.16s.2.rel, "PCoA", "bray")
pcoa$trace
pcoa$values %>% head()
pcoa$vectors[,1:2] %>% head()
pcoa_scores <- as.data.frame(pcoa$vectors[,1:2])
pcoa_scores


#####################
###     Envfit  

en = envfit(pcoa_scores, tab_soil, permutations = 1000, na.rm = TRUE)
en

df.envfit <- as.data.frame(en$vectors$arrows)
r2 <- en$vectors$r %>% as.vector()
r2
df.envfit$r2 <- r2
pval <- en$vectors$pvals %>% as.vector()
df.envfit$pval <- pval
df.envfit$sign <- ifelse(pval < 0.01, "yes", "no")
head(df.envfit)

write.table(df.envfit, "04_envfit_16s_table.tab", sep="\t", row.names=TRUE, col.names=TRUE,quote=FALSE,append=FALSE)


####### Getting table
data.scores <- as.data.frame(pcoa$vectors[,1:3])
length(rownames(dmeta.16s))
sum(rownames(data.scores) == rownames(dmeta.16s))
data.scores$Country = dmeta.16s$Country
data.scores$Site = dmeta.16s$site
head(data.scores)


en_coord_cont = as.data.frame(scores(en, "vectors")) * ordiArrowMul(en)
en_coord_cat = as.data.frame(scores(en, "factors")) * ordiArrowMul(en)
en_coord_cont$sign <- df.envfit$sign
en_coord_cont

# Percent variance explained by each axis
variance_explained <- 100 * pcoa$values$Relative_eig
variance_explained[1:2]
max(data.scores$Axis.1)
head(data.scores)


### colors
sites_cols <- c("#134b71","#1f77b4","#258dd5","#86c1ea","#d9ecf9","#ffbe84","#ff7f0e","#c15a00")

### permanova results
permanova <- data.frame("x" = c(max(data.scores$Axis.1)-0.3),
                        "y" = c(max(data.scores$Axis.2)+0.6),
                        "label" = paste0("Group: R2=0.13 ***\n",
                                         "Region: R2=0.22 ***\n",
                                         "Site: R2=0.46 ***"))


permanova

site_order <- c("AMOS","STFE","STET","ESSI","FORE","Santiago","FLOR1","FP2")
data.scores$Site <- factor(data.scores$Site, levels = site_order)


head(en_coord_cont)
en_coord_cont_sign <- en_coord_cont %>% filter(sign == "yes")

scale=0.7
gs <- ggplot(data = data.scores) + 
  geom_point(pch=21, data = data.scores, aes(x = Axis.1, y = Axis.2, fill = Site), size = 4, alpha = 1) + 
  scale_fill_manual(values = sites_cols)  + 
  geom_segment(aes(x = 0, y = 0, xend = Axis.1*scale, yend = Axis.2*scale), 
               data = en_coord_cont_sign, linewidth =0.3, alpha = 0.5, colour = "grey30",
               arrow=arrow(length = unit(0.1, "inches"))) +
 
  new_scale_color() +
  geom_text_repel(data = en_coord_cont_sign, aes(x = Axis.1*scale, y = Axis.2*scale),
            label = row.names(en_coord_cont_sign), color="black", size=6, force=0.001) + 

  labs(x = paste0("PCo 1"," (",round(variance_explained[1], 1),"%)"),
       y = paste0("PCo 2"," (",round(variance_explained[2], 1),"%)")) +
  
  
  geom_text(data=permanova, aes(x = x, y = y, label = label), hjust=0, size=6) +
  
  #scale_x_continuous(limits=c(-75,65)) +
  scale_y_continuous(limits=c(-0.7,1.2)) +
  theme(axis.title = element_text(size = 18), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.text = element_text(size = 16, colour = "grey30"),
        legend.position = c(0.14,0.78),
        )

gs


png("04_envfit_16s_site.png", w=1600,h=1600, res=300)
gs
dev.off()




#################################################################### envfit within clusters

head(sample_data(pseq.16s.2.rel))
pseq.country <- subset_samples(pseq.16s.2.rel, sample_data(pseq.16s.2.rel)$Country=="Canada")
pcoa <- ordinate(pseq.country, "PCoA", "bray")
pcoa_scores <- as.data.frame(pcoa$vectors[,1:2])
tab_soil.can <- tab_soil %>% filter(rownames(tab_soil) %in% sample_data(pseq.country)$long.name)
rownames(pcoa_scores) == rownames(tab_soil.can)

en = envfit(pcoa_scores, tab_soil.can, permutations = 1000, na.rm = TRUE)
df.envfit <- as.data.frame(en$vectors$arrows)
r2 <- en$vectors$r %>% as.vector()
df.envfit$r2 <- r2
pval <- en$vectors$pvals %>% as.vector()
df.envfit$pval <- pval

write.table(df.envfit, "04_envfit_16s_table_Canada.tab", sep="\t", row.names=TRUE, col.names=TRUE,quote=FALSE,append=FALSE)



################# 

pseq.country <- subset_samples(pseq.16s.2.rel, sample_data(pseq.16s.2.rel)$Country=="Mexico")
pcoa <- ordinate(pseq.country, "PCoA", "bray")
pcoa_scores <- as.data.frame(pcoa$vectors[,1:2])
tab_soil.can <- tab_soil %>% filter(rownames(tab_soil) %in% sample_data(pseq.country)$long.name)
rownames(pcoa_scores) == rownames(tab_soil.can)

en = envfit(pcoa_scores, tab_soil.can, permutations = 1000, na.rm = TRUE)
df.envfit <- as.data.frame(en$vectors$arrows)
r2 <- en$vectors$r %>% as.vector()
df.envfit$r2 <- r2
pval <- en$vectors$pvals %>% as.vector()
df.envfit$pval <- pval

write.table(df.envfit, "04_envfit_16s_table_Mexico.tab", sep="\t", row.names=TRUE, col.names=TRUE,quote=FALSE,append=FALSE)









######################
#---     ITS      ---#
######################

##########################
#---  TRANSFORMATION  ---#
##########################

pseq.its.2 = filter_taxa_min_frac(pseq.its.soil, 1, 0.05)
pseq.its.2.rel = microbiome::transform(pseq.its.2, 'compositional') # transform to relative abundance

####################################
### Soil

soil_metrics <- c("C_total","N_total","pH_H2O","P",
                  "K","Ca","Mg","Mn","Al","Fe",
                  "Na","CEC","Sable","Limon","Argile","C_N")

dim(dmeta.its.soil)
dmeta.its.soil$long.name
dmeta.its.soil %>% head()
tab_soil <- dmeta.its.soil %>% dplyr::select(all_of(soil_metrics))
tab_soil %>% head()
rownames(tab_soil) <- dmeta.its.soil$long.name

# Renaming
tab_soil <- tab_soil %>% dplyr::rename(
  Sand = Sable,
  Silt = Limon,
  Clay = Argile,
  `C:N ratio` = C_N,
  `pH` = pH_H2O,
  `total C` = C_total,
  `total N` = N_total,
)
head(tab_soil)




##################### PCoA
pcoa <- ordinate(pseq.its.2.rel, "PCoA", "bray")
pcoa$trace
pcoa$values %>% head()
pcoa$vectors[,1:2] %>% head()
pcoa_scores <- as.data.frame(pcoa$vectors[,1:2])
pcoa_scores


#####################
###     Envfit  

en = envfit(pcoa_scores, tab_soil, permutations = 1000, na.rm = TRUE)
en

df.envfit <- as.data.frame(en$vectors$arrows)
r2 <- en$vectors$r %>% as.vector()
r2
df.envfit$r2 <- r2
pval <- en$vectors$pvals %>% as.vector()
df.envfit$pval <- pval
df.envfit$sign <- ifelse(pval < 0.01, "yes", "no")
head(df.envfit)
df.envfit

write.table(df.envfit, "04_envfit_ITS_table.tab", sep="\t", row.names=TRUE, col.names=TRUE,quote=FALSE,append=FALSE)



####### Getting table
data.scores <- as.data.frame(pcoa$vectors[,1:3])
length(rownames(dmeta.its))
sum(rownames(data.scores) == rownames(dmeta.its))
data.scores$Country = dmeta.its$Country
data.scores$Site = dmeta.its$site
head(data.scores)


en_coord_cont = as.data.frame(scores(en, "vectors")) * ordiArrowMul(en)
en_coord_cat = as.data.frame(scores(en, "factors")) * ordiArrowMul(en)
en_coord_cont$sign <- df.envfit$sign
en_coord_cont

# Percent variance explained by each axis
variance_explained <- 100 * pcoa$values$Relative_eig
variance_explained[1:2]
max(data.scores$Axis.1)
head(data.scores)


### colors
sites_cols <- c("#134b71","#1f77b4","#258dd5","#86c1ea","#d9ecf9","#ffbe84","#ff7f0e","#c15a00")

### permanova results
permanova <- data.frame("x" = c(max(data.scores$Axis.1)-0.4),
                        "y" = c(max(data.scores$Axis.2)+0.4),
                        "label" = paste0("Group: R2=0.13 ***\n",
                                         "Region: R2=0.20 ***\n",
                                         "Site: R2=0.37 ***"))


permanova

site_order <- c("AMOS","STFE","STET","ESSI","FORE","Santiago","FLOR1","FP2")
data.scores$Site <- factor(data.scores$Site, levels = site_order)


head(en_coord_cont)
en_coord_cont_sign <- en_coord_cont %>% filter(sign == "yes")

scale=0.5
gs <- ggplot(data = data.scores) + 
  geom_point(pch=21, data = data.scores, aes(x = Axis.1, y = Axis.2, fill = Site), size = 4, alpha = 1) + 
  scale_fill_manual(values = sites_cols)  + 
  geom_segment(aes(x = 0, y = 0, xend = Axis.1*scale, yend = Axis.2*scale), 
               data = en_coord_cont_sign, linewidth =0.3, alpha = 0.5, colour = "grey30",
               arrow=arrow(length = unit(0.1, "inches"))) +
  
  new_scale_color() +
  geom_text_repel(data = en_coord_cont_sign, aes(x = Axis.1*scale, y = Axis.2*scale),
                  label = row.names(en_coord_cont_sign), color="black", size=6, force=0.001) + 
  
  labs(x = paste0("PCo 1"," (",round(variance_explained[1], 1),"%)"),
       y = paste0("PCo 2"," (",round(variance_explained[2], 1),"%)")) +
  
  
  geom_text(data=permanova, aes(x = x, y = y, label = label), hjust=0, size=6) +
  
  #scale_x_continuous(limits=c(-75,65)) +
  scale_y_continuous(limits=c(-0.7,1.1)) +
  theme(axis.title = element_text(size = 18), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.text = element_text(size = 16, colour = "grey30"),
        legend.position = "none",
  )

gs


png("04_envfit_ITS_site.png", w=1600,h=1600, res=300)
gs
dev.off()








#################################################################### envfit within clusters

head(sample_data(pseq.16s.2.rel))
pseq.country <- subset_samples(pseq.16s.2.rel, sample_data(pseq.16s.2.rel)$Country=="Canada")
pcoa <- ordinate(pseq.country, "PCoA", "bray")
pcoa_scores <- as.data.frame(pcoa$vectors[,1:2])
tab_soil.can <- tab_soil %>% filter(rownames(tab_soil) %in% sample_data(pseq.country)$long.name)
rownames(pcoa_scores) == rownames(tab_soil.can)

en = envfit(pcoa_scores, tab_soil.can, permutations = 1000, na.rm = TRUE)
df.envfit <- as.data.frame(en$vectors$arrows)
r2 <- en$vectors$r %>% as.vector()
df.envfit$r2 <- r2
pval <- en$vectors$pvals %>% as.vector()
df.envfit$pval <- pval

write.table(df.envfit, "04_envfit_16s_table_Canada.tab", sep="\t", row.names=TRUE, col.names=TRUE,quote=FALSE,append=FALSE)



################# 

pseq.country <- subset_samples(pseq.16s.2.rel, sample_data(pseq.16s.2.rel)$Country=="Mexico")
pcoa <- ordinate(pseq.country, "PCoA", "bray")
pcoa_scores <- as.data.frame(pcoa$vectors[,1:2])
tab_soil.can <- tab_soil %>% filter(rownames(tab_soil) %in% sample_data(pseq.country)$long.name)
rownames(pcoa_scores) == rownames(tab_soil.can)

en = envfit(pcoa_scores, tab_soil.can, permutations = 1000, na.rm = TRUE)
df.envfit <- as.data.frame(en$vectors$arrows)
r2 <- en$vectors$r %>% as.vector()
df.envfit$r2 <- r2
pval <- en$vectors$pvals %>% as.vector()
df.envfit$pval <- pval

write.table(df.envfit, "04_envfit_16s_table_Mexico.tab", sep="\t", row.names=TRUE, col.names=TRUE,quote=FALSE,append=FALSE)

