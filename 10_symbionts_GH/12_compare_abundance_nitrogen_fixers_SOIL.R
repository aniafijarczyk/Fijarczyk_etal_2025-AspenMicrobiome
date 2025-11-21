rm(list = ls())

library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(phyloseq)
library(microbiome)
library(readxl)
library(stringr)
#library(rstatix) # Welsh anova

library(MASS)
library(nlme)
library(emmeans) # to get Tukey test on gls model

source("run_NESTED_ANOVA.R")


#########################
#--------DATA-----------#
#########################


df <- read.csv("./03_nitro_read.csv", header=TRUE, sep=',')
head(df)
df.meta <- read.csv("./03_nitro_read_meta.csv", header=TRUE, sep=',')
head(df.meta)



rownames(df) <- df$NF
otu.table <- df %>% dplyr::select(-NF)
head(otu.table)
OTU = otu_table(otu.table, taxa_are_rows = TRUE)

df.tax <- df %>% dplyr::select(NF)
colnames(df.tax) <- c("Family")
df.tax.mat <- df.tax %>% as.matrix()
TAX = tax_table(df.tax.mat)
TAX


head(df.meta)
df.meta.format <- df.meta %>% dplyr::select(long.name, sample, Cluster, Population, Sample.type) %>% distinct()
rownames(df.meta.format) <- df.meta.format$long.name
head(df.meta.format)
samdata <- sample_data(df.meta.format)

# Importing to phyloseq
pseq.stands <- phyloseq(OTU, TAX, sam_data = samdata)
pseq.stands




#####################
#---  TRANSFORM  ---#
#####################

df.meta.format$Sample.type
pseq.16s.soil <- subset_samples(pseq.stands, sample_data(pseq.stands)$Sample.type=="bulk soil")
pseq.16s.soil

# Convert to tab

# Transform
pseq.stands.rel <- microbiome::transform(pseq.16s.soil, 'clr')
pseq.stands.rel@otu_table


# Convert to dataframe
dtab.stands <- as.data.frame(pseq.stands.rel@otu_table)
class(dtab.stands)
head(dtab.stands)

# Gather
dtab.stands$family <- rownames(dtab.stands)
dg.stands <- dtab.stands %>% gather(key = "sample", value = "CLR abd",-family)
dg.stands

# Add meta
#df.meta.format <- df.meta %>% dplyr::rename(long.name = long.name.sequencing_R.,
#                                            Country = Country_x)
#df.meta.format <- df.meta.format %>% dplyr::select(long.name, short.name, site, Cluster, Country, Lattitude, Longitude) %>% distinct()
dg.meta.stands <- merge(dg.stands, df.meta.format, by.x = "sample", by.y = "long.name", sort=FALSE)
head(dg.meta.stands)

# Add taxonomy
df.taxonomy <- as.data.frame(pseq.stands.rel@tax_table)
df.taxonomy$family <- rownames(df.taxonomy)

dm.stands <- merge(dg.meta.stands, df.taxonomy, by = 'family', sort=FALSE)
head(dm.stands)



#####################
#-----  PLOTS ------#
#####################

fams <- c('NF')

dsub.stands <- dm.stands %>% filter(family %in% fams) 
dsub.stands %>% head()
dsub.stands$Label <- factor(dsub.stands$Family,levels=c("NF"))
head(dsub.stands)

write.table(dsub.stands, file="nitro_04_clr_soil.tsv", sep="\t", row.names = F, col.names = T, append=F, quote=F)


#####################
#-----   SITE ------#
#####################

# Plotting average abundance for each site

head(dsub.stands)
dsub.stands$Cluster
dsub.stands$site <- factor(dsub.stands$Cluster, levels = c("AB","FORE","FP","SP"))
dsub.stands$Population <- factor(dsub.stands$Population, levels = c("NENA","MX"))

p1 <- ggplot(dsub.stands) +
  geom_boxplot(data=dsub.stands, aes(x = Cluster, y = `CLR abd`, fill = Population), outlier.shape = NA) +
  geom_point(data=dsub.stands, aes(x = Cluster, y = `CLR abd`), color="grey40", pch=21, size=1)
p1





#####################
#-----  TESTS  -----#
#####################


###############
### TESTS of difference between countries

head(dsub.stands)

# Transform CLR so it's positive
min_clr <- min(dsub.stands$`CLR abd`)
min_clr
dsub.stands$CLR_abd <- dsub.stands$`CLR abd` + abs(min_clr) + 0.01
hist(dsub.stands$`CLR abd`)
hist(dsub.stands$CLR_abd)
plot(x = dsub.stands$`CLR abd`, y = dsub.stands$CLR_abd)

df.models <- run_nested_models(dsub.stands, "CLR_abd", "Population", "Cluster")
df.models

write.table(df.models, file = "nitro_04_clr_soil_tests.tsv", sep="\t", row.names = F, col.names = T, append=F, quote=F)




################
### PLOT



dsub.stands$Cluster <- factor(dsub.stands$Cluster, levels = c("AB","FORE","FP","SP"))
dsub.stands$Population <- factor(dsub.stands$Population, levels = c("NENA","MX"))
head(df.models)
head(dsub.stands)

dLM <- df.models %>% filter(MOD == "ANOVA")
#dLM[dLM["statistic"]=="p_GR1","V1"]

p1C <- ggplot(dsub.stands) +
  geom_boxplot(data=dsub.stands, aes(x = Cluster, y = `CLR abd`, fill = Population), outlier.shape = NA) +
  geom_point(data=dsub.stands, aes(x = Cluster, y = `CLR abd`), color="grey40",
             pch=21, size=6, position=position_jitter(width=0.1)) +
  
  annotate("text", x = 0.5, y = -0.5, 
                            label = paste0("Pop: p=",signif(dLM[dLM["statistic"]=="p_GR1","V1"],digits=2),
                                           "\nSite: p=",signif(dLM[dLM["statistic"]=="p_GR2","V1"],digits=2)),
            hjust=0, size=6) +


  scale_fill_manual(values = c("white","grey"))  + 
  labs(y = "CLR abundance") +
  #ggtitle("Population") + 
  scale_y_continuous(limits=c(-2.3,-0.3)) +
  theme(axis.title.y = element_text(size = 20),
        axis.title.x = element_blank(),
        plot.title = element_text(size=20, face="bold", hjust=0.5),
        panel.background = element_blank(), 
        #panel.border = element_rect(fill = NA, colour = "grey30"), 
        panel.border = element_blank(), 
        axis.ticks = element_blank(), legend.key = element_blank(),
        axis.text.x = element_text(angle=90, hjust=1,vjust=0.5, size=18),
        axis.text.y = element_text(size=14),
        legend.position="none",
  )
p1C


#png('nitro_04_clr_soil.png', w=800, h=1200, res=300)
#p1C
#dev.off()

