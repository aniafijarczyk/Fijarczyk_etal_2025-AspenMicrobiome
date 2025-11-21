rm(list=ls())


library(ggplot2)
library(tidyr)
library(dplyr)
library(dartR)
library(cowplot)
library(RColorBrewer)
library(SNPRelate)
library(plotly)
library(readxl)
library(ggthemes)


brewer.pal(9, "BrBG")

######################
#-------DATA---------#
######################



### Metadata
dmeta <- read_excel("../../2022_MicrobiomeQuebecMexico/DATA/aspen_genomic_data/GBS.xlsx")
head(dmeta)
dim(dmeta)

df.geo <- read.csv("05_distances_geo.csv", sep=",", header=T, row.names=1)
head(df.geo)
dim(df.geo)

df.genet <- read.csv("05_distances_genetic.csv", sep=",", header=T, row.names=1)
head(df.genet)
dim(df.genet)

df.bac <- read.csv("05_distances_microbiome_16s.csv", sep=",", header=T, row.names=1)
head(df.bac)
dim(df.bac)

df.fun <- read.csv("05_distances_microbiome_its.csv", sep=",", header=T, row.names=1)
head(df.fun)
dim(df.fun)


rownames(df.geo) == rownames(df.genet)
rownames(df.geo) == rownames(df.bac)





### Distances from AMOS-001 - QUEBEC

amos.geo <- df.geo %>% dplyr::select(AMOS.001)
amos.genet <- df.genet %>% dplyr::select(AMOS.001)
amos.bac <- df.bac %>% dplyr::select(AMOS.001)
amos.fun <- df.fun %>% dplyr::select(AMOS.001)
amos.fun

amos.tot1 <- merge(amos.geo, amos.genet, by.x = 'row.names', by.y = 'row.names', sort=F, all = T)
amos.tot2 <- merge(amos.tot1, amos.bac, by.x = 'Row.names', by.y = 'row.names', sort=F, all = T)
amos.tot <- merge(amos.tot2, amos.fun, by.x = 'Row.names', by.y = 'row.names', sort=F, all = T)
head(amos.tot)

colnames(amos.tot) <- c("Sample","Geo","Genetic","Bacteria","Fungi")
dim(amos.tot)
head(amos.tot)
dm.tot <- merge(amos.tot, dmeta, by.x = "Sample", by.y = "Genotype", sort=F, all = T)
tail(dm.tot)

# Filling out some missing info
dm.tot[dm.tot$Sample == 'STET-001','Location'] <- "STET" 
dm.tot[dm.tot$Sample == 'STET-002','Location'] <- "STET"
dm.tot[dm.tot$Sample == 'STET-001','Region'] <- "Quebec" 
dm.tot[dm.tot$Sample == 'STET-002','Region'] <- "Quebec"
dm.tot[dm.tot$Sample == 'STET-001','Sample.code'] <- 'STET-001'
dm.tot[dm.tot$Sample == 'STET-002','Sample.code'] <- 'STET-002'

head(dm.tot)
tail(dm.tot)
mean.dist <- dm.tot %>% group_by(Location) %>% dplyr::summarize(geo_dist = mean(Geo, na.rm = TRUE))
head(mean.dist)

dm.loc <- merge(dm.tot, mean.dist, by = "Location", sort=F)
head(dm.loc)
dm.qc <- dm.loc %>% filter(Region == "Quebec")
dm.qc


### Distances from SAAA.2 - MEXICO

head(df.geo)

amos.geo <- df.geo %>% dplyr::select(SAAA.2)
amos.genet <- df.genet %>% dplyr::select(SAAA.2)
amos.bac <- df.bac %>% dplyr::select(SAAA.2)
amos.fun <- df.fun %>% dplyr::select(SAAA.2)
amos.fun

amos.tot1 <- merge(amos.geo, amos.genet, by.x = 'row.names', by.y = 'row.names', sort=F, all = T)
amos.tot2 <- merge(amos.tot1, amos.bac, by.x = 'Row.names', by.y = 'row.names', sort=F, all = T)
amos.tot <- merge(amos.tot2, amos.fun, by.x = 'Row.names', by.y = 'row.names', sort=F, all = T)
head(amos.tot)
colnames(amos.tot) <- c("Sample","Geo","Genetic","Bacteria","Fungi")
dim(amos.tot)
head(amos.tot)
dm.tot <- merge(amos.tot, dmeta, by.x = "Sample", by.y = "Genotype", sort=F, all = T)
tail(dm.tot)

#  Filling out some missing info
dm.tot[dm.tot$Sample == 'STET-001','Location'] <- "STET" 
dm.tot[dm.tot$Sample == 'STET-002','Location'] <- "STET"
dm.tot[dm.tot$Sample == 'STET-001','Region'] <- "Quebec" 
dm.tot[dm.tot$Sample == 'STET-002','Region'] <- "Quebec"
dm.tot[dm.tot$Sample == 'STET-001','Sample.code'] <- 'STET-001'
dm.tot[dm.tot$Sample == 'STET-002','Sample.code'] <- 'STET-002'
head(dm.tot)

mean.dist <- dm.tot %>% group_by(Location) %>% dplyr::summarize(geo_dist = mean(Geo, na.rm = TRUE))
head(mean.dist)

dm.loc <- merge(dm.tot, mean.dist, by = "Location", sort=F)
head(dm.loc)
dm.mx <- dm.loc %>% filter(Region == "Mexico")
head(dm.mx)

### Combine

dm.out <- rbind(dm.qc, dm.mx)
dim(dm.out)
head(dm.out)





##############################################

#          Genetic with microbiome

##############################################

dm.out.flt <- dm.out %>% filter(Genetic > 0)
head(dm.out.flt)


tab_blue <- "#1f77b4"
tab_orange <- "#ff7f0e"
head(dm.out.flt)

p5 <- ggplot(dm.out.flt) +
  geom_point(data=dm.out.flt, aes(x = Genetic, y = Bacteria, color = Region),size=3) +
  geom_smooth(data=dm.out.flt, aes(x = Genetic, y = Bacteria, color = Region), method="lm", se=FALSE) +
  scale_color_manual(values=c(tab_orange, tab_blue), name="Country") +

  geom_text(data=data.frame(x = 150, y = 0.3, label = "Mexico: rho=0.1, p=0.2"), aes(x=x, y=y, label=label), size=6, color = tab_orange, hjust=0) +
  geom_text(data=data.frame(x = 150, y = 0.2, label = "Canada: rho=-0.03, p=0.59"), aes(x=x, y=y, label=label), size=6, color = tab_blue, hjust=0) +
  
  labs(x= "Genetic distance", y = "Bacterial distance") +
  scale_x_continuous(limits=c(149,180)) +
  theme(panel.background = element_rect(fill=NA, color="grey20"),
        panel.grid = element_blank(),
        legend.position = "none",
        #legend.title = element_text(size=14),
        #legend.text = element_text(size=14),
        axis.text = element_text(size=14),
        axis.title = element_text(size=20))

p5

p6 <- ggplot(dm.out.flt) +
  geom_point(data=dm.out.flt, aes(x = Genetic, y = Fungi, color = Region),size=3) +
  geom_smooth(data=dm.out.flt, aes(x = Genetic, y = Fungi, color = Region), method="lm", se=FALSE) +
  scale_color_manual(values=c(tab_orange, tab_blue), name="Country") +

  geom_text(data=data.frame(x = 150, y = 0.3, label = "Mexico: rho=0.02, p=0.41"), aes(x=x, y=y, label=label), size=6, color = tab_orange, hjust=0) +
  geom_text(data=data.frame(x = 150, y = 0.2, label = "Canada: rho=-0.01, p=0.54"), aes(x=x, y=y, label=label), size=6, color = tab_blue, hjust=0) +
  
  labs(x= "Genetic distance", y = "Fungal distance") +
  scale_x_continuous(limits=c(149,180)) +
  theme(panel.background = element_rect(fill=NA, color="grey20"),
        panel.grid = element_blank(),
        legend.position = "none",
        #legend.title = element_text(size=14),
        #legend.text = element_text(size=14),
        axis.text = element_text(size=14),
        axis.title = element_text(size=20))
p6


plot_grid(p5, p6, ncol=1)

png("07_compare_distances_16S.png", w = 1400, h=1200, res=300)
p5
dev.off()

png("07_compare_distances_ITS.png", w = 1400, h=1200, res=300)
p6
dev.off()



