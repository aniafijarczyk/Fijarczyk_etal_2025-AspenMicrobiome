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




#################
#---FILTERING---#
#################

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









######################
#---     16s      ---#
######################



###############
#    PCoA     #
###############

pcoa <- ordinate(pseq.16s.2.rel, "PCoA", "bray")
pcoa$trace
pcoa$values %>% head()
pcoa$vectors %>% head()
pcoa$values$Eigenvalues
pcoa$values$Relative_eig
# Percent variance explained by each axis
variance_explained <- 100 * pcoa$values$Relative_eig
# Check first two axes
variance_explained[1:2]

####### Getting table
data.scores <- as.data.frame(pcoa$vectors[,1:3])
length(rownames(dmeta.16s))
head(dmeta.16s)
sum(rownames(data.scores) == rownames(dmeta.16s))
data.scores$Source = dmeta.16s$Sample.type
data.scores$Group = dmeta.16s$Population
data.scores$Site = dmeta.16s$Cluster
head(data.scores)





#####################
#   PCA  - sources  #
#####################


head(data.scores)
source_order <- c("roots","rhizosphere","bulk soil")
data.scores$Source_labels <- factor(data.scores$Source, levels = source_order, labels=c("roots","rhizo","potting\nmix"))
head(data.scores)

source_cols <- c("black","grey","#F0E442")

gs <- ggplot(data = data.scores) + 
  geom_point(pch=21, data = data.scores, aes(x = Axis.1, y = Axis.2, fill = Source_labels), size = 4, alpha = 1) + 
  scale_fill_manual(values = source_cols, name="Sample\ntype")  +
  
  ggtitle("\nBacteria/Archaea") + 
  
  labs(x = paste0("PCo 1 (",round(variance_explained[1],1),"%)"),
       y = paste0("PCo 2 (",round(variance_explained[2],1),"%)")) + 
  
  geom_text(data=data.frame("x" = c(min(data.scores$Axis.1)), 
                            "y" = c(max(data.scores$Axis.2+0.4)),
                            "label" = paste0("Sample type: R2=0.19**",
                                             "\nGroup: R2=0.07**",
                                             "\nSite: R2=0.17**")), aes(x = x, y = y, label = label), hjust=0, vjust=1,size=5) +

  theme(axis.title = element_text(size = 14), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(),
        legend.title = element_text(size = 12, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 12, colour = "grey30"),
        plot.title = element_text(size = 14, hjust=0.5, face="bold"),
  )

gs

#png("ordination_04_ordination_16S.source.png", w=1000,h=800, res=300)
#gs
#dev.off()


#####################
### PCA color by Site
#####################

sites_cols <- c("#134b71","#258dd5","#ffbe84","#c15a00")
head(data.scores)
data.scores$Source

df.soil <- data.scores %>% filter(Source == "bulk soil")
df.soil
site_order <- c("AB","FORE","FP","SP")
df.soil$Site <- factor(df.soil$Site, levels = site_order,
                       labels = c("AB","FORE","FP2","Santiago"))
df.soil

gs2 <- ggplot(data = df.soil) + 
  geom_point(pch=21, data = df.soil, aes(x = Axis.1, y = Axis.2, fill = Site), size = 4, alpha = 1) + 
  scale_fill_manual(values = sites_cols, name="Site")  +

  ggtitle("Potting mix\n bacteria/Archaea") + 
  
  labs(x = paste0("PCo 1 (",round(variance_explained[1],1),"%)"),
       y = paste0("PCo 2 (",round(variance_explained[2],1),"%)")) + 
  
  geom_text(data=data.frame("x" = c(min(df.soil$Axis.1)), 
                            "y" = c(max(df.soil$Axis.2)),
                            "label" = paste0("Group: R2=0.12 ns.",
                                             "\nSite: R2=0.35 *")), aes(x = x, y = y+0.2, label = label), hjust=0, vjust=1, size=5) +
  
  #scale_y_continuous(limits=c(min(df.soil$Axis.1),max(df.soil$Axis.1))) +
  
  theme(axis.title = element_text(size = 14), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(),
        legend.title = element_text(size = 12, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 12, colour = "grey30"),
        plot.title = element_text(size = 14, hjust=0.5, face="bold"),
        )

gs2


#png("ordination_04_ordination_16S.soil.png", w=1000,h=800, res=300)
#gs2
#dev.off()


png("ordination_04_ordination_16S.source.soil.png",  w=2200,h=1000, res=300)
plot_grid(gs, gs2, rel_widths = c(1,1))
dev.off()





#####################
### PCA    RHIZO  
#####################

brewer.pal(n=8,"BrBG")
sites_cols <- c("#134b71","#258dd5","#ffbe84","#c15a00")
head(data.scores)
data.scores$Source
df.rhizo <- data.scores %>% filter(Source == 'rhizosphere')
site_order <- c("AB","FORE","FP","SP")
df.rhizo$Site <- factor(df.rhizo$Site, levels = site_order,
                       labels = c("AB","FORE","FP2","Santiago"))
df.rhizo


gs3 <- ggplot(data = df.rhizo) + 
  geom_point(pch=21, data = df.rhizo, aes(x = Axis.1, y = Axis.2, fill = Site), size = 4, alpha = 1) + 
  scale_fill_manual(values = sites_cols, name="Site")  +
  
  ggtitle("Rhizosphere\nbacteria/Archaea") + 
  
  labs(x = paste0("PCo 1 (",round(variance_explained[1],1),"%)"),
       y = paste0("PCo 2 (",round(variance_explained[2],1),"%)")) + 
  
  geom_text(data=data.frame("x" = c(min(df.rhizo$Axis.1)), 
                            "y" = c(max(df.rhizo$Axis.2)+0.3),
                            "label" = paste0("Group: R2=0.13 ns.",
                                             "\nSite: R2=0.35 ns.")), aes(x = x, y = y, label = label), hjust=0,vjust=1, size=5) +

  theme(axis.title = element_text(size = 14), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(),
        legend.title = element_text(size = 12, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 12, colour = "grey30"),
        plot.title = element_text(size = 14, hjust=0.5, face="bold"),
  )

gs3


png("ordination_04_ordination_16S.rhizo.png", w=1100,h=900, res=300)
gs3
dev.off()



#####################
### PCA    ROOTS  
#####################

brewer.pal(n=8,"BrBG")
sites_cols <- c("#134b71","#258dd5","#ffbe84","#c15a00")

head(data.scores)
df.roots <- data.scores %>% filter(Source == 'roots')
site_order <- c("AB","FORE","FP","SP")
df.roots$Site <- factor(df.roots$Site, levels = site_order,
                        labels = c("AB","FORE","FP2","Santiago"))
df.roots


gs4 <- ggplot(data = df.roots) + 
  geom_point(pch=21, data = df.roots, aes(x = Axis.1, y = Axis.2, fill = Site), size = 4, alpha = 1) + 
  scale_fill_manual(values = sites_cols, name="Site")  +
  
  ggtitle("Roots\nbacteria/Archaea") + 
  
  labs(x = paste0("PCo 1 (",round(variance_explained[1],1),"%)"),
       y = paste0("PCo 2 (",round(variance_explained[2],1),"%)")) + 
  
  geom_text(data=data.frame("x" = c(min(df.roots$Axis.1)), 
                            "y" = c(max(df.roots$Axis.2)+0.3),
                            "label" = paste0("Group: R2=0.1 ns.",
                                             "\nSite: R2=0.25 ns.")), aes(x = x, y = y, label = label), hjust=0,vjust=1, size=5) +
  #scale_y_continuous(limits=c(min(df.roots$PC2),max(df.roots$PC2))) +
  
  theme(axis.title = element_text(size = 14), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(),
        legend.title = element_text(size = 12, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 12, colour = "grey30"),
        plot.title = element_text(size = 14, hjust=0.5, face="bold"),
  )

gs4


png("ordination_04_ordination_16S.roots.png", w=1100,h=900, res=300)
gs4
dev.off()





######################
#---     ITS      ---#
######################




###############
#    PCoA     #
###############

pcoa <- ordinate(pseq.its.2.rel, "PCoA", "bray")
pcoa$trace
pcoa$values %>% head()
pcoa$vectors %>% head()
pcoa$values$Eigenvalues
pcoa$values$Relative_eig
# Percent variance explained by each axis
variance_explained <- 100 * pcoa$values$Relative_eig
# Check first two axes
variance_explained[1:2]

####### Getting table
data.scores <- as.data.frame(pcoa$vectors[,1:3])
length(rownames(dmeta.its))
head(dmeta.its)
sum(rownames(data.scores) == rownames(dmeta.its))
data.scores$Source = dmeta.its$Sample.type
data.scores$Group = dmeta.its$Population
data.scores$Site = dmeta.its$Cluster
head(data.scores)







#####################
### PCA  - sources  #
#####################

head(data.scores)
data.scores$Source
source_order <- c("roots","rhizosphere","bulk soil")
data.scores$Source_label <- factor(data.scores$Source, levels = source_order, labels=c("roots","rhizo","potting\nmix"))

source_cols <- c("black","grey","#F0E442")


gs <- ggplot(data = data.scores) + 
  geom_point(pch=21, data = data.scores, aes(x = Axis.1, y = Axis.2, fill = Source_label), size = 4, alpha = 1) + 
  scale_fill_manual(values = source_cols, name="Sample\ntype")  +
  
  ggtitle("\nFungi") + 
  
  labs(x = paste0("PCo 1 (",round(variance_explained[1],1),"%)"),
       y = paste0("PCo 2 (",round(variance_explained[2],1),"%)")) + 
  
  geom_text(data=data.frame("x" = c(min(data.scores$Axis.1)), 
                            "y" = c(max(data.scores$Axis.2)+0.3),
                            "label" = paste0("Sample type: R2=0.1*",
                                             "\nGroup: R2=0.12**",
                                             "\nSite: R2=0.23**")), aes(x = x, y = y, label = label), hjust=0, vjust=1,size=5) +
  
  #scale_y_continuous(limits=c(min(data.scores$PC2),max(data.scores$PC2)+55+1)) +

  theme(axis.title = element_text(size = 14), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(),
        legend.title = element_text(size = 12, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 12, colour = "grey30"),
        plot.title = element_text(size = 14, hjust=0.5, face="bold"),
  )

gs

#png("ordination_04_ordination_ITS.source.png",  w=1000,h=800, res=300)
#gs
#dev.off()


#####################
### PCA color by Site
#####################

brewer.pal(n=8,"BrBG")
sites_cols <- c("#134b71","#258dd5","#ffbe84","#c15a00")
head(data.scores)

df.soil <- data.scores %>% filter(Source == 'bulk soil')
site_order <- c("AB","FORE","FP","SP")
df.soil$Site <- factor(df.soil$Site, levels = site_order,
                       labels = c("AB","FORE","FP2","Santiago"))
df.soil

gs2 <- ggplot(data = df.soil) + 
  geom_point(pch=21, data = df.soil, aes(x = Axis.1, y = Axis.2, fill = Site), size = 4, alpha = 1) + 
  scale_fill_manual(values = sites_cols, name="Site")  +
  
  ggtitle("Potting mix\nfungi") + 
  
  labs(x = paste0("PCo 1 (",round(variance_explained[1],1),"%)"),
       y = paste0("PCo 2 (",round(variance_explained[2],1),"%)")) + 
  
  geom_text(data=data.frame("x" = c(min(df.soil$Axis.1)), 
                            "y" = c(max(df.soil$Axis.2)+0.2),
                            "label" = paste0("Group: R2=0.18*",
                                             "\nSite: R2=0.35 ns.")), aes(x = x, y = y, label = label), hjust=0, vjust=1, size=5) +
  
  #scale_y_continuous(limits=c(min(df.soil$PC2),max(df.soil$PC2)+25+1)) +
  
  theme(axis.title = element_text(size = 14), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(),
        legend.title = element_text(size = 12, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 12, colour = "grey30"),
        plot.title = element_text(size = 14, hjust=0.5, face="bold"),
  )

gs2


#png("ordination_04_ordination_ITS.soil.png", w=550,h=400, res=150)
#gs2
#dev.off()


png("ordination_04_ordination_ITS.source.soil.png",  w=2200,h=1000, res=300)
plot_grid(gs, gs2, rel_widths = c(1,1))
dev.off()






#####################
### PCA    RHIZO  
#####################

brewer.pal(n=8,"BrBG")
sites_cols <- c("#134b71","#258dd5","#ffbe84","#c15a00")

head(data.scores)
data.scores$Source
df.rhizo <- data.scores %>% filter(Source == 'rhizosphere')
site_order <- c("AB","FORE","FP","SP")
df.rhizo$Site <- factor(df.rhizo$Site, levels = site_order,
                        labels = c("AB","FORE","FP2","Santiago"))
df.rhizo


gs3 <- ggplot(data = df.rhizo) + 
  geom_point(pch=21, data = df.rhizo, aes(x = Axis.1, y = Axis.2, fill = Site), size = 4, alpha = 1) + 
  scale_fill_manual(values = sites_cols, name="Site")  +
  
  ggtitle("Rhizosphere\nfungi") + 
  
  labs(x = paste0("PCo 1 (",round(variance_explained[1],1),"%)"),
       y = paste0("PCo 2 (",round(variance_explained[2],1),"%)")) + 
  
  geom_text(data=data.frame("x" = c(min(df.rhizo$Axis.1)), 
                            "y" = c(max(df.rhizo$Axis.2)+0.2),
                            "label" = paste0("Group: R2=0.15 *",
                                             "\nSite: R2=0.33 ns.")), aes(x = x, y = y, label = label),
            hjust=0, vjust=1, size=5) +
  
  #scale_y_continuous(limits=c(min(df.rhizo$Axis.1),max(df.rhizo$PC2)+25+1)) +

  theme(axis.title = element_text(size = 14), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(),
        legend.title = element_text(size = 12, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 12, colour = "grey30"),
        plot.title = element_text(size = 14, hjust=0.5, face="bold"),
  )

gs3


png("ordination_04_ordination_ITS.rhizo.png", w=1100,h=900, res=300)
gs3
dev.off()



#####################
### PCA    ROOTS  
#####################

brewer.pal(n=8,"BrBG")
sites_cols <- c("#134b71","#258dd5","#ffbe84","#c15a00")
head(data.scores)
df.roots <- data.scores %>% filter(Source == 'roots')
site_order <- c("AB","FORE","FP","SP")
df.roots$Site <- factor(df.roots$Site, levels = site_order,
                        labels = c("AB","FORE","FP2","Santiago"))
df.roots

gs4 <- ggplot(data = df.roots) + 
  geom_point(pch=21, data = df.roots, aes(x = Axis.1, y = Axis.2, fill = Site), size = 4, alpha = 1) + 
  scale_fill_manual(values = sites_cols, name="Site")  +
  
  ggtitle("Roots\nfungi") + 
  
  labs(x = paste0("PCo 1 (",round(variance_explained[1],1),"%)"),
       y = paste0("PCo 2 (",round(variance_explained[2],1),"%)")) + 
  
  geom_text(data=data.frame("x" = c(min(df.roots$Axis.1)), 
                            "y" = c(max(df.roots$Axis.2)+0.2),
                            "label" = paste0("Group: R2=0.12 ns.",
                                             "\nSite: R2=0.29 ns.")), aes(x = x, y = y, label = label),
            hjust=0, vjust=1, size=5) +
  
  #scale_y_continuous(limits=c(min(df.roots$PC2),max(df.roots$PC2)+15+1)) +  
  
  theme(axis.title = element_text(size = 14), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(),
        legend.title = element_text(size = 12, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 12, colour = "grey30"),
        plot.title = element_text(size = 14, hjust=0.5, face="bold"),
  )

gs4


png("ordination_04_ordination_ITS.roots.png", w=1100,h=900, res=300)
gs4
dev.off()


#png("ordination_04_ordination_ITS.rhizo.roots.png", w=1000,h=400, res=150)
#plot_grid(gs3, gs4, rel_widths = c(1,1))
#ev.off()

