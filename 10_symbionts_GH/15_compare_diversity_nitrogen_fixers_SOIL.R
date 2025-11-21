rm(list=ls())

library(openxlsx)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(chemometrics)
library(DAAG)
library(phyloseq)
library(microbiome)
library(cowplot)

library(MASS)
library(nlme)
library(emmeans) # to get Tukey test on gls model

source("run_NESTED_ANOVA.R")



#########################
#--------DATA-----------#
#########################



# Raw rarefied abundances, with a few outlier samples removed
load("../../2022_MicrobiomeAspenGreenhouse/01_data_prep/AVS_filtered_phyloseq_rarefied.RData")
sample_data(pseq.16s) %>% head()
otu_table(pseq.16s)
pseq.16s

tab_blue <- "#1f77b4"
tab_orange <- "#ff7f0e"
ametrics <- c('diversity_shannon', 'chao1', 'dominance_core_abundance','evenness_pielou')
site_order <- c("AB","FORE","FP","SP")

### Subsetting samples
dmeta <- read.csv("../../2022_MicrobiomeAspenGreenhouse/01_data_prep/dataPrep_01_overview_16s_META.tsv", sep="\t", header = T)
head(dmeta)
dmeta.soil <- dmeta %>% dplyr::filter(Sample.type == "bulk soil") %>% dplyr::filter(Cluster != "neg control extraction")
dim(dmeta.soil)
head(dmeta.soil)

head(sample_data(pseq.16s))
rseq <- subset_samples(pseq.16s, long.name %in% dmeta.soil$long.name)
pseq.16s
rseq


#####################################################
#--------     Subsetting  NF ASVs        -----------#
#####################################################

asv <- read.csv("03_nitro_read_NF_ASVs.txt", sep="\t", header=T)
asv$ASV
length(asv$ASV)
rseq.fams <- prune_taxa(asv$ASV, rseq)
rseq.fams



#============================================#
#             Alpha diversity                #
#============================================#


metrics <- c('diversity_shannon', 'chao1', 'dominance_core_abundance','evenness_pielou', 'diversity_gini_simpson')


meta <- data.frame(sample_data(rseq.fams))
alpha <- microbiome::alpha(rseq.fams, index = metrics)
dalpha <- merge(meta, alpha, by="row.names", sort=FALSE)
family <- c("NF")
dalpha$Family <- family
head(dalpha)

dalpha %>% head()
dalpha$diversity_gini_simpson[is.na(dalpha$diversity_gini_simpson)] <- 0

write.table(dalpha, file="nitro_05_alpha_soil.tsv", sep='\t', row.names = FALSE, col.names = TRUE, quote = FALSE, append=FALSE)






################
### BOXPLOTS ###
################

metrics
dg.dalpha <- dalpha %>% gather(key = 'measure', value = 'value', all_of(metrics))
head(dg.dalpha)
unique(dg.dalpha$Family)

dg.dalpha$measure <- factor(dg.dalpha$measure, levels = metrics)
dg.dalpha$Cluster <- factor(dg.dalpha$Cluster, levels = site_order)


p5 <- ggplot(dg.dalpha) + aes(y = value, x = Cluster, fill = Population) +
  geom_boxplot() +
  facet_wrap(Family~measure, scales = "free_y", ncol = 5)
p5


#png("nitro_05_alpha_nf_diversity.png", w=1800, h=1800, res=150)
#p5
#dev.off()



#####################
#-----  TESTS  -----#
#####################


###############
### TESTS of difference between countries

head(dalpha)

# Transform CLR so it's positive

dalpha$log_gini <- log(dalpha$diversity_gini_simpson)
min_stat <- min(dalpha$log_gini)
dalpha$log_gini_trans <- log(dalpha$diversity_gini_simpson) + abs(min_stat) + 0.01
hist(dalpha$diversity_gini_simpson)
hist(dalpha$log_gini)
plot(x = dalpha$diversity_gini_simpson, y = dalpha$log_gini_trans)

df.models <- run_nested_models(dalpha, "log_gini_trans", "Population", "Cluster")
df.models

write.table(df.models, file = "nitro_05_alpha_soil_tests.tsv", sep="\t", row.names = F, col.names = T, append=F, quote=F)





################
### PLOT



dalpha$Cluster <- factor(dalpha$Cluster, levels = c("AB","FORE","FP","SP"))
dalpha$Population <- factor(dalpha$Population, levels = c("NENA","MX"))
head(df.models)
head(dalpha)

dLM <- df.models %>% filter(MOD == "ANOVA")
dLM
#dLM[dLM["statistic"]=="p_GR1","V1"]

p1C <- ggplot(dalpha) +
  geom_boxplot(data=dalpha, aes(x = Cluster, y = log_gini, fill = Population), outlier.shape = NA) +
  geom_point(data=dalpha, aes(x = Cluster, y = log_gini), color="grey40",
             pch=21, size=6, position=position_jitter(width=0.1)) +
  
  annotate("text", x = 0.5, y = 0.05, 
           label = paste0("Pop: p=",signif(dLM[dLM["statistic"]=="p_GR1","V1"],digits=2),
                          "\nSite: p=",signif(dLM[dLM["statistic"]=="p_GR2","V1"],digits=2)),
           hjust=0, size=6) +
  
  
  scale_fill_manual(values = c("white","grey"))  + 
  labs(y = "Log Gini-Simpson") +
  #ggtitle("Population") + 
  scale_y_continuous(limits=c(-0.6,0.1)) +
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


#png('nitro_05_alpha_soil.png', w=800, h=1200, res=300)
#p1C
#dev.off()







