rm(list=ls())

library(ggplot2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)

#==================================#
#               DATA               #
#==================================#

df1 <- read.csv("11_get_shared_groups_16S_site_GENUS.tsv", sep="\t", header=T)
df2 <- read.csv("11_get_shared_groups_16S_site_ASV.tsv", sep="\t", header=T)

head(df1)




#==================================#
#           BARPLOTS               #
#==================================#

##################################################################################
## Try facet_grid
df1 <- read.csv("11_get_shared_groups_16S_site_GENUS.tsv", sep="\t", header=T)
df2 <- read.csv("11_get_shared_groups_16S_site_ASV.tsv", sep="\t", header=T)
df1$type <- "Genus"
df2$type <- "ASV"
d12 <- rbind(df1, df2)
n_max <- d12 %>% group_by(type) %>% dplyr::summarize(max_n = max(n))
head(d12)
head(n_max)
dm <- merge(d12, n_max, by = "type", sort=F)
head(dm)

lab_order <- c("Pan-population","Country-specific","Canada","Mexico",
               "Pan-region","Region-specific","Boreal","Cold-temperate","Warm-temperate",
               "Pan-site","Site-specific","AMOS","STFE","STET","ESSI","FORE","FLOR1","FP2","Santiago")
lab_labels <- c("Shared","Unique","Canada","Mexico",
                "Shared","Unique","Boreal","CTemp","WTemp",
                "Shared","Unique","AMOS","STFE","STET","ESSI","FORE","FLOR1","FP2","Santiago")

unique(df1$group_id)
dm$label <- factor(dm$label, levels=lab_order, labels=lab_labels)
dm$group_id <- factor(dm$group_id, levels = c("Population","Region","Site"))
dm$color_group <- factor(dm$color_group, levels = c("shared","unique","spec"))

p3 <- ggplot(dm) + aes(x = label, y = n, fill = color_group) +
  geom_bar(stat='identity') +
  facet_grid(type~group_id, space="free_x", scales="free") +
  scale_fill_manual(values=c("grey20","#CC79A7","lightgrey")) + 
  labs(x = "Group", y = "Count") + 
  
  geom_text(data=dm[dm$color_group %in% c("shared","unique"),],
            aes(x = label, y = n+(max_n*0.1), label = paste0(round(percent,1), "")),
            size=6) +
  
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size=18),
        axis.text.y = element_text(size=14),
        axis.title = element_text(size=18),
        #strip.background = element_blank(),
        strip.text = element_text(size=20),
        legend.position = "none",
        
  )

p3


png("12_plot_barplots_16S.png", w=2800,h=1800, res=300)
p3
dev.off()






##################################################################################
############### ITS
##################################################################################


df1 <- read.csv("11_get_shared_groups_ITS_site_GENUS.tsv", sep="\t", header=T)
df2 <- read.csv("11_get_shared_groups_ITS_site_ASV.tsv", sep="\t", header=T)
df1$type <- "Genus"
df2$type <- "ASV"
d12 <- rbind(df1, df2)
n_max <- d12 %>% group_by(type) %>% dplyr::summarize(max_n = max(n))
dm <- merge(d12, n_max, by = "type", sort=F)
head(dm)

lab_order <- c("Pan-population","Country-specific","Canada","Mexico",
               "Pan-region","Region-specific","Boreal","Cold-temperate","Warm-temperate",
               "Pan-site","Site-specific","AMOS","STFE","STET","ESSI","FORE","FLOR1","FP2","Santiago")
lab_labels <- c("Shared","Unique","Canada","Mexico",
                "Shared","Unique","Boreal","CTemp","WTemp",
                "Shared","Unique","AMOS","STFE","STET","ESSI","FORE","FLOR1","FP2","Santiago")

unique(df1$group_id)
dm$label <- factor(dm$label, levels=lab_order, labels=lab_labels)
dm$group_id <- factor(dm$group_id, levels = c("Population","Region","Site"))
dm$color_group <- factor(dm$color_group, levels = c("shared","unique","spec"))


p4 <- ggplot(dm) + aes(x = label, y = n, fill = color_group) +
  geom_bar(stat='identity') +
  facet_grid(type~group_id, space="free_x", scales="free") +
  scale_fill_manual(values=c("grey20","#CC79A7","lightgrey")) + 
  labs(x = "Group", y = "Count") + 
  
  geom_text(data=dm[dm$color_group %in% c("shared","unique"),],
            aes(x = label, y = n+(max_n*0.1), label = paste0(round(percent,1), "")),
            size=6) +
  
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size=18),
        axis.text.y = element_text(size=14),
        axis.title = element_text(size=18),
        #strip.background = element_blank(),
        strip.text = element_text(size=20),
        legend.position = "none",
  )

p4


png("12_plot_barplots_ITS.png", w=2800,h=1800, res=300)
p4
dev.off()
