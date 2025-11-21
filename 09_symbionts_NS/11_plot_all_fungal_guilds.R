rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(cowplot)


#########################
#--------DATA-----------#
#########################


dm <- read.csv("guilds_05_combine_guilds_inputs.tsv", sep="\t", header=T)
dt <- read.csv("guilds_05_combine_guilds_tests.tsv", sep="\t", header=T)

head(dm)
head(dt)


# Test correction
head(dt)
pvals <- dt %>% filter(statistic %in% c("p_GR1","p_GR2"))  %>% dplyr::select(statistic, Phylum, V1)
pvals
pvals$adj.p <- p.adjust(pvals$V1, method = "BH")
pvals_spread <- pvals %>% dplyr::select(statistic, Phylum, adj.p) %>% spread(key=statistic,adj.p) %>% rename(adj_GR1 = p_GR1, adj_GR2 = p_GR2)
pvals_spread

# Formatting test tables
dt_spread <- dt %>% spread(key=statistic, value=V1)
dt2 <- merge(dt_spread, pvals_spread, by = c("Phylum"))
head(dt2)

dt2$s1 <- ifelse(dt2$adj_GR1 < 0.05, "*", "ns.")
dt2$s2 <- ifelse(dt2$adj_GR1 < 0.01, "**", dt2$s1)
dt2$sign_GR1 <- ifelse(dt2$adj_GR1 < 0.001, "***", dt2$s2)

dt2$s3 <- ifelse(dt2$adj_GR2 < 0.05, "*", "ns.")
dt2$s4 <- ifelse(dt2$adj_GR2 < 0.01, "**", dt2$s3)
dt2$sign_GR2 <- ifelse(dt2$adj_GR2 < 0.001, "***", dt2$s4)
dt2

dg <- dm %>% group_by(metric, Phylum) %>% summarize(max_val = max(value), min_val = min(value), range=abs(min(value)-max(value))*0.5)

dg
dg$offset <- dg$max_val + dg$range
#dg$offset <- dg$max_val + 2
dtt <- merge(dt2, dg, by=c("metric","Phylum"), sort=F)
dtt




#########################
#--------PLOT-----------#
#########################
head(dm)
unique(dm$Phylum)

dm$Site <- factor(dm$site, levels = c("AMOS","STFE","STET","ESSI","FORE","Santiago","FLOR1","FP2"))
dm$Population <- factor(dm$Country, levels = c("Canada","Mexico"))


dm$Metric <- factor(dm$metric, levels = c("CLR.abd"),
                    labels = c("CLR\nabundance"))
dm$Label <- factor(dm$Phylum, levels=c("Ectomycorrhizal","Undefined Saprotroph","Ericoid Mycorrhizal","Endophyte",
                                                         "Plant Pathogen","Wood Saprotroph","Animal Pathogen","Fungal Parasite",
                                       "Arbuscular Mycorrhizal", "NF"),
                            labels=c("Ecto-\nmycorrhizal","Undefined\nsaprotroph","Ericoid\nmycorrhizal","\nEndophyte",
                                     "Plant\npathogen","Wood\nsaprotroph","Animal\npathogen","Fungal\nparasite",
                                     "Arbuscular\nMycorrhizal", "Nitrogen-fixing\nbacteria"))

head(dtt)
dtt$Metric <- factor(dtt$metric, levels = c("CLR.abd"),
                     labels = c("CLR\nabundance"))
dtt$Label <- factor(dtt$Phylum, levels=c("Ectomycorrhizal","Undefined Saprotroph","Ericoid Mycorrhizal","Endophyte",
                                         "Plant Pathogen","Wood Saprotroph","Animal Pathogen","Fungal Parasite",
                                         "Arbuscular Mycorrhizal", "NF"),
                    labels=c("Ecto-\nmycorrhizal","Undefined\nsaprotroph","Ericoid\nmycorrhizal","\nEndophyte",
                             "Plant\npathogen","Wood\nsaprotroph","Animal\npathogen","Fungal\nparasite",
                             "Arbuscular\nMycorrhizal", "Nitrogen-fixing\nbacteria"))


head(dm)
head(dtt)

p1C <- ggplot(dtt) +
  geom_boxplot(data=dm, aes(x = Site, y = value, fill = Population), outlier.shape = NA) +
  geom_point(data=dm, aes(x = Site, y = value), color="black",
             pch=21, size=3, position=position_jitter(width=0.1)) +
  scale_fill_manual(values = c("white","grey"))  +
  geom_vline(xintercept = c(2.5, 5.5), linetype="dashed", color="grey") +
  
  geom_text(data=dtt, aes(x = 0.5, y = offset, label = ifelse(MOD != "GLS",
                                                                paste0("Region: R2=",round(R2_GR1,2)," ",sign_GR1,
                                                              "\nSite: R2=",round(R2_GR2,2)," ",sign_GR2),
                                                              paste0("Region: gls ",sign_GR1,
                                                                     "\nSite: gls ",sign_GR2))),
            hjust=0, vjust=1) +
  
 
  labs(x = "Site", y = "CLR abundance") +
  #ggtitle("Population") + 
  #scale_y_continuous(limits=c(-2.3,-0.3)) +
  facet_wrap(~Label, scales="free_y", nrow=2) +
  theme(axis.title.y = element_text(size=20),
        axis.title.x = element_text(size=20),
        plot.title = element_text(size=20, face="bold", hjust=0.5),
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "grey30"), 
        #panel.border = element_blank(), 
        axis.ticks = element_blank(), legend.key = element_blank(),
        axis.text.x = element_text(angle=90, hjust=1,vjust=0.5, size=18),
        axis.text.y = element_text(size=14),
        legend.position="none",
        strip.background = element_blank(),
        strip.text = element_text(size=18)
  )
p1C



png('guilds_07_plot_all_guilds.png', w=3500, h=2000, res=300)
p1C
dev.off()





