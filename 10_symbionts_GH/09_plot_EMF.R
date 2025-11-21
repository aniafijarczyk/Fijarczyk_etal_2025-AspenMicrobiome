rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(cowplot)


#########################
#--------DATA-----------#
#########################


dm <- read.csv("guilds_05_combine_inputs.tsv", sep="\t", header=T)
dt <- read.csv("guilds_05_combine_tests.tsv", sep="\t", header=T)

head(dm)
head(dt)


# Test correction
head(dt)
pvals <- dt %>% filter(statistic %in% c("p_GR1","p_GR2"))  %>% dplyr::select(statistic, Sample.type, metric, V1)
pvals$adj.p <- p.adjust(pvals$V1, method = "BH")
pvals_spread <- pvals %>% dplyr::select(statistic, Sample.type, metric, adj.p) %>% tidyr::spread(key=statistic,adj.p) %>% 
  dplyr::rename(adj_GR1 = p_GR1, adj_GR2 = p_GR2)
pvals_spread

# Formatting test tables
dt_spread <- dt %>% spread(key=statistic, value=V1)
dt2 <- merge(dt_spread, pvals_spread, by = c("Sample.type","metric"))

dt2$s1 <- ifelse(dt2$adj_GR1 < 0.05, "*", "ns.")
dt2$s2 <- ifelse(dt2$adj_GR1 < 0.01, "**", dt2$s1)
dt2$sign_GR1 <- ifelse(dt2$adj_GR1 < 0.001, "***", dt2$s2)

dt2$s3 <- ifelse(dt2$adj_GR2 < 0.05, "*", "ns.")
dt2$s4 <- ifelse(dt2$adj_GR2 < 0.01, "**", dt2$s3)
dt2$sign_GR2 <- ifelse(dt2$adj_GR2 < 0.001, "***", dt2$s4)
dt2

dg <- dm %>% dplyr::group_by(metric) %>% dplyr::summarize(max_val = max(value), min_val = min(value), range=abs(min(value)-max(value))*0.2)
dg
dg$offset <- dg$max_val + dg$range
dtt <- merge(dt2, dg, by="metric", sort=F)
dtt

write.table(dtt, file="guilds_06_plot_main.tsv",sep="\t",row.names = F, col.names = T,quote=F,append=F)






#########################
#--------PLOT-----------#
#########################

dm$Site <- factor(dm$Cluster, levels = c("AB","FORE","FP","SP"), labels = c("AB","FORE","FP2","Santiago"))
dm$Population <- factor(dm$Population, levels = c("NENA","MX"))
dm$Source <- factor(dm$Sample.type, levels = c("bulk soil","rhizosphere","roots"),
                    labels = c("Potting mix","Rhizosphere","Roots"))
dm$Metric <- factor(dm$metric, levels = c("CLR.abd","log_gini","chao1"),
                    labels = c("CLR\nabundance","log Gini-\nSimpson","Richness"))


dtt$Source <- factor(dtt$Sample.type, levels = c("bulk soil","rhizosphere","roots"),
                    labels = c("Potting mix","Rhizosphere","Roots"))
dtt$Metric <- factor(dtt$metric, levels = c("CLR.abd","log_gini","chao1"),
                    labels = c("CLR\nabundance","log Gini-\nSimpson","Richness"))


head(dm)
head(dtt)



p1C <- ggplot(dtt) +
  geom_boxplot(data=dm, aes(x = Site, y = value, fill = Population), outlier.shape = NA) +
  geom_point(data=dm, aes(x = Site, y = value), color="black",
             pch=21, size=5, position=position_jitter(width=0.1)) +
  
  #annotate("text", x = 0.5, y = -0.5, 
  #         label = paste0("Pop: p=",signif(dLM[dLM["statistic"]=="p_GR1","V1"],digits=2),
  #                        "\nSite: p=",signif(dLM[dLM["statistic"]=="p_GR2","V1"],digits=2)),
  #         hjust=0, size=6) +
  
  
  geom_text(data=dtt, aes(x = 0.5, y = offset+5, label = ifelse(MOD %in% c("ANOVA","ANOVA_transformed"),
                                                                        paste0("Group: R2=",round(R2_GR1,2)," ",sign_GR1,
                                                                               "\nSite: R2=",round(R2_GR2,2)," ",sign_GR2),
                                                                        paste0("Group: gls ",sign_GR1,
                                                                               "\nSite: gls ",sign_GR2))),
            hjust=0, vjust=1, size=4) +
  
  
  scale_fill_manual(values = c("white","grey"))  + 
  labs(x = "Site") +
  #ggtitle("Population") + 
  #scale_y_continuous(limits=c(-2.3,-0.3)) +
  facet_grid(Metric~Source, scales="free_y") +
  theme(axis.title.y = element_blank(),
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
        strip.text = element_text(size=20)
  )
p1C



png('guilds_06_plot_main.png', w=1800, h=2000, res=300)
p1C
dev.off()






