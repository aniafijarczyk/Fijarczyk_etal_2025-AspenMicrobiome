rm(list=ls())


library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(cowplot)



dm_1 <- read.csv("01_alpha_16S.tsv", sep="\t", header=T)
dm_1$log_gini <- log(dm_1$diversity_gini_simpson + 0.01)
dm_1$marker = "16S"
dm_1 <- dm_1 %>% dplyr::select(long.name,Sample.type,Population,Cluster,marker,chao1,log_gini)
dt_1 <- read.csv("01_alpha_16S_ALL_TEST_RESULTS.tsv", sep="\t", header=T)

dm_2 <- read.csv("01_alpha_ITS.tsv", sep="\t", header=T)
dm_2$log_gini <- log(dm_2$diversity_gini_simpson + 0.01)
dm_2$marker <- "ITS"
dm_2 <- dm_2 %>% dplyr::select(long.name,Sample.type,Population,Cluster,marker,chao1,log_gini)
dt_2 <- read.csv("01_alpha_ITS_ALL_TEST_RESULTS.tsv", sep="\t", header=T)

dm_ <- rbind(dm_1, dm_2)
head(dm_)

dm <- dm_ %>% filter(Population %in% c("NENA","MX") & Sample.type %in% "roots") %>% 
  gather(key = "metric", value="value", c(chao1, log_gini))
head(dm)

dt_ <- rbind(dt_1, dt_2)
dt <- dt_ %>% filter(EVAL == "OK" & source %in% "roots")
head(dt)



# Test correction
head(dt)
pvals <- dt %>% filter(statistic %in% c("p_GR1","p_GR2"))  %>% dplyr::select(statistic, metric, marker, V1)
head(pvals)
pvals$adj.p <- p.adjust(pvals$V1, method = "BH")
pvals %>% head()
pvals %>% dplyr::select(statistic, metric, marker, adj.p) %>% spread(key=statistic,adj.p)
pvals_spread <- pvals %>% dplyr::select(statistic, metric, marker, adj.p) %>% spread(key=statistic,adj.p) %>% dplyr::rename(adj_GR1 = p_GR1, adj_GR2 = p_GR2)
pvals_spread




# Formatting test tables
dt_spread <- dt %>% spread(key=statistic, value=V1)
dt2 <- merge(dt_spread, pvals_spread, by = c("metric","marker"))
dt2

dt2$s1 <- ifelse(dt2$adj_GR1 < 0.05, "*", "ns.")
dt2$s2 <- ifelse(dt2$adj_GR1 < 0.01, "**", dt2$s1)
dt2$sign_GR1 <- ifelse(dt2$adj_GR1 < 0.001, "***", dt2$s2)

dt2$s3 <- ifelse(dt2$adj_GR2 < 0.05, "*", "ns.")
dt2$s4 <- ifelse(dt2$adj_GR2 < 0.01, "**", dt2$s3)
dt2$sign_GR2 <- ifelse(dt2$adj_GR2 < 0.001, "***", dt2$s4)
dt2

dg <- dm %>% dplyr::group_by(metric, marker) %>% dplyr::summarize(max_val = max(value, na.rm=TRUE), min_val = min(value, na.rm=TRUE),
                                                    range=abs(max(value, na.rm=TRUE)-min(value, na.rm=TRUE))*0.2)
dg
dg$offset <- dg$max_val + dg$range
dg$manual_offset <- c(650,170, 0.008,0.15)
dg
dt3 <- merge(dt2, dg, by=c("metric","marker"), sort=F)
head(dt3)

dt3





#########################
#--------PLOT-----------#
#########################


head(dm)
head(dt3)

unique(dm$marker)
dm$Loc <- factor(dm$Cluster, levels = c("AB","FORE","FP","SP"), labels = c("AB","FORE","FP2","Santiago"))
dm$Metric <- factor(dm$metric, levels = c("log_gini","chao1"),
                     labels = c("log Gini-\nSimpson","Richness"))

dt3$Metric <- factor(dt3$metric, levels = c("log_gini","chao1"),
                      labels = c("log Gini-\nSimpson","Richness"))


head(dm)
head(dt3)


##################### ##### 16S

dm_16 <- dm %>% filter(marker == "16S")
dtt1_16 <- dt3 %>% filter(marker == "16S")


p1C <- ggplot(dm_16) +
  geom_boxplot(data=dm_16, aes(x = Loc, y = value, fill=Population), outlier.shape = NA) +
  geom_point(data=dm_16, aes(x = Loc, y = value), color="black",
             pch=21, size=5, position=position_jitter(width=0.1)) +
  scale_fill_manual(values = c("grey","white"))  + 
  
  geom_vline(xintercept = c(2.5, 5.5), linetype="dashed", color="grey") +
  
  geom_text(data=dtt1_16, aes(x = 0.5, y = manual_offset, label = ifelse(MOD %in% c("ANOVA","ANOVA_transformed"),
                                                                        paste0("Group: R2=",round(R2_GR1,2)," ",sign_GR1,
                                                                               "\nSite: R2=",round(R2_GR2,2)," ",sign_GR2),
                                                                        paste0("Group: gls ",sign_GR1,
                                                                               "\nSite: gls ",sign_GR2))),
            hjust=0, vjust=1, size=5) +
  
  #ggtitle("Bacteria/archaea") +
  labs(x = "Site") +
  #scale_y_continuous(limits=c()) +
  facet_wrap(.~Metric, scales="free_y", ncol=2) +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(size=20),
        #plot.title = element_text(size=20, face="bold", hjust=0.5),
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "grey30"), 
        #panel.border = element_blank(), 
        axis.ticks = element_blank(), legend.key = element_blank(),
        axis.text.x = element_text(angle=90, hjust=1,vjust=0.5, size=18),
        axis.text.y = element_text(size=14),
        legend.position="none",
        strip.background = element_blank(),
        strip.text.y = element_blank(),
        strip.text.x = element_text(size=20),
  )
p1C



png('04_plot_markers_roots_16S.png', w=1800, h=1400, res=300)
p1C
dev.off()







##################### ##### ITS

dm_its <- dm %>% filter(marker == "ITS")
dtt1_its <- dt3 %>% filter(marker == "ITS")


p1A <- ggplot(dm_its) +
  geom_boxplot(data=dm_its, aes(x = Loc, y = value, fill=Population), outlier.shape = NA) +
  geom_point(data=dm_its, aes(x = Loc, y = value), color="black",
             pch=21, size=5, position=position_jitter(width=0.1)) +
  scale_fill_manual(values = c("grey","white"))  + 
  
  geom_vline(xintercept = c(2.5, 5.5), linetype="dashed", color="grey") +
  
  geom_text(data=dtt1_its, aes(x = 0.5, y = manual_offset, label = ifelse(MOD %in% c("ANOVA","ANOVA_transformed"),
                                                                         paste0("Group: R2=",round(R2_GR1,2)," ",sign_GR1,
                                                                                "\nSite: R2=",round(R2_GR2,2)," ",sign_GR2),
                                                                         paste0("Group: gls ",sign_GR1,
                                                                                "\nSite: gls ",sign_GR2))),
            hjust=0, vjust=1, size=5) +
  
  #ggtitle("Bacteria/archaea") +
  labs(x = "Site") +
  #scale_y_continuous(limits=c()) +
  facet_wrap(.~Metric, scales="free_y", ncol=2) +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(size=20),
        #plot.title = element_text(size=20, face="bold", hjust=0.5),
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "grey30"), 
        #panel.border = element_blank(), 
        axis.ticks = element_blank(), legend.key = element_blank(),
        axis.text.x = element_text(angle=90, hjust=1,vjust=0.5, size=18),
        axis.text.y = element_text(size=14),
        legend.position="none",
        strip.background = element_blank(),
        strip.text.y = element_blank(),
        strip.text.x = element_text(size=20),
  )
p1A



png('04_plot_markers_roots_ITS.png', w=1800, h=1400, res=300)
p1A
dev.off()


