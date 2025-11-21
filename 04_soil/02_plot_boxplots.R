rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(cowplot)


#########################
#--------DATA-----------#
#########################


dm <- read.csv("soil_01_data.tsv", sep="\t", header=T)
dt_ <- read.csv("soil_01_data_country_tests.tsv", sep="\t", header=T)
dt <- dt_ %>% filter(EVAL == "OK")

head(dm)
head(dt)

dm %>% filter(metric == "P") %>% arrange(-value) %>% head()
dm[(dm$metric == "P" & dm$value > 400), "value"] <- NaN
dm %>% filter(metric == "P") %>% arrange(-value) %>% head()

# Test correction
head(dt)
pvals <- dt %>% filter(statistic %in% c("p_GR1","p_GR2","p_GR3"))  %>% dplyr::select(statistic, metric, V1)
pvals$adj.p <- p.adjust(pvals$V1, method = "BH")
head(pvals)
pvals_spread <- pvals %>% dplyr::select(statistic, metric, adj.p) %>% spread(key=statistic,adj.p) %>% 
  dplyr::rename(adj_GR1 = p_GR1, adj_GR2 = p_GR2, adj_GR3 = p_GR3)
pvals_spread

# Formatting test tables
dt_spread <- dt %>% spread(key=statistic, value=V1)
dt2 <- merge(dt_spread, pvals_spread, by = c("metric"))

dt2$s1 <- ifelse(dt2$adj_GR1 < 0.05, "*", "ns.")
dt2$s2 <- ifelse(dt2$adj_GR1 < 0.01, "**", dt2$s1)
dt2$sign_GR1 <- ifelse(dt2$adj_GR1 < 0.001, "***", dt2$s2)

dt2$s3 <- ifelse(dt2$adj_GR2 < 0.05, "*", "ns.")
dt2$s4 <- ifelse(dt2$adj_GR2 < 0.01, "**", dt2$s3)
dt2$sign_GR2 <- ifelse(dt2$adj_GR2 < 0.001, "***", dt2$s4)

dt2$s5 <- ifelse(dt2$adj_GR3 < 0.05, "*", "ns.")
dt2$s6 <- ifelse(dt2$adj_GR3 < 0.01, "**", dt2$s5)
dt2$sign_GR3 <- ifelse(dt2$adj_GR3 < 0.001, "***", dt2$s6)


dt2

dg <- dm %>% dplyr::group_by(metric) %>% dplyr::summarize(max_val = max(value, na.rm=T), min_val = min(value, na.rm=T),
                                                          range=abs(min(value, na.rm=T)-max(value, na.rm=T))*0.9)
dg
dg$offset <- dg$max_val + dg$range
dtt <- merge(dt2, dg, by="metric", sort=F)
dtt




#########################
#--------PLOT-----------#
#########################

metrics <- c("C_total","N_total","C_N", "CEC", "pH_H2O",
             "P","K","Ca","Mg","Mn",
             "Al","Fe","Na")

metrics.labels <- c("Total C [%]","Total N [%]","C:N", "CEC [cmol(+)/kg]", "pH H2O",
             "P [mg/kg]","K [cmol(+)/kg]","Ca [cmol(+)/kg]","Mg [cmol(+)/kg]","Mn [cmol(+)/kg]",
             "Al [cmol(+)/kg]","Fe [cmol(+)/kg]","Na [cmol(+)/kg]")

head(dm)

dm$Site <- factor(dm$Site, levels = c("AMOS","STFE","STET","ESSI","FORE","Santiago","FLOR1","FP2"))
dm$Population <- factor(dm$Country, levels = c("Canada","Mexico"))

dm$Metric <- factor(dm$metric, levels = metrics,
                    labels = metrics.labels)

dtt$Metric <- factor(dtt$metric, levels = metrics,
                     labels = metrics.labels)


head(dm)
head(dtt)



p1C <- ggplot(dtt) +
  geom_boxplot(data=dm, aes(x = Site, y = value, fill = Population), outlier.shape = NA) +
  geom_point(data=dm, aes(x = Site, y = value), color="black",
             pch=21, size=2, position=position_jitter(width=0.1)) +
  
  scale_fill_manual(values = c("white","grey"))  + 
  
  geom_vline(xintercept = c(2.5, 5.5), linetype="dashed", color="grey") +
  
  geom_text(data=dtt, aes(x = 0.5, y = offset, label = ifelse(MOD != "GLS",
                                                                        paste0("Group: R2=",round(R2_GR1,digits=2)," ",sign_GR1,
                                                                               "\nRegion: R2=",round(R2_GR2,digits=2)," ",sign_GR2,
                                                                               "\nSite: R2=",round(R2_GR3,digits=2)," ",sign_GR3),
                                                                        paste0("Group: gls ",sign_GR1,
                                                                               "\nRegion: gls ",sign_GR2,
                                                                               "\nSite: gls ",sign_GR3))),
            hjust=0, vjust=1, size=4) +
  
  

  labs(x = "Site") +
  #ggtitle("Population") + 
  #scale_y_continuous(limits=c(-2.3,-0.3)) +
  #facet_grid(Metric~., scales="free_y") +
  facet_wrap(.~Metric, scales="free_y") +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(size=16),
        plot.title = element_text(size=20, face="bold", hjust=0.5),
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "grey30"), 
        #panel.border = element_blank(), 
        axis.ticks = element_blank(), legend.key = element_blank(),
        axis.text.x = element_text(angle=90, hjust=1,vjust=0.5, size=18),
        axis.text.y = element_text(size=14),
        legend.position="none",
        strip.background = element_blank(),
        strip.text = element_text(size=16)
  )
p1C







png('soil_03_plot_main_country.png', w=2600, h=2800, res=300)
p1C
dev.off()






