rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(cowplot)


#########################
#--------DATA-----------#
#########################


dm_ <- read.csv("guilds_05_combine_inputs.tsv", sep="\t", header=T)
dt <- read.csv("guilds_05_combine_tests.tsv", sep="\t", header=T)

head(dm_)
head(dt)

# Formatting table
head(dm_)
#dm <- dm_ %>% gather(key = "Group", value = "loc", c("Country","region"))
dm <- dm_
head(dm)
dm %>% filter(sample == "X018.AMOS.001.S.Soil.ITS")


###########################
### Some outliers to remove
outliers = c("X18S.AMF.036.STFE.007.S.Sol","18S-AMF-036-STFE-007-S-Sol")

# Masking values of diversity that was 0 (NaN)
dm[(dm$metric == "log_gini" & dm$value < -4), "value"] <- NaN
dm[(dm$metric == "log_gini" & dm$sample %in% outliers), "value"] <- NaN

head(dm)


# Test correction
dt1 <- dt %>% filter(Group == "region")
head(dt1)
pvals <- dt1 %>% filter(statistic %in% c("p_GR1","p_GR2"))  %>% dplyr::select(statistic, metric, Group, Phylum, V1)
head(pvals)
pvals$adj.p <- p.adjust(pvals$V1, method = "BH")
pvals %>% head()
pvals_spread <- pvals %>% dplyr::select(statistic, metric, Group, Phylum, adj.p) %>% spread(key=statistic,adj.p) %>% rename(adj_GR1 = p_GR1, adj_GR2 = p_GR2)
pvals_spread

# Formatting test tables
dt_spread <- dt1 %>% spread(key=statistic, value=V1)
dt2 <- merge(dt_spread, pvals_spread, by = c("metric","Group","Phylum"))
dt2

dt2$s1 <- ifelse(dt2$adj_GR1 < 0.05, "*", "ns.")
dt2$s2 <- ifelse(dt2$adj_GR1 < 0.01, "**", dt2$s1)
dt2$sign_GR1 <- ifelse(dt2$adj_GR1 < 0.001, "***", dt2$s2)

dt2$s3 <- ifelse(dt2$adj_GR2 < 0.05, "*", "ns.")
dt2$s4 <- ifelse(dt2$adj_GR2 < 0.01, "**", dt2$s3)
dt2$sign_GR2 <- ifelse(dt2$adj_GR2 < 0.001, "***", dt2$s4)
dt2

dg <- dm %>% group_by(metric, Phylum) %>% summarize(max_val = max(value, na.rm=TRUE), min_val = min(value, na.rm=TRUE),
                                                    range=abs(max(value, na.rm=TRUE)-min(value, na.rm=TRUE))*0.2)
dg
dg$offset <- dg$max_val + dg$range
dg$manual_offset <- c(14,10,1, 180,100,40, 1.5,0.15,2.2)
dg
dtt1 <- merge(dt2, dg, by=c("metric","Phylum"), sort=F)
head(dtt1)

dtt1

head(dm)
write.table(dtt1, file="guilds_08_plot_main.tsv",sep="\t",row.names = F, col.names = T,quote=F,append=F)
head(dtt1)

#########################
#--------PLOT-----------#
#########################

dm1 <- dm
head(dm1)
head(dtt1)

unique(dm1$Phylum)
#dm$site <- factor(dm$site, levels = c("AB","FORE","FP","SP"), labels = c("AB","FORE","FP2","SP"))
dm1$Loc <- factor(dm1$site, levels = c("AMOS","STFE","STET","ESSI","FORE","Santiago","FLOR1","FP2"),
                  labels <- c("AMOS","STFE","STET","ESSI","FORE","Santiago","FLOR1","FP2"))
dm1$Phylum <- factor(dm1$Phylum, levels = c("Ectomycorrhizal","Arbuscular Mycorrhizal","NF"),
                    labels = c("EMF","AMF","Nitrogen-fixing"))
dm1$Metric <- factor(dm1$metric, levels = c("CLR.abd","log_gini","chao1"),
                    labels = c("CLR\nabundance","log Gini-\nSimpson","Richness"))


dtt1$Phylum <- factor(dtt1$Phylum, levels = c("Ectomycorrhizal","Arbuscular Mycorrhizal","NF"),
                      labels = c("EMF","AMF","Nitrogen-fixing"))
dtt1$Metric <- factor(dtt1$metric, levels = c("CLR.abd","log_gini","chao1"),
                    labels = c("CLR\nabundance","log Gini-\nSimpson","Richness"))


head(dm1)
head(dtt1)

dtt1 %>% filter(metric=="log_gini") %>% filter(Phylum == "EMF") %>% t()
dtt1 %>% filter(metric=="chao1") %>% filter(Phylum == "EMF") %>% t()




##################### EMF
dm_emf <- dm1 %>% filter(Phylum == "EMF")
dt_emf <- dtt1 %>% filter(Phylum == "EMF")

dt_emf %>% filter(metric=="CLR.abd") %>% t()


p1C <- ggplot(dt_emf) +
  geom_boxplot(data=dm_emf, aes(x = Loc, y = value, fill=Country), outlier.shape = NA) +
  geom_point(data=dm_emf, aes(x = Loc, y = value), color="black",
             pch=21, size=5, position=position_jitter(width=0.1)) +
  scale_fill_manual(values = c("white","grey"))  + 
  
  geom_vline(xintercept = c(2.5, 5.5), linetype="dashed", color="grey") +

  geom_text(data=dt_emf, aes(x = 0.5, y = manual_offset, label = ifelse(MOD != "GLS",
                                                                      paste0("Region: R2=",round(R2_GR1,2)," ",sign_GR1,
                                                                             "\nSite: R2=",round(R2_GR2,2)," ",sign_GR2),
                                                                      paste0("Region: gls ",sign_GR1,
                                                                             "\nSite: gls ",sign_GR2))),
            hjust=0, vjust=1, size=5) +
  
  geom_line(data=data.frame("x"=c(1.5,3.8),"y"=c(7,7), Metric="CLR\nabundance"), aes(x=x, y=y), color="black") +
  geom_line(data=data.frame("x"=c(4.2,7),"y"=c(7,7), Metric="CLR\nabundance"), aes(x=x, y=y), color="black") +
  geom_line(data=data.frame("x"=c(1.5,7),"y"=c(8,8), Metric="CLR\nabundance"), aes(x=x, y=y), color="black") +

  geom_text(data=data.frame("x"=c(2.65),"y"=c(7.3),"label"=c("*"), Metric="CLR\nabundance"), aes(x=x,y=y,label=label),size=8) +
  geom_text(data=data.frame("x"=c(5.8),"y"=c(7.5),"label"=c("ns."), Metric="CLR\nabundance"), aes(x=x,y=y,label=label),size=5) +
  geom_text(data=data.frame("x"=c(4.2),"y"=c(8.3), "label"=c("***"), Metric="CLR\nabundance"), aes(x=x,y=y,label=label),size=8) +  

  ggtitle("EMF") +
  labs(x = "Site") +
  #scale_y_continuous(limits=c(-2.3,-0.3)) +
  facet_grid(Metric~., scales="free_y") +
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
        strip.text.y = element_text(size=20),
        strip.text.x = element_blank()
  )
p1C



png('guilds_08_plot_main_EMF.png', w=1400, h=2800, res=300)
p1C
dev.off()


##################### ##################### ##################### ##################### 
##################### AMF
dm_amf <- dm1 %>% filter(Phylum == "AMF")
dt_amf <- dtt1 %>% filter(Phylum == "AMF")

dt_amf %>% filter(metric=="log_gini") %>% t()
dt_amf %>% filter(metric=="chao1") %>% t()

p1A <- ggplot(dt_amf) +
  geom_boxplot(data=dm_amf, aes(x = Loc, y = value, fill=Country), outlier.shape = NA) +
  geom_point(data=dm_amf, aes(x = Loc, y = value), color="black",
             pch=21, size=5, position=position_jitter(width=0.1)) +
  scale_fill_manual(values = c("white","grey"))  + 
  
  geom_vline(xintercept = c(2.5, 5.5), linetype="dashed", color="grey") +
  
  geom_text(data=dt_amf, aes(x = 0.5, y = manual_offset, label = ifelse(MOD != "GLS",
                                                                        paste0("Region: R2=",round(R2_GR1,2)," ",sign_GR1,
                                                                               "\nSite: R2=",round(R2_GR2,2)," ",sign_GR2),
                                                                        paste0("Region: gls ",sign_GR1,
                                                                               "\nSite: gls ",sign_GR2))),
            hjust=0, vjust=1, size=5) +
  
  geom_line(data=data.frame("x"=c(1.5,3.8),"y"=c(0.2,0.2), Metric="log Gini-\nSimpson"), aes(x=x, y=y), color="black") +
  geom_line(data=data.frame("x"=c(4.2,7),"y"=c(0.2,0.2), Metric="log Gini-\nSimpson"), aes(x=x, y=y), color="black") +
  geom_line(data=data.frame("x"=c(1.5,7),"y"=c(0.55,0.56), Metric="log Gini-\nSimpson"), aes(x=x, y=y), color="black") +
  
  geom_text(data=data.frame("x"=c(2.65),"y"=c(0.3),"label"=c("*"), Metric="log Gini-\nSimpson"), aes(x=x,y=y,label=label),size=8) +
  geom_text(data=data.frame("x"=c(5.8),"y"=c(0.3),"label"=c("***"), Metric="log Gini-\nSimpson"), aes(x=x,y=y,label=label),size=8) +
  geom_text(data=data.frame("x"=c(4.2),"y"=c(0.70), "label"=c("ns."), Metric="log Gini-\nSimpson"), aes(x=x,y=y,label=label),size=5) +  

  geom_line(data=data.frame("x"=c(1.5,3.8),"y"=c(110, 110), Metric="Richness"), aes(x=x, y=y), color="black") +
  geom_line(data=data.frame("x"=c(4.2,7),"y"=c(110, 110), Metric="Richness"), aes(x=x, y=y), color="black") +
  geom_line(data=data.frame("x"=c(1.5,7),"y"=c(130,130), Metric="Richness"), aes(x=x, y=y), color="black") +
  
  geom_text(data=data.frame("x"=c(2.65),"y"=c(112),"label"=c("*"), Metric="Richness"), aes(x=x,y=y,label=label),size=8) +
  geom_text(data=data.frame("x"=c(5.8),"y"=c(112),"label"=c("***"), Metric="Richness"), aes(x=x,y=y,label=label),size=8) +
  geom_text(data=data.frame("x"=c(4.2),"y"=c(138), "label"=c("ns."), Metric="Richness"), aes(x=x,y=y,label=label),size=5) +  
  
  
  
    
  ggtitle("AMF") +
  labs(x = "Site") +
  #scale_y_continuous(limits=c(-2.3,-0.3)) +
  facet_grid(Metric~., scales="free_y") +
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
        strip.text.y = element_text(size=20),
        strip.text.x = element_blank()
  )
p1A



png('guilds_08_plot_main_AMF.png', w=1400, h=2800, res=300)
p1A
dev.off()






##################### ##################### ##################### ##################### 
##################### NF
dm_nf <- dm1 %>% filter(Phylum == "Nitrogen-fixing")
dt_nf <- dtt1 %>% filter(Phylum == "Nitrogen-fixing")

dt_nf %>% filter(metric=="log_gini") %>% t()
dt_nf %>% filter(metric=="chao1") %>% t()

p1B <- ggplot(dt_nf) +
  geom_boxplot(data=dm_nf, aes(x = Loc, y = value, fill=Country), outlier.shape = NA) +
  geom_point(data=dm_nf, aes(x = Loc, y = value), color="black",
             pch=21, size=5, position=position_jitter(width=0.1)) +
  scale_fill_manual(values = c("white","grey"))  + 
  
  geom_vline(xintercept = c(2.5, 5.5), linetype="dashed", color="grey") +
  
  geom_text(data=dt_nf, aes(x = 0.5, y = manual_offset, label = ifelse(MOD != "GLS",
                                                                        paste0("Region: R2=",round(R2_GR1,2)," ",sign_GR1,
                                                                               "\nSite: R2=",round(R2_GR2,2)," ",sign_GR2),
                                                                        paste0("Region: gls ",sign_GR1,
                                                                               "\nSite: gls ",sign_GR2))),
            hjust=0, vjust=1, size=5) +

  geom_line(data=data.frame("x"=c(1.5,3.8),"y"=c(-1, -1), Metric="CLR\nabundance", Grouping="Region"), aes(x=x, y=y), color="black") +
  geom_line(data=data.frame("x"=c(4.2,7),"y"=c(-0.6,-0.6), Metric="CLR\nabundance", Grouping="Region"), aes(x=x, y=y), color="black") +
  geom_line(data=data.frame("x"=c(1.5,7),"y"=c(-0.2,-0.2), Metric="CLR\nabundance", Grouping="Region"), aes(x=x, y=y), color="black") +
  
  geom_text(data=data.frame("x"=c(2.65),"y"=c(-0.9),"label"=c("***"), Metric="CLR\nabundance", Grouping="Region"), aes(x=x,y=y,label=label),size=8) +
  geom_text(data=data.frame("x"=c(5.8),"y"=c(-0.5),"label"=c("***"), Metric="CLR\nabundance", Grouping="Region"), aes(x=x,y=y,label=label),size=8) +
  geom_text(data=data.frame("x"=c(4.2),"y"=c(-0.1), "label"=c("**"), Metric="CLR\nabundance", Grouping="Region"), aes(x=x,y=y,label=label),size=8) +  
  
    
  geom_line(data=data.frame("x"=c(1.5,3.8),"y"=c(0, 0), Metric="log Gini-\nSimpson"), aes(x=x, y=y), color="black") +
  geom_line(data=data.frame("x"=c(4.2,7),"y"=c(0,0), Metric="log Gini-\nSimpson"), aes(x=x, y=y), color="black") +
  geom_line(data=data.frame("x"=c(1.5,7),"y"=c(0.4,0.4), Metric="log Gini-\nSimpson"), aes(x=x, y=y), color="black") +
  
  geom_text(data=data.frame("x"=c(2.65),"y"=c(0.2),"label"=c("ns."), Metric="log Gini-\nSimpson"), aes(x=x,y=y,label=label),size=5) +
  geom_text(data=data.frame("x"=c(5.8),"y"=c(0.1),"label"=c("**"), Metric="log Gini-\nSimpson"), aes(x=x,y=y,label=label),size=8) +
  geom_text(data=data.frame("x"=c(4.2),"y"=c(0.5), "label"=c("***"), Metric="log Gini-\nSimpson"), aes(x=x,y=y,label=label),size=8) +  
  
  geom_line(data=data.frame("x"=c(1.5,3.8),"y"=c(20, 20), Metric="Richness"), aes(x=x, y=y), color="black") +
  geom_line(data=data.frame("x"=c(4.2,7),"y"=c(24,24), Metric="Richness"), aes(x=x, y=y), color="black") +
  geom_line(data=data.frame("x"=c(1.5,7),"y"=c(28,28), Metric="Richness"), aes(x=x, y=y), color="black") +
  
  geom_text(data=data.frame("x"=c(2.65),"y"=c(21),"label"=c("**"), Metric="Richness"), aes(x=x,y=y,label=label),size=8) +
  geom_text(data=data.frame("x"=c(5.8),"y"=c(25),"label"=c("*"), Metric="Richness"), aes(x=x,y=y,label=label),size=8) +
  geom_text(data=data.frame("x"=c(4.2),"y"=c(29), "label"=c("***"), Metric="Richness"), aes(x=x,y=y,label=label),size=8) +  
  
  
  
  
  ggtitle("NF") +
  labs(x = "Site") +
  #scale_y_continuous(limits=c(-2.3,-0.3)) +
  facet_grid(Metric~., scales="free_y") +
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
        strip.text.y = element_text(size=20),
        strip.text.x = element_blank()
  )
p1B



png('guilds_08_plot_main_NF.png', w=1400, h=2800, res=300)
p1B
dev.off()
