rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(cowplot)



#=======================#
#         DATA          #
#=======================#

dm_ <- read.csv("03_alpha_stats.tsv", sep="\t", header=T)
dt1 <- read.csv("03_alpha_16S_country_tests.tsv", sep="\t", header=T)
dt1_ <- read.csv("03_alpha_16S_region_tests.tsv", sep="\t", header=T)

dt2 <- read.csv("03_chao1_16S_country_tests.tsv", sep="\t", header=T)
dt2_ <- read.csv("03_chao1_16S_region_tests.tsv", sep="\t", header=T)

head(dm_)
head(dt)



# Formatting table
head(dm_)
dm_$log_gini <- log(dm_$diversity_gini_simpson + 0.01)
dm_2 <- dm_ %>% filter(source=="Soil" & marker == "16") %>% dplyr::select(long.name, site, Country, region, log_gini, chao1) 
dm <- dm_2 %>% gather(key = "metric", value = "value", c(log_gini, chao1))
head(dm)


head(dt1)
df1 <- dt1 %>% filter(EVAL == "OK")
df1$metric <- "log_gini"
df1$Group <- "site"

df1_ <- dt1_ %>% filter(EVAL == "OK")
df1_$metric <- "log_gini"
df1_$Group <- "region"

df2 <- dt2 %>% filter(EVAL == "OK")
df2$metric <- "chao1"
df2$Group <- "site"

df2_ <- dt2_ %>% filter(EVAL == "OK")
df2_$metric <- "chao1"
df2_$Group <- "region"

dt <- rbind(df1_, df2_)
head(dt)



# Test correction
head(dt)
pvals <- dt %>% filter(statistic %in% c("p_GR1","p_GR2"))  %>% dplyr::select(statistic, metric, Group, V1)
pvals %>% head()
pvals$adj.p <- p.adjust(pvals$V1, method = "BH")
pvals_spread <- pvals %>% dplyr::select(statistic, metric, Group, adj.p) %>% spread(key=statistic,adj.p) %>% dplyr::rename(adj_GR1 = p_GR1, adj_GR2 = p_GR2)
head(pvals_spread)

# Formatting test tables
dt_spread <- dt %>% spread(key=statistic, value=V1)
head(dt_spread)
ds2 <- merge(dt_spread, pvals_spread, by = c("metric","Group"))
ds2

ds2$s1 <- ifelse(ds2$adj_GR1 < 0.05, "*", "ns.")
ds2$s2 <- ifelse(ds2$adj_GR1 < 0.01, "**", ds2$s1)
ds2$sign_GR1 <- ifelse(ds2$adj_GR1 < 0.001, "***", ds2$s2)

ds2$s3 <- ifelse(ds2$adj_GR2 < 0.05, "*", "ns.")
ds2$s4 <- ifelse(ds2$adj_GR2 < 0.01, "**", ds2$s3)
ds2$sign_GR2 <- ifelse(ds2$adj_GR2 < 0.001, "***", ds2$s4)
ds2

head(dm)
dg <- dm %>% dplyr::group_by(metric) %>% summarize(max_val = max(value), min_val = min(value), range=abs(max(value)-min(value))*0.2)
dg
dg$offset <- dg$max_val + dg$range
dg$manual_offset <- c(3100, 0.025)
dg
dtt1 <- merge(ds2, dg, by=c("metric"), sort=F)
head(dtt1)
dtt1








#########################
#--------PLOT-----------#
#########################

dm1 <- dm
head(dm1)
head(dtt1)

#dm11 <- dm1 %>% filter(Group == "site")
#dtt11 <- dtt1 %>% filter(Group == "site")

# c("AMOS","STFE","STET","ESSI","FORE","Santiago","FLOR1","FP2")
#dm$site <- factor(dm$site, levels = c("AB","FORE","FP","SP"), labels = c("AB","FORE","FP2","SP"))
dm1$Loc <- factor(dm1$site, levels = c("AMOS","STFE","STET","ESSI","FORE","Santiago","FLOR1","FP2"),
                  labels <- c("AMOS","STFE","STET","ESSI","FORE","Santiago","FLOR1","FP2"))
#dm11$Grouping <- factor(dm11$Group, levels = c("site","region"),
#                       labels = c("Population","Region"))
dm1$Metric <- factor(dm1$metric, levels = c("log_gini","chao1"),
                     labels = c("log Gini-\nSimpson","Richness"))


#dtt1$Grouping <- factor(dtt1$Group, levels = c("site","region"),
#                        labels = c("Population","Region"))
dtt1$Metric <- factor(dtt1$metric, levels = c("log_gini","chao1"),
                      labels = c("log Gini-\nSimpson","Richness"))


head(dm1)
head(dtt1)

dtt1 %>% filter(metric=="log_gini") %>% t()
dtt1 %>% filter(metric=="chao1") %>% t()

p1C <- ggplot(dtt1) +
  geom_boxplot(data=dm1, aes(x = Loc, y = value, fill=Country), outlier.shape = NA) +
  geom_point(data=dm1, aes(x = Loc, y = value), color="black",
             pch=21, size=5, position=position_jitter(width=0.1)) +
  scale_fill_manual(values = c("white","grey"))  + 
  
  geom_vline(xintercept = c(2.5, 5.5), linetype="dashed", color="grey") +
  
  
  geom_text(data=dtt1, aes(x = 0.5, y = manual_offset, label = ifelse(MOD %in% c("ANOVA","ANOVA_transformed"),
                                                                      paste0("Region: R2=",round(R2_GR1,2)," ",sign_GR1,
                                                                             "\nSite: R2=",round(R2_GR2,2)," ",sign_GR2),
                                                                      paste0("Region: gls ",sign_GR1,
                                                                             "\nSite: gls ",sign_GR2))),
            hjust=0, vjust=1, size=5) + 
  
  geom_line(data=data.frame("x"=c(1.5,3.8),"y"=c(0.01,0.01), Metric="log Gini-\nSimpson", Grouping="site"), aes(x=x, y=y), color="black") +
  geom_line(data=data.frame("x"=c(4.2,7),"y"=c(0.01,0.01), Metric="log Gini-\nSimpson", Grouping="site"), aes(x=x, y=y), color="black") +
  geom_line(data=data.frame("x"=c(1.5,7),"y"=c(0.013,0.013), Metric="log Gini-\nSimpson", Grouping="site"), aes(x=x, y=y), color="black") +
  
  geom_text(data=data.frame("x"=c(2.65),"y"=c(0.012),"label"=c("ns."), Metric="log Gini-\nSimpson", Grouping="site"), aes(x=x,y=y,label=label),size=5) +
  geom_text(data=data.frame("x"=c(5.8),"y"=c(0.012),"label"=c("ns."), Metric="log Gini-\nSimpson", Grouping="site"), aes(x=x,y=y,label=label),size=5) +
  geom_text(data=data.frame("x"=c(4.2),"y"=c(0.014), "label"=c("***"), Metric="log Gini-\nSimpson", Grouping="site"), aes(x=x,y=y,label=label),size=8) +  
  
  
  labs(x = "Site") +
  #scale_y_continuous(limits=c(-2.3,-0.3)) +
  facet_grid(Metric~., scales="free") +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(size=20),
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




  
png('07_plot_main_16S.png', w=1200, h=1800, res=300)
p1C
dev.off()

