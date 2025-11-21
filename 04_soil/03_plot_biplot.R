rm(list=ls())

library(tidyverse)
library(cowplot)
library(ggrepel)
library(readxl)
library(RColorBrewer)


tab_blue <- "#1f77b4"
tab_orange <- "#ff7f0e"


#########################
#--------DATA-----------#
#########################


# Data with soil without 2 outliers
data_soil_v2 <- read.csv("soil_01_data_wide.tsv", sep="\t", header=TRUE)
data_soil_v2 %>% head()
df <- data_soil_v2 %>% dplyr::select(C_total,N_total,pH_H2O,P,K,Ca,Mg,Mn,Al,Fe,Na,CEC,C_N)
rownames(df) <- data_soil_v2$short.name
head(df)
dim(df)

#########################
#------- PCA -----------#
#########################

pca.run <- prcomp(df, center=TRUE, scale.=TRUE)
summary(pca.run)


pca.x <- as.data.frame(pca.run$x)
head(pca.x)
class(pca.x)
dim(pca.x)

# Sample with biplot
biplot(prcomp(df, center=TRUE, scale.=TRUE))

write.table(pca.x, file = "soil_04_biplot_PCA.tsv", sep="\t", row.names = T, col.names = T, append=F, quote=F)





### Prep table for plotting

diag(pca.run$sdev)
# Position of arrows
evecs <- pca.run$rotation
evecs
loaded = evecs %*% diag(pca.run$sdev) %>% as.data.frame()
loaded %>% head()
rownames(loaded)
colnames(df)
labels <- c("C","N","pH","P","K","Ca","Mg","Mn","Al","Fe","Na","CEC","C:N")
loaded$lab <- labels
head(loaded)

pca.x <- as.data.frame(pca.run$x)
head(pca.x)
pca.x.size <- merge(pca.x, df, by = 'row.names', sort=FALSE)
head(pca.x.size)
pca.x.size$Country <- data_soil_v2$Country
pca.x.size$Site <- data_soil_v2$Site
colnames(pca.x.size)
head(pca.x.size)


var <- as.data.frame(t(summary(pca.run)$importance))
var$`Proportion of Variance`[1]
var$`Proportion of Variance`
### Basic


#sites_cols <- c("#5E4FA2","#3288BD","#66C2A5","#ABDDA4","#E6F598","#F46D43","#D53E4F","#9E0142")
sites_cols <- c("#134b71","#1f77b4","#258dd5","#86c1ea","#d9ecf9","#ffbe84","#ff7f0e","#c15a00")
site_order <- c("AMOS","STFE","STET","ESSI","FORE","Santiago","FLOR1","FP2")
pca.x.size$Site <- factor(pca.x.size$Site, levels = site_order)

arrow_scale <- 3
r0 <- ggplot(pca.x.size) + aes(x = PC1, y = PC2) +
  #geom_vline(xintercept = 0, color = "grey80", linetype = "dashed") +
  #geom_hline(yintercept = 0, color = "grey80", linetype = "dashed") +
  geom_point(pch=21, aes(fill=Site), size=5) +
  scale_fill_manual(values = sites_cols)  + 
  geom_segment(data = loaded, aes(x = 0, y = 0, xend = V1*arrow_scale*0.75, yend = V2*arrow_scale*0.75),
               arrow = arrow(length = unit(0.3, "cm")), color = "grey30", linewidth=0.5) +
  geom_text_repel(data = loaded, aes(x = V1*arrow_scale, y = V2*arrow_scale, label = lab), 
            color="black", size=7, force=0.75) +
  labs(x = paste0("PC1"," (",round(var$`Proportion of Variance`[1]*100, 1),"%)"),
       y = paste0("PC2"," (",round(var$`Proportion of Variance`[2]*100, 1),"%)")) +
  theme(panel.background = element_rect(fill=NA, color="grey20"),
        panel.grid = element_blank(),
        legend.text = element_text(size=16),
        #legend.position="none",
        legend.title = element_text(size=18),
        legend.key = element_rect(color=NA),
        axis.text = element_text(size=18),
        axis.title = element_text(size=18))
r0



png("soil_04_biplot.png", w=2200, h=1600, res=300)
r0
dev.off()

#pdf("soil_04_biplot.pdf", w=8, h=8)
#r0
#dev.off()

