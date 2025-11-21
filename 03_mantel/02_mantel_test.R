rm(list=ls())

library(ggplot2)
library(tidyr)
library(dplyr)
library(dartR)
library(cowplot)
library(RColorBrewer)
library(plotly)
library(readxl)
library(ggthemes)
library(vegan)



######################
#-------DATA---------#
######################


### Metadata
dmeta <- read_excel("../../2022_MicrobiomeQuebecMexico/DATA/aspen_genomic_data/GBS.xlsx")
head(dmeta)
dim(dmeta)

df.geo <- read.csv("05_distances_geo.csv", sep=",", header=T, row.names=1)
head(df.geo)
dim(df.geo)

df.genet <- read.csv("05_distances_genetic.csv", sep=",", header=T, row.names=1)
head(df.genet)
dim(df.genet)

df.bac <- read.csv("05_distances_microbiome_16s.csv", sep=",", header=T, row.names=1)
head(df.bac)
dim(df.bac)

df.fun <- read.csv("05_distances_microbiome_its.csv", sep=",", header=T, row.names=1)
head(df.fun)
dim(df.fun)

rownames(df.geo) == rownames(df.genet)
rownames(df.geo) == rownames(df.bac)


# distance matrices need to be sorted according to row and column names




########################
#-------colors---------#
########################

col1 <- brewer.pal(3, "BrBG")[1]
col2 <- "black"
col3 <- brewer.pal(3, "BrBG")[3]





##############################################

#               Mantel all

##############################################

# 16S Canada
# 16S Mexico

# ITS Canada
# ITS Mexico


canada <- c("ESSI","AMOS","STFE","STET","FORE")
mexico <- c("SAAA","MEXA","MEXP")



### Canada

### 16S
### Samples common to genetic and microbial datasets
head(df.genet)
tab.genet <- data.frame("rowname" = intersect(rownames(df.genet), rownames(df.bac)),
                        "colname" = intersect(colnames(df.genet), colnames(df.bac)))
tab.genet$site <- sub("\\..*", "", tab.genet$colname)

# Samples unique to Canada
df.samples <- tab.genet %>% filter(site %in% canada)
samples.cols <- df.samples$colname
samples.rows <- df.samples$rowname

# Subsetting matrices
can.genet <- df.genet[samples.rows, samples.cols]
can.bac <- df.bac[samples.rows, samples.cols]
dist.can.genet <- as.dist(can.genet, upper=F, diag=F)
dist.can.bac <- as.dist(can.bac, upper=F, diag=F)



### ITS
### Samples common to genetic and microbial datasets
tab.genet <- data.frame("rowname" = intersect(rownames(df.genet), rownames(df.fun)),
                        "colname" = intersect(colnames(df.genet), colnames(df.fun)))
tab.genet$site <- sub("\\..*", "", tab.genet$colname)

# Samples unique to Canada
df.samples <- tab.genet %>% filter(site %in% canada)
samples.cols <- df.samples$colname
samples.rows <- df.samples$rowname

# Subsetting matrices
can.genet <- df.genet[samples.rows, samples.cols]
can.fun <- df.fun[samples.rows, samples.cols]
disti.can.genet <- as.dist(can.genet, upper=F, diag=F)
disti.can.fun <- as.dist(can.fun, upper=F, diag=F)


### Mantel test - Canada


mantel.bac.can <- mantel(dist.can.genet, dist.can.bac, permutations = 10000, na.rm=T, method="spearman")
mantel.fun.can <- mantel(disti.can.genet, disti.can.fun, permutations = 10000, na.rm=T, method="spearman")

mantel.results.can <- data.frame("test" = c("Bacteria","Fungi"),
                             "region" = c("Canada","Canada"),
                             "stat" = c(mantel.bac.can$statistic,
                                        mantel.fun.can$statistic),
                             "pvalue" = c(mantel.bac.can$signif,
                                          mantel.fun.can$signif))
mantel.results.can








### Mexico

### 16S
### Samples common to genetic and microbial datasets
head(df.genet)
tab.genet <- data.frame("rowname" = intersect(rownames(df.genet), rownames(df.bac)),
                        "colname" = intersect(colnames(df.genet), colnames(df.bac)))
tab.genet$site <- sub("\\..*", "", tab.genet$colname)

# Samples unique to Mexico
df.samples <- tab.genet %>% filter(site %in% mexico)
samples.cols <- df.samples$colname
samples.rows <- df.samples$rowname

# Subsetting matrices
can.genet <- df.genet[samples.rows, samples.cols]
can.bac <- df.bac[samples.rows, samples.cols]
dist.can.genet <- as.dist(can.genet, upper=F, diag=F)
dist.can.bac <- as.dist(can.bac, upper=F, diag=F)



### ITS
### Samples common to genetic and microbial datasets
tab.genet <- data.frame("rowname" = intersect(rownames(df.genet), rownames(df.fun)),
                        "colname" = intersect(colnames(df.genet), colnames(df.fun)))
tab.genet$site <- sub("\\..*", "", tab.genet$colname)

# Samples unique to Mexico
df.samples <- tab.genet %>% filter(site %in% mexico)
samples.cols <- df.samples$colname
samples.rows <- df.samples$rowname

# Subsetting matrices
can.genet <- df.genet[samples.rows, samples.cols]
can.fun <- df.fun[samples.rows, samples.cols]
disti.can.genet <- as.dist(can.genet, upper=F, diag=F)
disti.can.fun <- as.dist(can.fun, upper=F, diag=F)


### Mantel test - Mexico


mantel.bac.mx <- mantel(dist.can.genet, dist.can.bac, permutations = 10000, na.rm=T, method="spearman")
mantel.fun.mx <- mantel(disti.can.genet, disti.can.fun, permutations = 10000, na.rm=T, method="spearman")

mantel.results.mx <- data.frame("test" = c("Bacteria","Fungi"),
                                 "region" = c("Mexico","Mexico"),
                                 "stat" = c(mantel.bac.mx$statistic,
                                            mantel.fun.mx$statistic),
                                 "pvalue" = c(mantel.bac.mx$signif,
                                              mantel.fun.mx$signif))
mantel.results.mx



### Combine Mantel

mantel.res <- rbind(mantel.results.can, mantel.results.mx)
mantel.res$p.adj <- p.adjust(mantel.res$pvalue, "BH")


write.table(mantel.res, file = "06_compare_distances_mantel.tsv", sep="\t", row.names = F, col.names = T, quote=F, append=F)




