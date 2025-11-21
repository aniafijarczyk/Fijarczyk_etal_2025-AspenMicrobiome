rm(list=ls())


library(openxlsx)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(phyloseq)
library(microbiome)
library(cowplot)



#=======================#
#--------DATA-----------#
#=======================#


dat_soil <- read.xlsx(xlsxFile = "../../2022_MicrobiomeQuebecMexico/DATA/results/soil.xlsx")
head(dat_soil)

# Climate
dclim <- read.csv("../13_climate/02_format_table.tsv", sep="\t", header=T)
head(dclim)

df.soil <- merge(dat_soil, dclim, by.x = "Site", by.y="SampleID", sort=F)
head(df.soil)

### Sampling year
unique(df.soil$Site)
df.year <- data.frame("Site" = c("ESSI","AMOS","STFE","FORE","STET","Santiago","FLOR1","FP2" ),
                      "Year" = c("2018","2018","2018","2018","2019","2019","2019","2019"))
head(df.year)

data_soil <- merge(df.soil, df.year, by = "Site", sort=F)
head(data_soil)

# Converting numeric columns to numbers
data_soil[8:23] %>% head()
data_soil[8:23] <- data_soil[8:23] %>% mutate_if(is.character, as.numeric)
data_soil$C_N <- data_soil$C_total/data_soil$N_total
data_soil$Year <- as.numeric(data_soil$Year)
data_soil %>% head()


# Checking for outliers
hist(data_soil$N_total, breaks = 100)
hist(data_soil$C_total, breaks = 100)
hist(data_soil$C_N, breaks = 100)
hist(data_soil$P, breaks = 100)
hist(data_soil$MAT, breaks = 100)
hist(data_soil$PPT_WS, breaks = 100)

dim(data_soil)
head(data_soil)




########################
#  Combining datasets  #
########################

# Data with diversity estimates
meta.alpha.1 <- read.csv("03_alpha_stats.tsv", sep='\t', header=TRUE)
dim(meta.alpha.1)
meta.alpha.1$marker
meta.alpha.1 <- meta.alpha.1 %>% dplyr::filter(marker == "16") %>% dplyr::select(source, site, Country, region, short.code, diversity_gini_simpson, diversity_shannon)
meta.alpha.1 %>% head()

#meta.alpha.2 <- read.csv("04_alpha_ITS_stats.tsv", sep='\t', header=TRUE)
meta.alpha.2 <- meta.alpha.1 %>% dplyr::filter(marker == "ITS") %>% dplyr::select(source, site, Country, region, short.code, diversity_gini_simpson, diversity_shannon)
meta.alpha.2$marker
meta.alpha.2 <- meta.alpha.2 %>% dplyr::select(source, short.code, diversity_gini_simpson, diversity_shannon)
meta.alpha.2 <- meta.alpha.2 %>% dplyr::rename("diversity_gini_simpson_its" = "diversity_gini_simpson",
                               "diversity_shannon_its" = "diversity_shannon")
meta.alpha.2 %>% head()

meta.alpha <- merge(meta.alpha.1, meta.alpha.2, by = c("source","short.code"))
head(meta.alpha)
dim(meta.alpha)

# Orders of metrics
site_order <- c("AMOS","STFE","STET","ESSI","FORE","Santiago","FLOR1","FP2")

meta.alpha.soil <- meta.alpha %>% filter(source == 'Soil')
meta.alpha.soil %>% head()
meta.alpha.soil %>% dim()

df_soil <- data_soil %>% dplyr::select(short.code, C_N, C_total, N_total, P, pH_H2O, 
                                       K, Mg, Ca, Mn, Al, Fe, Na, CEC, Sable, Limon, Argile, Classe.texturale,
                                       MAT, PPT_WS, TD, Year)
head(df_soil)
dim(df_soil)

#df <- merge(meta.alpha.soil, df_soil, by = 'short.code', sort=FALSE, all.x = T)
df <- merge(meta.alpha.soil, df_soil, by = 'short.code', sort=FALSE)
df <- df %>% dplyr::rename("Site"="site")
dim(df)
dim(df_soil)
df %>% dplyr::group_by(Country, Site) %>% dplyr::summarise(n = n())
dim(df)
df$Site <- as.factor(df$Site)
df$Country <- as.factor(df$Country)
df$num <- c(1:length(df$short.code))
df %>% head()
dim(df)


df %>% arrange(diversity_gini_simpson) %>% head()

write.table(df, file = "15_diversity_vs_soil_tables.tsv", sep="\t", col.names = T, row.names = F, append=F, quote=F)

df %>% head()



