rm(list = ls())

library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(phyloseq)
library(microbiome)
library(readxl)
library(stringr)
library(RColorBrewer)
library(plyr)
library(openxlsx)
library(ggrepel)
library(vegan)
library(ggnewscale)


source("run_NESTED_ANOVA_EXT_3LEV.R")



#########################
#--------DATA-----------#
#########################


########
### SOIL

data_soil <- read.xlsx("../../2022_MicrobiomeQuebecMexico/DATA/results/soil.xlsx")
data_soil %>% head()
data_soil$C_N <- as.numeric(data_soil$C_total)/as.numeric(data_soil$N_total)
metrics <- c("C_total","N_total","pH_H2O","P","K","Ca","Mg","Mn","Al","Fe","Na","CEC","C_N")
data_soil <- data_soil %>% mutate_at(metrics, as.numeric)
dim(data_soil)
data_soil %>% head()

data_soil$region1 <- ifelse(data_soil$Site %in% c("AMOS","STFE"),"Boreal","Rest")
data_soil$region2 <- ifelse(data_soil$Site %in% c("STET","ESSI","FORE"),"Cold_temperate",data_soil$region1)
data_soil$region <- ifelse(data_soil$Site %in% c("Santiago","FLOR1","FP2"),"Warm_temperate",data_soil$region2)
data_soil$region <- factor(data_soil$region, levels=c("Boreal","Cold_temperate","Warm_temperate"))
data_soil %>% head()

write.table(data_soil, file = "soil_01_data_wide.tsv", sep="\t", row.names = F, col.names = T, append=F, quote=F)


dg_soil <- data_soil %>% gather(key = "metric", value="value", all_of(metrics), -short.name, -Site,-Country, -region)
dg_soil <- dg_soil %>% dplyr::select(short.name, Site, Country, region, metric, value)
head(dg_soil)

write.table(dg_soil, file = "soil_01_data.tsv", sep="\t", row.names = F, col.names = T, append=F, quote=F)






### Colors
source_cols <- c("#0072B2","#D55E00","#F0E442")
country_cols <- c("#E84A5F","#DDCC77")
site_cols <- brewer.pal(9,"BrBG")
phylum_cols <- c("#332288","#117733","#44AA99","#88CCEE","#DDCC77","#CC6677","#AA4499","#882255")
phylum_cols_wong <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

tab_blue <- "#1f77b4"
tab_orange <- "#ff7f0e"
site_order <- c("AMOS","STFE","STET","ESSI","FORE","Santiago","FLOR1","FP2")







#####################
#-----  TESTS ------#
#####################



###############
### TESTS of difference between countries

head(dg_soil)


outputs <- list()
for (fam in unique(dg_soil$metric)) {
  print(fam)
  df_ <- dg_soil[dg_soil$metric == fam,]
  df.models <- run_nested_models(df_, "value", "Country", "region", "Site")
  df.models$metric <- fam
  outputs[[fam]] <- df.models
}

length(outputs)
outputs[[1]]
dOUT <- bind_rows(outputs)
dOUT %>% filter(MOD == "GLS" & EVAL == "OK")
dOUT %>% filter(metric == "P")



write.table(dOUT, file = "soil_01_data_country_tests.tsv", sep="\t", row.names = F, col.names = T, append=F, quote=F)








##########################################################################################
### TESTS of difference between REGIONS

head(dg_soil)


outputs <- list()
for (fam in unique(dg_soil$metric)) {
  print(fam)
  df_ <- dg_soil[dg_soil$metric == fam,]
  df.models <- run_nested_models(df_, "value", "region", "Site")
  df.models$metric <- fam
  outputs[[fam]] <- df.models
}

length(outputs)
outputs[[1]]
dOUT <- bind_rows(outputs)
dOUT %>% filter(MOD == "GLS" & EVAL == "OK")


write.table(dOUT, file = "soil_01_data_region_tests.tsv", sep="\t", row.names = F, col.names = T, append=F, quote=F)






###############
### TESTS of difference between countries - 3 levels

head(dg_soil)
df_ <- dg_soil[dg_soil$metric == "C_total",]
mod <- lm(value ~ Country / region / Site, df_)
anova(mod)
df.models <- run_nested_models(df_, "value", "Country", "region", "Site")
df.models




# Run for all metrics

outputs <- list()
for (fam in unique(dg_soil$metric)) {
  print(fam)
  df_ <- dg_soil[dg_soil$metric == fam,]
  df.models <- run_nested_models(df_, "value", "Country", "region","Site")
  df.models$metric <- fam
  outputs[[fam]] <- df.models
}

length(outputs)
outputs[[1]]
dOUT <- bind_rows(outputs)
dOUT %>% filter(MOD == "GLS" & EVAL == "OK")
dOUT %>% filter(metric == "P")



write.table(dOUT, file = "soil_01_data_country_tests.tsv", sep="\t", row.names = F, col.names = T, append=F, quote=F)







