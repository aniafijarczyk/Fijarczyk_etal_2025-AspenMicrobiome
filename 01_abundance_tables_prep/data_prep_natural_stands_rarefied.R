rm(list = ls())
setwd("D:/NRCan/2024_Aspen/01_data_prep")


library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(phyloseq)
library(microbiome)
library(readxl)
library(stringr)

####################################
#---- Installing packages----------#
####################################

#if(!requireNamespace("BiocManager")){
#  install.packages("BiocManager")
#}
#BiocManager::install("phyloseq")
#packageVersion('phyloseq')
# [1] ‘1.40.0’

### Installing microbiome

#library(BiocManager)
#BiocManager::install("microbiome")


#########################
#--------DATA-----------#
#########################
# Following this tutorial: http://joey711.github.io/phyloseq/import-data


#########
#--16S--#
#########


otu.file.16s <- "../../2022_MicrobiomeQuebecMexico/DATA/16S/Full/Peuplier_REPRISE_16S_oct2022.ASV_table_rarefied_6370_dnNA.tsv"
df <- read.csv(otu.file.16s, sep="\t", header=TRUE, comment.char = "")
head(df)
colnames(df)

# Removing poor samples
#poor <- read.csv("../01_data_explore/dataExplore_03_abundances_16s_poor_samples_below_3000_above_190000.tsv", sep="\t", header=TRUE)
#samples.to.exclude <- as.vector(poor$long.name)
#df <- dfr %>% dplyr::select(-all_of(samples.to.exclude))
#dim(df)

# Getting otu_table
rownames(df) <- df$X.OTU.ID
head(df)
otu.table <- df %>% dplyr::select(-X.OTU.ID, -taxonomy)
head(otu.table)
OTU = otu_table(otu.table, taxa_are_rows = TRUE)
OTU

# Getting tax_table

getTaxTab <- function(dataframe) {
  
  tax.table <- dataframe %>% dplyr::select(taxonomy)
  head(tax.table)
  
  row_domain <- c()
  row_phylum <- c()
  row_class <- c()
  row_order <- c()
  row_family <- c()
  row_genus <- c()
  row_species <- c()
  
  for (i in 1:length(tax.table$taxonomy)) {
    row = tax.table$taxonomy[i]
    
    # Getting d from row
    unlist_row <- unlist(strsplit(row, ";"))
    d_ele <- NA
    p_ele <- NA
    c_ele <- NA
    o_ele <- NA
    f_ele <- NA
    g_ele <- NA
    s_ele <- NA  
    
    for (ele in unlist_row) {
      if (str_detect(ele, "d__")) { d_ele <- gsub(" ","", ele, fixed=TRUE) } 
      if (str_detect(ele, "p__")) { p_ele <- gsub(" ","", ele, fixed=TRUE) } 
      if (str_detect(ele, "c__")) { c_ele <- gsub(" ","", ele, fixed=TRUE) } 
      if (str_detect(ele, "o__")) { o_ele <- gsub(" ","", ele, fixed=TRUE) } 
      if (str_detect(ele, "f__")) { f_ele <- gsub(" ","", ele, fixed=TRUE) } 
      if (str_detect(ele, "g__")) { g_ele <- gsub(" ","", ele, fixed=TRUE) } 
      if (str_detect(ele, "s__")) { s_ele <- gsub(" ","", ele, fixed=TRUE) } 
    }  
    
    if (!is.na(d_ele)) { row_domain <- c(row_domain, d_ele)
    } else { row_domain <- c(row_domain, NA) }
    
    if (!is.na(p_ele)) { row_phylum <- c(row_phylum, p_ele)
    } else { row_phylum <- c(row_phylum, NA) }
   
    if (!is.na(c_ele)) { row_class <- c(row_class, c_ele)
    } else { row_class <- c(row_class, NA) }
    
    if (!is.na(o_ele)) { row_order <- c(row_order, o_ele)
    } else { row_order <- c(row_order, NA) }
    
    if (!is.na(f_ele)) { row_family <- c(row_family, f_ele)
    } else { row_family <- c(row_family, NA) }
    
    if (!is.na(g_ele)) { row_genus <- c(row_genus, g_ele)
    } else { row_genus <- c(row_genus, NA) }
    
    if (!is.na(s_ele)) { row_species <- c(row_species, s_ele)
    } else { row_species <- c(row_species, NA) }  
  }
  
  tax.table$Domain <- row_domain
  tax.table$Phylum <- row_phylum
  tax.table$Class <- row_class
  tax.table$Order <- row_order
  tax.table$Family <- row_family
  tax.table$Genus <- row_genus
  tax.table$Species <- row_species
  
  tax.matrix.spread <- tax.table %>% dplyr::select(-taxonomy) %>% as.matrix()
  #head(tax.matrix.spread)
  return(tax.matrix.spread)
  
}

tax.matrix.spread <- getTaxTab(df)

TAX = tax_table(tax.matrix.spread)

# Getting sample_data
dmeta <- read_excel("../../2022_MicrobiomeQuebecMexico/METADATA/dataExplore_03_abundances_combined_metadata.xlsx")
dmeta.16s <- dmeta %>% filter(marker == '16') %>% as.data.frame()
colnames(dmeta.16s)[5] <- "long.name"
rownames(dmeta.16s) <- dmeta.16s$long.name
head(dmeta.16s)
dmeta.16s.flt <- dmeta.16s %>% filter(long.name %in% colnames(OTU))
dim(dmeta.16s.flt)

# Sorting samples
length(colnames(OTU))
dmeta.16s.flt$long.name
ds <- data.frame("long.name"=colnames(OTU))
head(ds)
dm <- merge(ds, dmeta.16s.flt, by = "long.name", sort=F)
head(dm)
dim(dm)
dm$long.name
colnames(OTU)
rownames(dm) <- dm$long.name

samdata <- sample_data(dm)



# Importing to phyloseq
rseq.16s <- phyloseq(OTU, TAX, sam_data = samdata)
rseq.16s




#########
#--ITS--#
#########


otu.file.its <- "../../2022_MicrobiomeQuebecMexico/DATA/ITS/Full/Peuplier_REPRISE_ITS_oct2022.ASV_table_rarefied_11982_dnNA.tsv"
df <- read.csv(otu.file.its, sep="\t", header=TRUE, comment.char = "")
head(df)
dim(df)

colnames(df)
write.table(data.frame("cols" = colnames(df)), file="check_columns_2.txt", sep="\t", quote=F, append=F)

# Getting otu_table
rownames(df) <- df$X.OTU.ID
head(df)
otu.table <- df %>% dplyr::select(-X.OTU.ID, -taxonomy)
head(otu.table)
dim(otu.table)
OTU = otu_table(otu.table, taxa_are_rows = TRUE)
dim(OTU)


# Getting tax_table
tax.matrix.spread <- getTaxTab(df)
TAX = tax_table(tax.matrix.spread)


# Getting sample_data
dmeta <- read_excel("../../2022_MicrobiomeQuebecMexico/METADATA/dataExplore_03_abundances_combined_metadata.xlsx")
dmeta.its <- dmeta %>% filter(marker == 'ITS') %>% as.data.frame()
head(dmeta.its)
colnames(dmeta.its)[5] <- "long.name"
rownames(dmeta.its) <- dmeta.its$long.name
head(dmeta.its)
dim(dmeta.its)
colnames(OTU)
dmeta.its.flt <- dmeta.its %>% filter(long.name %in% colnames(OTU))
dim(dmeta.its.flt)

# Sorting samples
length(colnames(OTU))
dmeta.its.flt$long.name
ds <- data.frame("long.name"=colnames(OTU))
head(ds)
dm <- merge(ds, dmeta.its.flt, by = "long.name", sort=F)
head(dm)
dim(dm)
dm$long.name
colnames(OTU)
rownames(dm) <- dm$long.name
dim(dm)

samdata <- sample_data(dm)



# Importing to phyloseq
rseq.its <- phyloseq(OTU, TAX, sam_data = samdata)


############
#--SAVING--#
############


# Save phyloseq object
save(rseq.16s, rseq.its, file = "AVS_filtered_phyloseq_rarefied.RData")







### 


head(dm)
dm$long.name
dm %>% group_by(site) %>% dplyr::summarize(n = n())









plot_bar(pseq.16s, fill = "Domain")

### CLR transform

pseq.16s.clr <- microbiome::transform(pseq.16s, 'clr')
pseq.16s.clr@otu_table

