rm(list = ls())

library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(phyloseq)
library(microbiome)
library(readxl)
library(stringr)


#########################
#--------DATA-----------#
#########################

#########
#--16S--#
#########

otu.file.16s <- "../../2022_MicrobiomeQuebecMexico/DATA/16S/Full/Peuplier_REPRISE_16S_oct2022.ASV_table_norarefaction_dnNA.tsv"
df <- read.csv(otu.file.16s, sep="\t", header=TRUE, comment.char = "")
head(df)


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
tax.matrix.spread

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
pseq.16s <- phyloseq(OTU, TAX, sam_data = samdata)
pseq.16s




#########
#--ITS--#
#########


otu.file.its <- "../../2022_MicrobiomeQuebecMexico/DATA/ITS/Full/Peuplier_REPRISE_ITS_oct2022.ASV_table_norarefaction_dnNA.tsv"
df <- read.csv(otu.file.its, sep="\t", header=TRUE, comment.char = "")
head(df)
dim(df)
df[1:10, 1:10]


# Removing poor samples
# Getting otu_table
rownames(df) <- df$X.OTU.ID
head(df)
otu.table <- df %>% dplyr::select(-X.OTU.ID, -taxonomy, -X139.control.ITS)
head(otu.table)
df_check <- as.data.frame(apply(otu.table, 2, sum))
colnames(df_check) <- c("reads")
df_check$ID <- rownames(df_check)
df_check %>% arrange(reads) %>% head()

poor_samples <- df_check %>% filter(reads < 1000) %>% pull(ID)
poor_samples

# Getting otu_table
rownames(df) <- df$X.OTU.ID
head(df)
otu.table <- df %>% dplyr::select(-X.OTU.ID, -taxonomy, -X139.control.ITS) %>% dplyr::select(!any_of(poor_samples))
head(otu.table)
#write.table(data.frame("cols" = colnames(otu.table)), file="check_columns_1.txt", sep="\t", quote=F, append=F)

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
pseq.its <- phyloseq(OTU, TAX, sam_data = samdata)


############
#--SAVING--#
############


# Save phyloseq object
save(pseq.16s, pseq.its, file = "AVS_filtered_phyloseq.RData")



### Save as tables
dasv <- as.data.frame(otu_table(pseq.16s))
dasv[1:10, 1:10]
write.table(dasv, "dataPrep_01_overview_16S_ASV.tsv", sep="\t", col.names = T,row.names = T,append=F, quote=F)

dtax <- as.data.frame(tax_table(pseq.16s))
class(dtax)
write.table(dtax, "dataPrep_01_overview_16S_TAX.tsv", sep="\t", col.names = T,row.names = T,append=F, quote=F)

dmeta <- data.frame(sample_data(pseq.16s))
class(dmeta)
head(dmeta)
write.table(dmeta, "dataPrep_01_overview_16S_META.tsv", sep="\t", col.names = T,row.names = T,append=F, quote=F)


### ITS
dasv <- as.data.frame(otu_table(pseq.its))
dasv[1:10, 1:10]
write.table(dasv, "dataPrep_01_overview_ITS_ASV.tsv", sep="\t", col.names = T,row.names = T,append=F, quote=F)

dtax <- as.data.frame(tax_table(pseq.its))
class(dtax)
write.table(dtax, "dataPrep_01_overview_ITS_TAX.tsv", sep="\t", col.names = T,row.names = T,append=F, quote=F)

dmeta <- data.frame(sample_data(pseq.its))
class(dmeta)
head(dmeta)
write.table(dmeta, "dataPrep_01_overview_ITS_META.tsv", sep="\t", col.names = T,row.names = T,append=F, quote=F)





