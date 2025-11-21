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


otu.file <- "../DATA/Aspen_18S_may2024/all/Aspen_18S_may2024.ASV_table_rarefied_19500_dnNA.tsv"
df <- read.csv(otu.file, sep="\t", header=TRUE, comment.char = "")
dim(df)
head(df)

colnames(df)

# Getting otu_table
rownames(df) <- df$X.OTU.ID
head(df)
otu.table <- df %>% dplyr::select(-X.OTU.ID, -taxonomy)
head(otu.table)
OTU = otu_table(otu.table, taxa_are_rows = TRUE)


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
dmeta <- read_excel("../DATA/Aspen_18S_may2024/metadata_18S.xlsx")
head(dmeta)
dmeta.marker <- dmeta %>% filter(marker == '18S') %>% as.data.frame()

# Changing names
colnames(OTU)
dmeta.marker$sampleid
long.names <- paste0("X", lapply(dmeta.marker$sampleid, chartr, old="-", new=".") %>% unlist())
dmeta.marker$long.name <- long.names

rownames(dmeta.marker) <- dmeta.marker$long.name
head(dmeta.marker)

samdata <- sample_data(dmeta.marker)

# Importing to phyloseq
pseq.18s <- phyloseq(OTU, TAX, sam_data = samdata)

pseq.18s.AM <- subset_taxa(pseq.18s, Class == "c__Glomeromycetes")
pseq.18s.AM

############
#--SAVING--#
############


# Save phyloseq object
save(pseq.18s,pseq.18s.AM, file = "AVS_phyloseq_rarefied.RData")




