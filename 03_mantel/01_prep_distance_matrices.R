rm(list=ls())

library(ggplot2)
library(tidyr)
library(dplyr)
library(dartR)
library(cowplot)
library(RColorBrewer)
library(SNPRelate)
library(plotly)
library(readxl)
library(ggthemes)


brewer.pal(9, "BrBG")

######################
#-------DATA---------#
######################


#vcf.fn <- "../DATA/aspen_genomic_data/microbiome.recode.sort.filter.vcf.recode.vcf"
vcf.fn <- "../../2022_MicrobiomeQuebecMexico/DATA/aspen_genomic_data/microbiome.recode.sort.filter.recode.vcf"
file.size(vcf.fn)
snpgdsVCF2GDS(vcf.fn, "microbiome.recode.sort.filter.gds")
snpgdsSummary("microbiome.recode.sort.filter.gds")


# Opening gds genofile

#snpgdsClose(genofile)
genofile <- snpgdsOpen("microbiome.recode.sort.filter.gds")
genofile
head(genofile)
str(genofile)


### Metadata
dmeta <- read_excel("../../2022_MicrobiomeQuebecMexico/DATA/aspen_genomic_data/GBS.xlsx")
head(dmeta)
dim(dmeta)


########################
#-------colors---------#
########################

col1 <- brewer.pal(3, "BrBG")[1]
col2 <- "black"
col3 <- brewer.pal(3, "BrBG")[3]





########################
#------DISTANCES-------#
########################


####################
#-----GENETIC------#
####################

### Euclidean distance

### Convert GDS to genind
snpset <- read.gdsn(index.gdsn(genofile, "snp.id"))
samples <- read.gdsn(index.gdsn(genofile, "sample.id"))
samples

mat1 <- snpgdsGetGeno(genofile, snp.id=snpset, snpfirstdim=FALSE)
dim(mat1)
mat1[1:10, 1:10]
ncol(mat1)
dgeno <- as.data.frame(mat1)
ncol(dgeno)
colnames(dgeno) <- snpset
rownames(dgeno) <- samples
dgeno[1:10, 1:10]

### Eculidean distance
distgen.eucl <- dist(dgeno, method = "euclidean", diag = TRUE, upper = TRUE, p = 2)
df.eucl <- as.data.frame(as.matrix(distgen.eucl))
class(df.eucl)
head(df.eucl)

write.csv(df.eucl, "05_distances_genetic.csv", col.names=T, row.names=T, quote=FALSE)



####################
#----GEOGRAPHIC----#
####################

head(dmeta)
dgeo <- dmeta %>% dplyr::select(Lat, Long)
rownames(dgeo) <- dmeta$Genotype
dmeta$Genotype
dist.geo <- dist(dgeo, diag=T, upper=T)
df.geo <- as.data.frame(as.matrix(dist.geo))
head(df.geo)
dim(df.geo)

write.csv(df.geo, "05_distances_geo.csv",  col.names=T, row.names=T, quote=FALSE)



####################
#----MICROBIOME----#
####################


### 16S
df.16s <- read.csv("../27_ordination/ordination_02_16S_BCdistance.tsv", sep="\t",row.names=1)
df.16s
head(df.16s)


# vector with new names
load("../01_data_prep/AVS_filtered_phyloseq.RData")
df.sampe <- data.frame(sample_data(pseq.16s)) %>% filter(source=="Soil")
dvec <- df.sampe$Sample.code
names(dvec) <- df.sampe$long.name
dvec


# renaming distance matrix
rownames(df.16s) <- ifelse(
  rownames(df.16s) %in% names(dvec),
  dvec[rownames(df.16s)],
  rownames(df.16s)  # keep unchanged if not in dictionary
)

colnames(df.16s) <- ifelse(
  colnames(df.16s) %in% names(dvec),
  dvec[colnames(df.16s)],
  colnames(df.16s)  # keep unchanged if not in dictionary
)

head(df.16s)


write.csv(df.16s, "05_distances_microbiome_16s.csv",  col.names=T, row.names=T, quote=FALSE)



### ITS
df.its <- read.csv("../27_ordination/ordination_02_ITS_BCdistance.tsv", sep="\t",row.names=1)
head(df.its)


# vector with new names
df.sampe <- data.frame(sample_data(pseq.its)) %>% filter(source=="Soil")
fvec <- df.sampe$Sample.code
names(fvec) <- df.sampe$long.name
fvec

# renaming distance matrix
rownames(df.its) <- ifelse(
  rownames(df.its) %in% names(fvec),
  fvec[rownames(df.its)],
  rownames(df.its)  # keep unchanged if not in dictionary
)

colnames(df.its) <- ifelse(
  colnames(df.its) %in% names(fvec),
  fvec[colnames(df.its)],
  colnames(df.its)  # keep unchanged if not in dictionary
)

head(df.its)




head(df.its)
write.csv(df.its, "05_distances_microbiome_its.csv",  col.names=T, row.names=T, quote=FALSE)







