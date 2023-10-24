#----load packages-------

library(DESeq2)
library(ggplot2)
library(factoextra)
library(FactoMineR)
library(reshape2)
library(ggplot2)
library(Vennerable)
library(rio)
library(ggsci)
library(RColorBrewer)
library(viridis)
library(tidyverse)
library(hrbrthemes)
library(pheatmap)
library(tidyverse)
library(scales)
col <- c('#2a9d8f','#e9c46a','#e63946')
show_col(col)
npg.col <- c("#E7959B","#DB5E67","#CFD99F",'#C5EC41','#E5CD94',"#FC0C00","#7C4374","#339E55","#000376","#2A82C6","#8C6D36","#CB70C0",
             "#EBB854",'#FC8D37',"#63753A","#6D2811","#DD9AD2","#68AADE","#3B397C","#9D9AE5","#B8CF6E","#949494","#BF4F8F","#844346")
options(scipen = 6)
sample.col <- rev(RColorBrewer::brewer.pal(6,"RdBu"))

#----load data----

### load gene info
pig.gene.info <- import('~/project/99.gaodengfeng/01.multi.embory/06.analysis/00.datainfo/01.pig/01.mydata/sus.105.xlsx',which = 'pig_105_clean') # gene annotation info
pig.gene.info <- pig.gene.info[pig.gene.info$biotype == 'protein_coding' | pig.gene.info$biotype == 'lncRNA' | pig.gene.info$biotype == 'pseudogene' | pig.gene.info$biotype == 'processed_pseudogene',]
which(pig.gene.info$genename =='DEFB1')
pig.gene.info[10314,2] <- 'AMELY.X'
pig.gene.info[14647,2] <- 'AMELY.Y'
pig.gene.info[12197,2] <- 'GZMA.1'
pig.gene.info[12198,2] <- 'GZMA.2'
pig.gene.info = pig.gene.info[-c(14831),]

### input raw count
samples <- read.table('../../../01.data/17rnaseq/rna.sample.txt',sep = '\t') 
data <- read.table('../../../04.quant/MSC.ref105.txt',header = T,sep = '\t',comment.char = '#',row.names = 1) # import sample data
data <- data[pig.gene.info[pig.gene.info$biotype=='protein_coding',]$geneid, ]
data <- na.omit(data)
data.clean <- data[,c(6:ncol(data))]
colnames(data.clean) <- samples$V1
data.clean <- na.omit(data.clean)
meta <- data[,c(1:5)]

### calculate TPM 
kb <- meta$Length / 1000
rpk <- data.clean / kb
tpm <- data.frame(t(t(rpk)/colSums(rpk) * 1000000))
rm(kb,rpk)
