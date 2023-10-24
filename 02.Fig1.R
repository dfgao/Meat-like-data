library(ComplexHeatmap)
orth <- import('../../pigref105_human.unique.xlsx',which = 'Sheet1')

#----load data-------

## PEF 
pef <- read.table('../../../04.quant/PEF.ref105.txt',header = T,sep = '\t',comment.char = '#',row.names = 1)
pef <- pef[pig.gene.info[pig.gene.info$biotype=='protein_coding',]$geneid, ]
pef <- na.omit(pef)
pef.clean <- pef[,c(6:ncol(pef))]
colnames(pef.clean) <- c('PEF.1','PEF.2')
pef.clean <- na.omit(pef.clean)

## calculate TPM 
kb <- meta$Length / 1000
rpk <- pef.clean / kb
pef.tpm <- data.frame(t(t(rpk)/colSums(rpk) * 1000000))
rm(kb,rpk)

## get data for analysis
pg.pef.count <- cbind(data.clean[,c(10:12,7:9)],pef.clean)
pg.pef.tpm <- cbind(tpm[,c(10:12,7:9)],pef.tpm)

# Fig.1e ------
pg.pef.count <- pg.pef.count %>% rownames_to_column(var = 'geneid') %>% left_join( y = pig.gene.info[,c(1,2)], by = 'geneid')
rownames(pg.pef.count) <- pg.pef.count$genename
pg.pef.tpm <- pg.pef.tpm %>% rownames_to_column(var = 'geneid') %>% left_join( y = pig.gene.info[,c(1,2)], by = 'geneid')
rownames(pg.pef.tpm) <- pg.pef.tpm$genename

plu.gene <- c('NANOG','SOX2','ENSSSCG00000001393','OTX2','LIN28A','TCF3','FGF2','LEFTY2','SMARCAD1','MYST3','SETDB1','JARID2','RIF1','REST')
tpm.plu <- na.omit(pg.pef.tpm[plu.gene,-c(1,ncol(pg.pef.tpm))])

ann <- data.frame(Cell_type = c(rep('pgEpiSCs.Low',3), rep('pgEpiSCs.High',3), rep('pEF',2)))
rownames(ann) <- colnames(tpm.plu)

ha_top.col <- list(Cell_type = c(sample.col[1:2],'#e9c46a'))
names(ha_top.col$Cell_type) <- factor(c('pgEpiSCs.Low', 'pgEpiSCs.High','pEF'))

pheatmap(t(scale(t(tpm.plu))),
         border_color =  NA, 
         cluster_rows = T,
         cluster_cols = F,
         annotation_col = ann,
         annotation_colors = ha_top.col,
         color = scales::alpha(rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"),alpha=T,bias=.8)(256)),alpha = 1),
         angle_col = '315',
         breaks = unique(c(seq(-2,2, length=256)))
         )
# Fig.1f --------


