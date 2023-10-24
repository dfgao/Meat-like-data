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
# Fig.1d --------
## get tpm genes > 0.5 in each group 
pick.genes <- unique(c(rownames(pg.pef.tpm[(matrixStats::rowMedians(as.matrix(pg.pef.tpm[,c(2:9)]),cols = 1:3) > 0.5) ,]), rownames(pg.pef.tpm[(matrixStats::rowMedians(as.matrix(pg.pef.tpm[,c(2:9)]),cols = 4:6) > 0.5) ,]), rownames(pg.pef.tpm[(matrixStats::rowMedians(as.matrix(pg.pef.tpm[,c(2:9)]),cols = 7:8) > 0.5) ,])))
group_list_pg <- factor(c( rep('pgEpiSCs.Low',3), rep('pgEpiSCs.High',3),rep('PEF',2)),levels = c('pgEpiSCs.Low','pgEpiSCs.High','PEF'))

d_p <- function(tpm,group_list){
  tpm.t = as.data.frame(t(tpm))
  tpm.t = cbind(tpm.t,group_list)
  tpm.pca <- PCA(tpm.t[,-ncol(tpm.t)],graph = FALSE)
  fviz_screeplot(tpm.pca, addlabels = TRUE, ylim = c(0, 70),title='Dim choose')
  fviz_pca_ind(tpm.pca,
               mean.point=F,
               # repel = T,
               label = "none",
               geom.ind = c("point",'text'),
               fill.ind = tpm.t$group_list,
               palette = c(sample.col[1:2],'#e9c46a'),
               legend.title = "Groups",
               pointsize = 6,
               pointshape = 21,
               col.ind = "black",
               title = 'Samples PCA'
  ) + theme(axis.text = element_text(size = 12), axis.title = element_text(size = 15))
}
d_p(pg.pef.tpm[pick.genes,c(2:9)],group_list_pg)

# Fig.1f -------
library(EnhancedVolcano)

## SF
options(scipen = 10)
cancer.SF <- c('SRSF1','ENSSSCG00000036592','SRSF3','SRSF5','SRSF6','SRSF10','TRA2B','ENSSSCG00000000288','ENSSSCG00000036350','HNRNPH1','HNRNPK','HNRNPL','HNRNPM','ENSSSCG00000013421','RBM5','RBM10','RBM39','ESRP1','RBFOX2','QKI')

cancer.SF.tpm <- tpm[cancer.SF,]
cancer.SF.tpm.pg <- pg.pef.tpm[cancer.SF,]
library(ggrastr)

EnhancedVolcano(pg.high_pef,
                lab = pg.high_pef$geneid,
                selectLab = cancer.SF,
                x = 'log2FoldChange',
                y = 'padj',
                FCcutoff = 2,
                pCutoff = 0.05,
                ylab = bquote(~-Log[10]~ 'padj'),
                pointSize = 1.0,
                labSize = 2.0,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                parseLabels = T,
                col = c('gray80', '#2a9d8f', '#e9c46a', 'red3'),
                colAlpha = 2/5,
                legendPosition = 'bottom',
                legendLabSize = 12,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = .5,
                max.overlaps = 300,
                colConnectors = '#14213d',
                title = 'Cancer Splice Factors',
                subtitle = 'P.H vs PEF',
                raster = T,
                arrowheads = F
) +
  ggplot2::coord_cartesian(xlim=c(-10, 10)) +
  ggplot2::scale_x_continuous(breaks=c(seq(-10,-2,4), seq(-2,2,1),seq(2,10,4))) 
  # coord_flip()
  # scale_x_log10() +
  # scale_y_log10()


EnhancedVolcano(pg.low_pef,
                lab = pg.low_pef$geneid,
                selectLab = cancer.SF,
                x = 'log2FoldChange',
                y = 'padj',
                FCcutoff = 2,
                pCutoff = 0.05,
                ylab = bquote(~-Log[10]~ 'padj'),
                pointSize = 1.0,
                labSize = 2.0,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                parseLabels = T,
                col = c('gray80', '#2a9d8f', '#e9c46a', 'red3'),
                colAlpha = 2/5,
                legendPosition = 'bottom',
                legendLabSize = 10,
                legendIconSize = 3.0,
                drawConnectors = TRUE,
                widthConnectors = .5,
                max.overlaps = 300,
                colConnectors = '#14213d',
                title = 'Cancer Splice Factors',
                subtitle = 'P.L vs PEF',
                raster = T,
                arrowheads = F
                )+
  ggplot2::coord_cartesian(xlim=c(-10, 10)) +
  ggplot2::scale_x_continuous(breaks=c(seq(-10,-2,4), seq(-2,2,1),seq(2,10,4))) 
# coord_flip()
