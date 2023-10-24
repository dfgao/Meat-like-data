# Fig.S3a --------
pg.tpm_cor <- cor(pg.pef.tpm[,c(2:9)],method = 'pearson',use = 'pairwise.complete.ob')

pheatmap(pg.tpm_cor,
         border_color = NA, 
         clustering_method = 'ward.D',
         color = scales::alpha(colorRampPalette(colors = c('#00509d','gray80','#f35b04'),alpha=T,bias=1)(256),alpha = 1),
         angle_col = '315',
         display_numbers = T,
         annotation_col = ann,
         annotation_colors = ha_top.col,
         number_color = 'white')

# Fig.S3b -------
group_list_pg <- factor(c( rep('pgEpiSCs.Low',3), rep('pgEpiSCs.High',3),rep('PEF',2)),levels = c('pgEpiSCs.Low','pgEpiSCs.High','PEF'))

d_p <- function(tpm,group_list){
  tpm.t = as.data.frame(t(tpm))
  tpm.t = cbind(tpm.t,group_list)
  tpm.pca <- PCA(tpm.t[,-ncol(tpm.t)],graph = FALSE)
  fviz_screeplot(tpm.pca, addlabels = TRUE, ylim = c(0, 70),title='Dim choose')
}
d_p(pg.pef.tpm[pick.genes,c(2:9)],group_list_pg)

# Fig.S3c ------
# prcomp score
tpm.t <- data.frame(t(pg.pef.tpm[pick.genes,2:9]))
a <- prcomp(tpm.t,scale=TRUE)
pca.data <- data.frame(Sample=rownames(a$x),
                       X=a$x[,1],
                       Y=a$x[,2])
pca.var <- a$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100,1)
ggplot(data=pca.data, aes(x=X, y=Y, label=Sample))+
  geom_text() +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep=""))+
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep=""))+
  theme_bw()+
  ggtitle("My PCA Graph")

# only PC1
prcomp.score <- a$rotation[,1]
prcomp.score <- sort(prcomp.score, decreasing = F)
prcomp.score <- as.data.frame(prcomp.score)
colnames(prcomp.score) <- 'PC1'
prcomp.score$Rank <- seq(1,nrow(prcomp.score))

prcomp.score <- prcomp.score %>%
  rownames_to_column(var="Symbol") %>% 
  mutate(highlight = "")

prcomp.score[c(na.omit(match(rownames(tpm.plu), prcomp.score$Symbol))), 4] <- rownames(tpm.plu)
p1 <- ggplot(prcomp.score, aes(x=Rank, y=PC1)) +
  geom_point_rast() +
  geom_text_repel(
    aes(label = highlight),
    size = 3,
    min.segment.length = 0,
    seed = 100,
    box.padding = 0.5,
    max.overlaps = Inf,
    arrow = arrow(length = unit(0.010, "npc")),
    nudge_x = 800,
    nudge_y = .001,
    color = "grey30"
  ) + 
  geom_rug(col="steelblue",alpha=0.1, size=.1,sides = 'l') +
  theme_bw()
p1

# Fig.S3d & S3g------
## pg.low-PEF
group = pg.pef.count[pick.genes,c(2:4,8:9)]
group_list.use <- c(rep('pgEpiSCs.Low',3),rep('PEF',2))
colData = data.frame(row.names = colnames(group),group_list = group_list.use)

dds <- DESeqDataSetFromMatrix(countData = group,colData = colData,design = ~group_list)
dds <- DESeq(dds)
plotDispEsts(dds)
RES <- results(dds, contrast = c('group_list','pgEpiSCs.Low','PEF'))

pg.low_pef <- RES %>%
  data.frame() %>%
  rownames_to_column(var="geneid") %>%
  as_tibble() %>%
  arrange(padj)
pg.low_pef$change <- factor(ifelse(pg.low_pef$padj < 0.05 ,
                                   ifelse(pg.low_pef$log2FoldChange > 2,'pg.Low.UP', 
                                          ifelse(pg.low_pef$log2FoldChange < -2,'pg.Low.DOWN','NOT change')),'NOT change'))
pg.low_pef.sig <- dplyr::filter(pg.low_pef, padj < 0.05) %>%
  dplyr::arrange(padj)

## pg.high-PEF
group = pg.pef.count[pick.genes,c(5:9)]
group_list.use <- c(rep('pgEpiSCs.High',3),rep('PEF',2))
colData = data.frame(row.names = colnames(group),group_list = group_list.use)

dds <- DESeqDataSetFromMatrix(countData = group,colData = colData,design = ~group_list)
dds <- DESeq(dds)
plotDispEsts(dds)
RES <- results(dds, contrast = c('group_list','pgEpiSCs.High','PEF'))

pg.high_pef <- RES %>%
  data.frame() %>%
  rownames_to_column(var="geneid") %>%
  as_tibble() %>%
  arrange(padj)
pg.high_pef$change <- factor(ifelse(pg.high_pef$padj < 0.05 ,
                                   ifelse(pg.high_pef$log2FoldChange > 2,'pg.High.UP', 
                                          ifelse(pg.high_pef$log2FoldChange < -2,'pg.High.DOWN','NOT change')),'NOT change'))
pg.high_pef.sig <- dplyr::filter(pg.high_pef, padj < 0.05) %>%
  dplyr::arrange(padj)


## venn for LFC > 2 
DEGvenn = list(pgEpiSC.Low_pEF = pg.low_pef.sig[(pg.low_pef.sig$log2FoldChange > 2 | pg.low_pef.sig$log2FoldChange < -2),]$geneid, 
               pgEpiSC.High_pEF = pg.high_pef.sig[(pg.high_pef.sig$log2FoldChange > 2 | pg.high_pef.sig$log2FoldChange < -2),]$geneid)
DEGvenn_res = Venn(DEGvenn)
plot(DEGvenn_res,doWeights = F,type="circles")

## get orth

pg.low_pef.sig <- pg.low_pef.sig %>% 
  left_join(y = pig.gene.info[,c(1:2)], by = c('geneid' = 'genename'))
pg.low_pef.sig <- pg.low_pef.sig %>% 
  left_join(y = orth[,c(1:4)], by = c('geneid.y' = 'geneid'))

pg.high_pef.sig <- pg.high_pef.sig %>% 
  left_join(y = pig.gene.info[,c(1:2)], by = c('geneid' = 'genename'))
pg.high_pef.sig <- pg.high_pef.sig %>% 
  left_join(y = orth[,c(1:4)], by = c('geneid.y' = 'geneid'))

pg.low_pef <- pg.low_pef %>% 
  left_join(y = pig.gene.info[,c(1:2)], by = c('geneid' = 'genename'))
pg.low_pef <- pg.low_pef %>% 
  left_join(y = orth[,c(1:4)], by = c('geneid.y' = 'geneid'))

pg.high_pef <- pg.high_pef %>% 
  left_join(y = pig.gene.info[,c(1:2)], by = c('geneid' = 'genename'))
pg.high_pef <- pg.high_pef %>% 
  left_join(y = orth[,c(1:4)], by = c('geneid.y' = 'geneid'))

DEGs <- list(pg.low_pef = pg.low_pef, pg.high_pef = pg.high_pef)
export(DEGs,file = '../../part4_new/DEGs.xlsx',row.names=T)



## venn for down FC 2
library(venn)
DEGvenn_down = list(pgEpiSC.Low_down = pg.low_pef.sig[(pg.low_pef.sig$log2FoldChange < -2 ),]$geneid, 
               pgEpiSC.High_down = pg.high_pef.sig[(pg.high_pef.sig$log2FoldChange < -2 ),]$geneid)
DEGvenn_down_res = Venn(DEGvenn_down)
plot(DEGvenn_down_res,doWeights = F,type="circles")

## venn for up FC 2
DEGvenn_up = list(pgEpiSC.Low_up = pg.low_pef.sig[(pg.low_pef.sig$log2FoldChange > 2 ),]$geneid, 
                    pgEpiSC.High_up = pg.high_pef.sig[(pg.high_pef.sig$log2FoldChange > 2 ),]$geneid)
DEGvenn_up_res = Venn(DEGvenn_up)
plot(DEGvenn_up_res,doWeights = F,type="circles")

venn(DEGvenn_up, 
     zcolor=sample.col[1:2],
     opacity = .7,
     box=F,
     sncs=1.2,
     ilcs=1.2)

DEG.res <- list(DEGvenn_down_pg.low=DEGvenn_down_res@IntersectionSets[["10"]], 
                DEGvenn_down_pg.high=DEGvenn_down_res@IntersectionSets[["01"]], 
                DEGvenn_down_common=DEGvenn_down_res@IntersectionSets[["11"]], 
                DEGvenn_up_pg.low=DEGvenn_up_res@IntersectionSets[["10"]], 
                DEGvenn_up_pg.high=DEGvenn_up_res@IntersectionSets[["01"]], 
                DEGvenn_up_common=DEGvenn_up_res@IntersectionSets[["11"]]
                )

export(DEG.res,file = '../../part4_new/DEG.res.change.LFC2.xlsx',row.names=T)
export(pg.pef.tpm,file = '../../part4_new/pg.pef.tpm.xlsx',row.names=T)


## venn for down FC 2 for orth
DEGvenn_down = list(pgEpiSC.Low_down = na.omit(pg.low_pef.sig[(pg.low_pef.sig$log2FoldChange < -2 ),]$human.genename), 
                    pgEpiSC.High_down = na.omit(pg.high_pef.sig[(pg.high_pef.sig$log2FoldChange < -2 ),]$human.genename))
DEGvenn_down_res = Venn(DEGvenn_down)
plot(DEGvenn_down_res,doWeights = F,type="circles")

## venn for up FC 2 for orth
DEGvenn_up = list(pgEpiSC.Low_up = na.omit(pg.low_pef.sig[(pg.low_pef.sig$log2FoldChange > 2 ),]$human.genename), 
                  pgEpiSC.High_up = na.omit(pg.high_pef.sig[(pg.high_pef.sig$log2FoldChange > 2 ),]$human.genename))
DEGvenn_up_res = Venn(DEGvenn_up)
plot(DEGvenn_up_res,doWeights = F,type="circles")

DEG.res <- list(DEGvenn_down_pg.low=DEGvenn_down_res@IntersectionSets[["10"]], 
                DEGvenn_down_pg.high=DEGvenn_down_res@IntersectionSets[["01"]], 
                DEGvenn_down_common=DEGvenn_down_res@IntersectionSets[["11"]], 
                DEGvenn_up_pg.low=DEGvenn_up_res@IntersectionSets[["10"]], 
                DEGvenn_up_pg.high=DEGvenn_up_res@IntersectionSets[["01"]], 
                DEGvenn_up_common=DEGvenn_up_res@IntersectionSets[["11"]]
)

export(DEG.res,file = '../../part4_new/DEG.res.change.orth.xlsx',row.names=T)

pheatmap(
  t(scale(t(pg.pef.tpm[unique(c(pg.low_pef.sig[(pg.low_pef.sig$log2FoldChange < -2 ),]$geneid, 
                                pg.high_pef.sig[(pg.high_pef.sig$log2FoldChange < -2 ),]$geneid, 
                                pg.low_pef.sig[(pg.low_pef.sig$log2FoldChange > 2 ),]$geneid,
                                pg.high_pef.sig[(pg.high_pef.sig$log2FoldChange > 2 ),]$geneid)),c(2:9)]))),
  border_color = NA, 
  clustering_method = 'ward.D',
  annotation_col = ann,
  annotation_colors = ha_top.col,
  color = scales::alpha(rev(colorRampPalette(RColorBrewer::brewer.pal(8,"RdBu"),alpha=T,bias=1)(256)),alpha = 1),
  angle_col = '315',
  show_rownames = F)

# Fig.S3e --------
library(ggtern)

for.tern.tpm <- data.frame(pgEpiSC.low = matrixStats::rowMeans2(as.matrix(pg.pef.tpm[pick.genes, 2:4])),
                           pgEpiSC.high = matrixStats::rowMeans2(as.matrix(pg.pef.tpm[pick.genes, 5:7])),
                           pEF = matrixStats::rowMeans2(as.matrix(pg.pef.tpm[pick.genes, 8:9])))

for.tern.tpm$gene <- pick.genes
for.tern.tpm$average <- rowMeans(for.tern.tpm[,1:3])
for.tern.tpm$gene.type <- 'Other'
for.tern.tpm[c(na.omit(match(rownames(tpm.plu), for.tern.tpm$gene))), 6] <- 'Pluripotency'

for.tern.tpm$gene.type <- factor(for.tern.tpm$gene.type,levels = c('Pluripotency','Other'))
for.tern.tpm <- for.tern.tpm[order(for.tern.tpm$gene.type),]

ggtern(data=for.tern.tpm,aes(x=pgEpiSC.low,y=pgEpiSC.high,z=pEF))+    
  geom_point_rast(data = for.tern.tpm[for.tern.tpm$gene.type == 'Other',], aes(size = average), color = "grey80",
             alpha = 0.8, show.legend = FALSE) +
  theme_rgbw(base_size = 12 )+   
  labs(title = "3 cell type genes")+ 
  theme(plot.title = element_text(size=15,hjust = 0.5)) +
  theme_showarrows() + 
  stat_density_tern(aes(fill=..level.., alpha=..level..),geom='polygon',bdl = 0.08, h = 2,alpha=0.8,expand=0.75,bins = 6) +
  scale_fill_gradient2(high = "#71D0F5FF") +
  geom_point(data = for.tern.tpm[for.tern.tpm$gene.type == 'Pluripotency',], aes(size = average, color = '#B2182B'), 
             show.legend = F)+    
  geom_mask() 

# Fig.S3f -------
can.genes <- c('SERPINE1','	MT2A','UBC','ACTB','TUBA1B','HMGB2','EIF2S3','RPS6','PCNA','STMN1','HMGB1','H2AFZ','VIM','FN1','FOS','DDIT4','NFKBIA','JUNB','TUBB','PKM','LGALS3','BIRC5','TOP2A','NPM1','RPL35','HNRNPH1','SOX4','SOX11','VEGFA','EGFR','MYC','CD44','FOXR2','TP53')

tpm.can <- na.omit(pg.pef.tpm[can.genes,-c(1,ncol(pg.pef.tpm))])
pheatmap(t(scale(t(tpm.can))),border_color = NA, clustering_method = 'ward.D2',cluster_cols = F,color = scales::alpha(rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"),alpha=T,bias=.8)(256)),alpha = 1),angle_col = '315',breaks = unique(c(seq(-2,2, length=256))))

EnhancedVolcano(pg.low_pef,
                lab = pg.low_pef$geneid,
                selectLab =can.genes,
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
                title = 'Cancer Hallmarkers Genes',
                subtitle = 'P.L vs PEF',
                arrowheads = F,
                raster = T
) +
  ggplot2::coord_cartesian(xlim=c(-10, 10)) +
  ggplot2::scale_x_continuous(breaks=c(seq(-10,-2,4), seq(-2,2,1),seq(2,10,4))) 
# coord_flip()

EnhancedVolcano(pg.high_pef,
                lab = pg.high_pef$geneid,
                selectLab =can.genes,
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
                title = 'Cancer Hallmarkers Genes',
                subtitle = 'P.H vs PEF',
                arrowheads = F,
                raster = T
) +
  ggplot2::coord_cartesian(xlim=c(-10, 10)) +
  ggplot2::scale_x_continuous(breaks=c(seq(-10,-2,4), seq(-2,2,1),seq(2,10,4))) 

# Fig.S3h --------
go.bp <- import_list('../../part4_new/LFC2_gobp.xlsx')
go.bp.inte <- rbind(go.bp$PEF.GOBP,go.bp$PG.GOBP)
go.bp.inte$Z_score <- (go.bp.inte$Count -mean(go.bp.inte$Count))/sd(go.bp.inte$Count)
go.bp.inte$Type <- c(rep('PEF',20),rep('PG',20))
go.bp.inte$GeneRatio <- go.bp.inte$Count/go.bp.inte$Background

go.vis <- go.bp.inte[1:20,1:2] #pEF
colnames(go.vis) <- c('Terms','LogP')
go.vis$lp <- -go.vis$LogP
go.vis <- go.vis[order(go.vis$lp),]
go.vis$Terms <- factor(go.vis$Terms,levels = go.vis$Terms)

ggplot(go.vis, aes(x=Terms, y=lp)) +
  geom_segment( aes(x=Terms, xend=Terms, y=0, yend=lp), color="#e9c46a",size=2,lty = 1) +
  geom_point( color="#e63946", size=4, alpha=0.6) +
  theme_light() +
  coord_flip() + 
  xlab( 'Top 20 GO BP cluster terms') +
  ylab(bquote(~-Log[10]~ 'P'))+
  theme(
    axis.text = element_text(size = 12,colour = 'black'),
    axis.title.x = element_text(size = 15,colour = 'black'),
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank()
  )


go.vis <- go.bp.inte[21:40,1:2] #pg
colnames(go.vis) <- c('Terms','LogP')
go.vis$lp <- -go.vis$LogP
go.vis <- go.vis[order(go.vis$lp),]
go.vis$Terms <- factor(go.vis$Terms,levels = go.vis$Terms)

ggplot(go.vis, aes(x=Terms, y=lp)) +
  geom_segment( aes(x=Terms, xend=Terms, y=0, yend=lp), color="skyblue",size=2,lty = 1) +
  geom_point( color="blue", size=4, alpha=0.6) +
  theme_light() +
  coord_flip() + 
  xlab( 'Top 20 GO BP cluster terms') +
  ylab(bquote(~-Log[10]~ 'P'))+
  theme(
    axis.text = element_text(size = 12,colour = 'black'),
    axis.title.x = element_text(size = 15,colour = 'black'),
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank()
  )

