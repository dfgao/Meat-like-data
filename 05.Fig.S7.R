# Fig.S7a -------
kmeans.info <- data.frame(Gene.type=c('Exp.PCG','Not.exp.PCG'), count = c(14546,6734))
kmeans.info$fraction = kmeans.info$count / sum(kmeans.info$count)
kmeans.info$ymax = cumsum(kmeans.info$fraction)
kmeans.info$ymin = c(0, head(kmeans.info$ymax, n=-1))
kmeans.info$labelPosition <- (kmeans.info$ymax + kmeans.info$ymin) / 2
kmeans.info$label <- paste0(kmeans.info$Gene.type, "\n number: ", kmeans.info$count)

ggplot(kmeans.info, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Gene.type)) +
  geom_rect() +
  geom_label( x=3.5, aes(y=labelPosition, label=label), size=5) +
  # scale_fill_manual(values = c('#e63946','#2a9d8f'),alpha = .2) +
  scale_fill_simpsons() +
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  theme_void() +
  theme(legend.position = "none")

# Fig.S7b -------
group_list_pg <- factor(c( rep('pgEpiSCs.Low',3), rep('pgEpiSCs.MPC',2),rep('pgEpiSCs.MC',3)),levels = c('pgEpiSCs.Low','pgEpiSCs.MPC','pgEpiSCs.MC'))

d_p <- function(tpm,group_list){
  tpm.t = as.data.frame(t(tpm))
  tpm.t = cbind(tpm.t,group_list)
  tpm.pca <- PCA(tpm.t[,-ncol(tpm.t)],graph = FALSE)
  fviz_screeplot(tpm.pca, addlabels = TRUE, ylim = c(0, 70),title='Dim choose')
}
d_p(pg.pgmc.tpm.clean,group_list_pg)

# Fig.S7e -------
library(ComplexHeatmap)
library(circlize)

gene_pos <- match(rownames(diff.genes.tpm),pg.pgmc.tpm.clean.cluster$gene) %>% na.omit()
row_anno <- rowAnnotation(gene=anno_mark(at=gene_pos,labels = rownames(diff.genes.tpm)))
rownames(pg.pgmc.tpm.clean.cluster) <- pg.pgmc.tpm.clean.cluster$gene

annotation_row_matrix <- data.frame(pg.pgmc.tpm.clean.cluster[,c(10)])
rownames(annotation_row_matrix) <- rownames(pg.pgmc.tpm.clean.cluster)
colnames(annotation_row_matrix) <- 'Kmeans.cluster'
annotation_col_matrix <- data.frame(Group=c(rep('pgEpiSC',3), rep('pgEpiSC.MPC',2), rep('pgEpiSC.MC',3)))
rownames(annotation_col_matrix) = colnames(pg.pgmc.tpm.clean.cluster)[2:9]

show_col(pal_npg("nrc")(10))
ha_left <- list(Kmeans.cluster=pal_npg("nrc")(10), Group = col)
names(ha_left$Kmeans.cluster) = levels(annotation_row_matrix$Kmeans.cluster)
names(ha_left$Group) = c('pgEpiSC','pgEpiSC.MPC','pgEpiSC.MC')

ComplexHeatmap::pheatmap(t(scale(t(pg.pgmc.tpm.clean.cluster[,c(2:9)]))),
         border_color = NA, 
         clustering_method = 'ward.D2',
         show_rownames = F,
         cluster_cols = F,cluster_rows = T,
         color = scales::alpha(rev(colorRampPalette(RColorBrewer::brewer.pal(5,"RdBu"),alpha=T,bias=.8)(256)),alpha = .8),
         angle_col = '315',
         annotation_row = annotation_row_matrix,
         annotation_colors = ha_left,
         annotation_col = annotation_col_matrix,
         breaks = unique(c(seq(-2,2, length=256))),
         right_annotation = row_anno)

# Fig.S7g ------
data.clean <- data.clean %>% rownames_to_column(var = 'geneid') %>% left_join( y = pig.gene.info[,c(1,2)], by = 'geneid')
rownames(data.clean) <- data.clean$genename

group = data.clean[mc.pick.genes,c(11:13,5:7)]
group_list.use <- c(rep('pgEpiSC',3),rep('MC',3))
colData = data.frame(row.names = colnames(group),group_list = group_list.use)

dds <- DESeqDataSetFromMatrix(countData = group,colData = colData,design = ~group_list)
dds <- DESeq(dds)
plotDispEsts(dds)
RES <- results(dds, contrast = c('group_list','MC','pgEpiSC'))

pg.low_mc <- RES %>%
  data.frame() %>%
  rownames_to_column(var="geneid") %>%
  as_tibble() %>%
  arrange(padj)
pg.low_mc$change <- factor(ifelse(pg.low_mc$padj < 0.05, ifelse(pg.low_mc$log2FoldChange > 0,'MC.UP','MC.DOWN'),'NOT change'))
pg.low_mc.sig <- dplyr::filter(pg.low_mc, padj < 0.05) %>%
  dplyr::arrange(padj)

## pg.low-pg.mc
group = data.clean[mc.pick.genes,c(11:16)]
group_list.use <- c(rep('pgEpiSC',3),rep('pg.MC',3))
colData = data.frame(row.names = colnames(group),group_list = group_list.use)

dds <- DESeqDataSetFromMatrix(countData = group,colData = colData,design = ~group_list)
dds <- DESeq(dds)
plotDispEsts(dds)
RES <- results(dds, contrast = c('group_list','pg.MC','pgEpiSC'))

pg.low_pgmc <- RES %>%
  data.frame() %>%
  rownames_to_column(var="geneid") %>%
  as_tibble() %>%
  arrange(padj)
pg.low_pgmc$change <- factor(ifelse(pg.low_pgmc$padj < 0.05, ifelse(pg.low_pgmc$log2FoldChange > 0,'pgMC.UP','pgMC.DOWN'),'NOT change'))
pg.low_pgmc.sig <- dplyr::filter(pg.low_pgmc, padj < 0.05) %>%
  dplyr::arrange(padj)


## venn for LFC > 1.5 
DEGvenn = list(pgEpiSC.MC = pg.low_mc.sig[(pg.low_mc.sig$log2FoldChange > log2(1.5) | pg.low_mc.sig$log2FoldChange < -log2(1.5)),]$geneid, 
               pgEpiSC.pgMC = pg.low_pgmc.sig[(pg.low_pgmc.sig$log2FoldChange > log2(1.5) | pg.low_pgmc.sig$log2FoldChange < -log2(1.5)),]$geneid)
DEGvenn_res = Venn(DEGvenn)
plot(DEGvenn_res,doWeights = T,type="circles")

## get orth

pg.low_mc.sig <- pg.low_mc.sig %>% 
  left_join(y = pig.gene.info[,c(1:2)], by = c('geneid' = 'genename'))
pg.low_mc.sig <- pg.low_mc.sig %>% 
  left_join(y = orth[,c(1:4)], by = c('geneid.y' = 'geneid'))

pg.low_pgmc.sig <- pg.low_pgmc.sig %>% 
  left_join(y = pig.gene.info[,c(1:2)], by = c('geneid' = 'genename'))
pg.low_pgmc.sig <- pg.low_pgmc.sig %>% 
  left_join(y = orth[,c(1:4)], by = c('geneid.y' = 'geneid'))

DEG.sig <- list(pg.low_mc.sig = pg.low_mc.sig, pg.low_pgmc.sig = pg.low_pgmc.sig)
export(DEG.sig,file = '../../MSC/kmeans/pg-MC-pgMC,deg.sig.xlsx',row.names=T)
