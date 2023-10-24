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

