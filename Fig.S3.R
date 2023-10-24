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
