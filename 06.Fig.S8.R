# Fig.S8
cluster.tpm.tmp <- pg.pgmc.tpm.clean.cluster[c(pg.pgmc.tpm.clean.cluster$cluster=='C_1(1232 genes)' | pg.pgmc.tpm.clean.cluster$cluster=='C_7(1343 genes)'),c(2:9)]
cluster.tpm.tmp <- pg.pgmc.tpm.clean.cluster[c(pg.pgmc.tpm.clean.cluster$cluster=='C_2(1523 genes)' | pg.pgmc.tpm.clean.cluster$cluster=='C_3(1253 genes)'),c(2:9)]
cluster.tpm.tmp <- pg.pgmc.tpm.clean.cluster[c(pg.pgmc.tpm.clean.cluster$cluster=='C_5(1349 genes)' | pg.pgmc.tpm.clean.cluster$cluster=='C_8(1830 genes)'),c(2:9)]
cluster.tpm.tmp <- pg.pgmc.tpm.clean.cluster[c(pg.pgmc.tpm.clean.cluster$cluster=='C_10(2666 genes)' ),c(2:9)]

pheatmap(t(scale(t(cluster.tpm.tmp))), 
         border_color = NA, 
         clustering_method = 'ward.D',
         show_rownames = F,
         cluster_cols = F,cluster_rows = T,
         color = scales::alpha(rev(colorRampPalette(RColorBrewer::brewer.pal(5,"RdBu"),alpha=T,bias=.8)(256)),alpha = .8),
         angle_col = '315',
         # annotation_row = pick.genes.for.cluster.tpm[,c(18,19)],
         breaks = unique(c(seq(-2,2, length=256))),fontsize_row = 8,
)

mc.gobp <- import('../../MSC/kmeans/cluster.top20.GOBP.xlsx')
mc.gobp$GeneRatio <- mc.gobp$Count/mc.gobp$Background
mc.gobp.tmp <- mc.gobp[1:20,]
mc.gobp.tmp2 <- mc.gobp.tmp[order(mc.gobp.tmp$GeneRatio),]
mc.gobp.tmp2$Description <- factor(mc.gobp.tmp2$Description, levels = mc.gobp.tmp2$Description)

ggplot(data = mc.gobp.tmp2, aes(x = GeneRatio, y = Description)) + 
  geom_point(aes(size = Count,color = -LogQ)) +
  scale_color_gradient(low = col[1],high = col[3], name=expression(-log[10](padj))) +
  scale_size(range  =  c(0, 6)) +
  labs( y = 'Top 20 GO BP cluster 1/7') +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10,face = "bold")) 
