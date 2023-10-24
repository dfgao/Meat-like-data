# Fig.S1a --------
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
