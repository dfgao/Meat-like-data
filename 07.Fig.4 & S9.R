# Fig.4c ------
library(mixOmics)
library(pheatmap)
library(ComplexHeatmap)
library(tidyverse)
library(rio)

# data input------
meta_intensity_all <- import('~/project/01.zhugaoxiang/04.pig.RNA.MSC/17_metabolome_report_0402/report1/3.Metabolite_identify/meta_intensity_all.clean.txt')
rownames(meta_intensity_all) <- meta_intensity_all$name

meta.ratio.pc_pgmc <- read.csv('~/project/01.zhugaoxiang/04.pig.RNA.MSC/17_metabolome_report_0402/report1/5.Metabolite_Diff/pgEpiSCs.Low_pgEpiSCs.MC/pgEpiSCs.Low.vs.pgEpiSCs.MC_all_All.csv',header = T,check.names = F,allowEscapes = T,fill = T)
rownames(meta.ratio.pc_pgmc) <- meta.ratio.pc_pgmc$name

meta.ratio.pc_pgmc.diff <- read.csv('~/project/01.zhugaoxiang/04.pig.RNA.MSC/17_metabolome_report_0402/report1/5.Metabolite_Diff/pgEpiSCs.Low_pgEpiSCs.MC/pgEpiSCs.Low.vs.pgEpiSCs.MC_all_Diff.csv',header = T,check.names = F,allowEscapes = T,fill = T)

library(EnhancedVolcano)
meta.ratio.pc_pgmc$log2FC <- as.numeric(meta.ratio.pc_pgmc$log2FC)
EnhancedVolcano(meta.ratio.pc_pgmc,
                lab = meta.ratio.pc_pgmc$ID,
                # selectLab = cancer.SF,
                x = 'log2FC',
                y = 'p.value',
                FCcutoff = 0.263,
                pCutoff = 0.05,
                # ylab = bquote(~-Log[10]~ 'padj'),
                pointSize = 2.0,
                labSize = 3.0,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                parseLabels = T,
                col = c('gray80', '#2a9d8f', '#e9c46a', 'red3'),
                colAlpha = 4/5,
                legendPosition = 'bottom',
                legendLabSize = 10,
                legendIconSize = 5.0,
                drawConnectors = TRUE,
                widthConnectors = .5,
                max.overlaps = 10,
                colConnectors = '#14213d',
                title = 'pgEpiSCs vs pg.MC',
                subtitle = NULL,
) 

# Fig.4d ------
test <- meta.ratio.pc_pgmc.diff[,7:18]
colnames(test) <- rownames(ann)

pheatmap(t(scale(t(test))),
         border_color = NA, 
         clustering_method = 'ward.D',
         color = scales::alpha(rev(colorRampPalette(RColorBrewer::brewer.pal(5,"RdBu"),alpha=T,bias=1)(256)),alpha = .8),
         angle_col = '315',
         display_numbers = F,
         annotation_col = ann,
         annotation_colors = ha_top.col,
         number_color = 'white')

# Fig.4e ------
pg.down.kegg <- import('~/project/01.zhugaoxiang/04.pig.RNA.MSC/17_metabolome_report_0402/report1/5.Metabolite_Diff/pgEpiSCs.Low_pgEpiSCs.MC/KEGG/pgEpiSCs.Low.vs.pgEpiSCs.MC_all_DOWN_kegg_enrichment.txt',sep = '\t',header = T,check.names = F,allowEscapes = T,fill = T)

plot.data <- pg.down.kegg[pg.down.kegg$EnrichDirect == 'Over',]%>% slice_max(order_by = x, n = 10)
plot.data$ratio <- plot.data$x/plot.data$y
plot.data <- plot.data[,c(2,11)]

colnames(plot.data) <- c('Terms','Ratio')

plot.data <- plot.data[order(plot.data$Ratio),]
plot.data$Terms <- factor(plot.data$Terms,levels = plot.data$Terms)

ggplot(plot.data, aes(x=Terms, y=Ratio)) +
  geom_segment( aes(x=Terms, xend=Terms, y=0, yend=Ratio), color="#F77F00",size=2,lty = 1) +
  geom_point( color="red3", size=4, alpha=0.6) +
  theme_light() +
  coord_flip() +
  xlab("") + ylab('Ratio')+
  theme(
    axis.text = element_text(size = 12,colour = 'black'),
    axis.title.x = element_text(size = 15,colour = 'black'),
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank()
  )

# Fig.4f ------
diff.genes <- c('PGK1','PGAM1','ENSSSCG00000001930','LDHA','PDHA1','PDK4','PDP1','HK1','HK2','SLC2A11','SLC2A12','SLC2A13','SLC16A3','PGM1','GPI')
diff.genes.tpm <- clean_tpm[match(diff.genes, rownames(clean_tpm)),] %>% na.omit()

ha_left <- list(gene.type = npg.col[13:15],Group = c(col[c(1,3)]))
names(ha_left$Group) = c('pgEpiSC','pgEpiSC.MC')

gene.type.diff <- data.frame(gene.type = c(rep('ATP-generating',4), rep('Acetyl-CoA generating',3), rep('Oxidation and Glycolysis',8)))
rownames(gene.type.diff) <- rownames(diff.genes.tpm)
names(ha_left$gene.type) = unique(gene.type.diff$gene.type)

annotation_col_matrix.mc.diff <- data.frame(Group = annotation_col_matrix.mc[-c(4,5,9:11),])
rownames(annotation_col_matrix.mc.diff) <- rownames(annotation_col_matrix.mc)[c(1:3,6:8)]

ComplexHeatmap::pheatmap(t(scale(t(diff.genes.tpm[,c(10:12,13:15)]))),
                         border_color = NA,
                         clustering_method = 'median',
                         show_rownames = T,
                         cluster_cols = F,cluster_rows = F,
                         color = scales::alpha(rev(colorRampPalette(RColorBrewer::brewer.pal(5,"RdBu"),alpha=T,bias=1)(256)),alpha = .8),
                         angle_col = '315',
                         annotation_row = gene.type.diff,
                         annotation_col = annotation_col_matrix.mc.diff,
                         annotation_colors = ha_left,
                         # breaks = unique(c(seq(0,12, length=256))),
                         fontsize_row = 10
)

# Fig.4g -------
tca <- c('Succinic acid')
tca.ratio <- meta.ratio.pc_pgmc.diff[match(tca,meta.ratio.pc_pgmc.diff$name), c(2,7:18)] %>% na.omit()

test <- melt(tca.ratio,id.vars = 'name')
test$Cell.type <- factor(c(rep('pgEpiSCs',6),rep('pg.MC',6)), levels = c('pgEpiSCs','pg.MC'))

ggplot(test, aes(x=name, y=value, fill=Cell.type)) + 
  geom_boxplot() +
  labs(x = 'Compounds', y = 'Abundance ratio') + 
  scale_fill_manual(values = c("#2a9d8f", "#e63946")) +
  # geom_jitter(color="black", size=0.4, alpha=0.9) +
  facet_wrap(~name, scale="free") + 
  theme_bw() +
  theme(axis.title = element_text(size = 13)) +
  ggtitle('TCA cycle')

glu <- c('Glucose 6-phosphate')
glu.ratio <- meta.ratio.pc_pgmc.diff[match(glu,meta.ratio.pc_pgmc.diff$name), c(2,7:18)] %>% na.omit()

test <- melt(glu.ratio,id.vars = 'name')
test$Cell.type <- factor(c(rep('pgEpiSCs',36),rep('pg.MC',36)), levels = c('pgEpiSCs','pg.MC'))

ggplot(test, aes(x=name, y=value, fill=Cell.type)) + 
  geom_boxplot() +
  labs(x = 'Compounds', y = 'Abundance ratio') + 
  scale_fill_manual(values = c("#2a9d8f", "#e63946")) +
  # geom_jitter(color="black", size=0.4, alpha=0.9) +
  facet_wrap(~name, scale="free") + 
  theme_bw() +
  theme(axis.title = element_text(size = 13)) +
  ggtitle('Glycolysis')

