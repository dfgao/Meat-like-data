# Fig.3a -----
pg.mc.tpm <- tpm[,-c(8:10)]
pg.pgmc.tpm <- tpm[,c(1,11:19)]
mc.pick.genes <- unique(c(rownames(pg.pgmc.tpm[(matrixStats::rowMedians(as.matrix(pg.pgmc.tpm[,c(2:(ncol(pg.pgmc.tpm)-1))]),cols = 1:3) > 0.5) ,]), rownames(pg.pgmc.tpm[(matrixStats::rowMedians(as.matrix(pg.pgmc.tpm[,c(2:(ncol(pg.pgmc.tpm)-1))]),cols = 4:6) > 0.5) ,]), rownames(pg.pgmc.tpm[(matrixStats::rowMedians(as.matrix(pg.pgmc.tpm[,c(2:(ncol(pg.pgmc.tpm)-1))]),cols = 7:8) > 0.5) ,])))
pg.pgmc.tpm.clean <- pg.pgmc.tpm[mc.pick.genes,c(2:(ncol(pg.pgmc.tpm)-1))]
pg.pgmc.tpm.clean <- pg.pgmc.tpm.clean[,c(1:3,7,8,4:6)]

group_list_pg <- factor(c( rep('pgEpiSCs.Low',3), rep('pgEpiSCs.MPC',2),rep('pgEpiSCs.MC',3)),levels = c('pgEpiSCs.Low','pgEpiSCs.MPC','pgEpiSCs.MC'))

d_p <- function(tpm,group_list){
  tpm.t = as.data.frame(t(tpm))
  tpm.t = cbind(tpm.t,group_list)
  tpm.pca <- PCA(tpm.t[,-ncol(tpm.t)],graph = FALSE)
  # fviz_screeplot(tpm.pca, addlabels = TRUE, ylim = c(0, 70),title='Dim choose')
  fviz_pca_ind(tpm.pca,
               mean.point=F,
               # repel = T,
               label = "none",
               geom.ind = c("point",'text'),
               fill.ind = tpm.t$group_list,
               palette = rev(RColorBrewer::brewer.pal(3,"RdBu")),
               legend.title = "Groups",
               pointsize = 6,
               pointshape = 21,
               col.ind = "black",
               title = 'Samples PCA'
  ) + theme(axis.text = element_text(size = 12), axis.title = element_text(size = 15))
}
d_p(pg.pgmc.tpm.clean,group_list_pg)

# Fig.3b -----
library(ggtern)

# pluripotency
for.tern.tpm <- data.frame(pgEpiSC.low = matrixStats::rowMeans2(as.matrix(pg.pgmc.tpm.clean[ 1:3])),
                           pgEpiSC.MPC = matrixStats::rowMeans2(as.matrix(pg.pgmc.tpm.clean[ 4:5])),
                           pgEpiSC.MC = matrixStats::rowMeans2(as.matrix(pg.pgmc.tpm.clean[ 6:8])))

for.tern.tpm$gene <- rownames(pg.pgmc.tpm.clean)
for.tern.tpm$average <- rowMeans(for.tern.tpm[,1:3])

for.tern.tpm$gene.type <- 'Other'
for.tern.tpm[c(na.omit(match(plu.gene[c(1,3,5,8)], for.tern.tpm$gene))), 6] <- 'Pluripotency'

for.tern.tpm$gene.type <- factor(for.tern.tpm$gene.type,levels = c('Pluripotency','Other'))
for.tern.tpm <- for.tern.tpm[order(for.tern.tpm$gene.type),]

ggtern(data=for.tern.tpm,aes(x=pgEpiSC.low,y=pgEpiSC.MPC,z=pgEpiSC.MC))+    
  geom_point_rast(data = for.tern.tpm[for.tern.tpm$gene.type == 'Other',], aes(size = average), color = "grey80",
             alpha = 0.8, show.legend = FALSE) +
  theme_rgbw(base_size = 12 )+   
  labs(title = "Pluripotency genes")+ 
  theme(plot.title = element_text(size=15,hjust = 0.5)) +
  theme_showarrows() + 
  stat_density_tern(aes(fill=..level.., alpha=..level..),geom='polygon',bdl = 0.08, h = 2,alpha=0.8,expand=0.75,bins = 6) +
  scale_fill_gradient2(high = "#4895ef") +
  geom_point(data = for.tern.tpm[for.tern.tpm$gene.type == 'Pluripotency',],
             size = 4,
             color = 'red3',
             show.legend = F)+    
  geom_mask() 


# early myo

for.tern.tpm$gene.type <- 'Other'
for.tern.tpm[c(na.omit(match(c('PAX7','MYOD1','TCF3','SMARCAD1'), for.tern.tpm$gene))), 6] <- 'Early myogenesis'

for.tern.tpm$gene.type <- factor(for.tern.tpm$gene.type,levels = c('Early myogenesis','Other'))
for.tern.tpm <- for.tern.tpm[order(for.tern.tpm$gene.type),]

ggtern(data=for.tern.tpm,aes(x=pgEpiSC.low,y=pgEpiSC.MPC,z=pgEpiSC.MC))+    
  geom_point_rast(data = for.tern.tpm[for.tern.tpm$gene.type == 'Other',], aes(size = average), color = "grey80",
             alpha = 0.8, show.legend = FALSE) +
  theme_rgbw(base_size = 12 )+   
  labs(title = "Early myogenesis genes")+ 
  theme(plot.title = element_text(size=15,hjust = 0.5)) +
  theme_showarrows() + 
  stat_density_tern(aes(fill=..level.., alpha=..level..),geom='polygon',bdl = 0.08, h = 2,alpha=0.8,expand=0.75,bins = 6) +
  scale_fill_gradient2(high = "#4895ef") +
  geom_point(data = for.tern.tpm[for.tern.tpm$gene.type == 'Early myogenesis',],
             size = 3,
             color = 'red3',
             show.legend = F)+    
  geom_mask() 


# late myo

for.tern.tpm$gene.type <- 'Other'
for.tern.tpm[c(na.omit(match(c('MYH2','MYMK','MYH3','MYF5'), for.tern.tpm$gene))), 6] <- 'Late myogenesis'

for.tern.tpm$gene.type <- factor(for.tern.tpm$gene.type,levels = c('Late myogenesis','Other'))
for.tern.tpm <- for.tern.tpm[order(for.tern.tpm$gene.type),]

ggtern(data=for.tern.tpm,aes(x=pgEpiSC.low,y=pgEpiSC.MPC,z=pgEpiSC.MC))+    
  geom_point(data = for.tern.tpm[for.tern.tpm$gene.type == 'Other',], aes(size = average), color = "grey80",
             alpha = 0.8, show.legend = FALSE) +
  theme_rgbw(base_size = 12 )+   
  labs(title = "Late myogenesis genes")+ 
  theme(plot.title = element_text(size=15,hjust = 0.5)) +
  theme_showarrows() + 
  stat_density_tern(aes(fill=..level.., alpha=..level..),geom='polygon',bdl = 0.08, h = 2,alpha=0.8,expand=0.75,bins = 6) +
  scale_fill_gradient2(high = "#4895ef") +
  geom_point(data = for.tern.tpm[for.tern.tpm$gene.type == 'Late myogenesis',], 
             size = 3,
             color = 'red3',
             show.legend = F)
  # ggrepel::geom_label_repel(data = for.tern.tpm[for.tern.tpm$gene.type == 'Late myogenesis',])+
  # geom_mask() +
  # geom_text(data = for.tern.tpm[for.tern.tpm$gene.type == 'Pluripotency',],
  #           aes(label = gene),
  #           check_overlap = F,
  #           color = 'black')

# integrate
for.tern.tpm$gene.type <- 'Other'
# for.tern.tpm[c(na.omit(match(rownames(pick.genes.for.cluster.tpm), for.tern.tpm$gene))), 6] <- c(pick.genes.for.cluster$gene.type[c(1,3:19,21:47)])
for.tern.tpm[c(na.omit(match(rownames(pick.genes.for.cluster.tpm), for.tern.tpm$gene))), 6] <- c(pick.genes.for.cluster$gene.type)

for.tern.tpm$gene.type <- factor(for.tern.tpm$gene.type,levels = c('Pluripotency','Early myogenesis','Late myogenesis','Other'))
for.tern.tpm <- for.tern.tpm[order(for.tern.tpm$gene.type),]

ggtern(data=for.tern.tpm,aes(x=pgEpiSC.low,y=pgEpiSC.MPC,z=pgEpiSC.MC))+    
  geom_point_rast(data = for.tern.tpm[for.tern.tpm$gene.type == 'Other',], aes(size = average), color = "grey90",
             alpha = 0.8, show.legend = FALSE) +
  theme_rgbw(base_size = 12 )+   
  labs(title = "3 cell types genes")+ 
  theme(plot.title = element_text(size=15,hjust = 0.5)) +
  theme_showarrows() + 
  stat_density_tern(aes(fill=..level.., alpha=..level..),geom='polygon',bdl = 0.08, h = 2,alpha=0.8,expand=0.75,bins = 6) +
  scale_fill_gradient2(high = "#a8dadc") +
  geom_point_rast(data = for.tern.tpm[for.tern.tpm$gene.type == 'Pluripotency',], 
             color="blue",
             # fill="#69b3a2",
             shape=16,
             alpha=.8,
             size=4,
             stroke = 1,
             show.legend = F) + 
  geom_point_rast(data = for.tern.tpm[for.tern.tpm$gene.type == 'Early myogenesis',], 
             color="red3",
             # fill="#69b3a2",
             shape=15,
             alpha=.8,
             size=3,
             stroke = 1,
             show.legend = F) + 
  geom_point_rast(data = for.tern.tpm[for.tern.tpm$gene.type == 'Late myogenesis',], 
             color="#2d6a4f",
             # fill="#69b3a2",
             shape=17,
             alpha=.8,
             size=3,
             stroke = 1,
             show.legend = F)

# Fig.3c -----
tpm.t <- data.frame(t(pg.pgmc.tpm.clean))
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

prcomp.score[c(na.omit(match(unique(c(rownames(diff.genes.tpm), 'PAX7','MUSK','TNNT1','MYOD1','PITX3','MYH7','MYH11','MEF2A','MEF2C','MYMK','MYOG','DES','MYH3','MYH8','RUNX1')), prcomp.score$Symbol))), 4] <- unique(c(rownames(diff.genes.tpm), 'PAX7','MUSK','TNNT1','MYOD1','PITX3','MYH7','MYH11','MEF2A','MEF2C','MYMK','MYOG','DES','MYH3','MYH8','RUNX1'))
p1 <- ggplot(prcomp.score, aes(x=Rank, y=PC1)) +
  geom_point_rast() +
  geom_text_repel(
    aes(label = highlight),
    size = 1,
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

# Fig.3d & Fig.S7c & Fig.S7d ------
pg.pgmc.tpm.clean.median <- data.frame(pg.median = matrixStats::rowMedians(as.matrix(pg.pgmc.tpm.clean),cols = 1:3),
                                       pg.mpc.median = matrixStats::rowMedians(as.matrix(pg.pgmc.tpm.clean),cols = 4:5),
                                       pg.mc.median = matrixStats::rowMedians(as.matrix(pg.pgmc.tpm.clean),cols = 6:8))
rownames(pg.pgmc.tpm.clean.median) <- rownames(pg.pgmc.tpm.clean)
pg.pgmc.tpm.clean.median.scale <- data.frame(t(scale(t(pg.pgmc.tpm.clean.median))))

nk=2:100
Wss<-sapply(nk,function(k){
  kmeans(pg.pgmc.tpm.clean.median.scale,centers = k,iter.max = 100)$tot.withinss})
plot(nk,Wss,type = "o",xlab="Number of k",ylab="Within sum of squares",col='red3',)
abline(v=10,col='black')
legend("topleft", legend = 'k = 10', 
       col= "red3",
       pch = 15, bty = "n", pt.cex = 2, cex = 1.2,  horiz = F, inset =  0.1)



fit <- kmeans(x = pg.pgmc.tpm.clean.median.scale[,c(1:3)],10,iter.max = 100)
kmeansBIC(fit)
table(fit$cluster)
pg.pgmc.tpm.clean.median.scale$kmeans.cluster <- fit$cluster
pg.pgmc.tpm.clean.median.scale$gene <- rownames(pg.pgmc.tpm.clean.median)

kmectocluster <- data.frame(cluster = fit$cluster)

pg.pgmc.tpm.clean.median.log2 <- log2(pg.pgmc.tpm.clean.median+1)
pg.pgmc.tpm.clean.median.log2$gene <- rownames(pg.pgmc.tpm.clean.median.log2)

kmectocluster <- 
  kmectocluster %>% 
  rownames_to_column(var = 'gene') %>%
  left_join(y = pg.pgmc.tpm.clean.median.log2,
            by = 'gene')
kmectocluster$cluster <- paste0('C_',kmectocluster$cluster)
kmectocluster$cluster <- factor(kmectocluster$cluster,levels = c(paste0('C_',seq(1,10))))
kmectocluster <- kmectocluster %>% group_by(cluster) %>%  mutate(cluster=paste0(cluster,"(",n()," genes)"))
table(kmectocluster$cluster)
kmectocluster$cluster <- factor(kmectocluster$cluster,levels = c('C_1(1232 genes)', 'C_2(1523 genes)', 'C_3(1253 genes)', 'C_4(1333 genes)',  'C_5(1349 genes)', 'C_6(644 genes)',  'C_7(1343 genes)','C_8(1830 genes)',  'C_9(1373 genes)','C_10(2666 genes)'))

kmectocluster %>% 
  select('pg.median','pg.mpc.median','pg.mc.median','cluster') %>%
  gather(key = 'stage',value = 'exp', -cluster) %>%
  mutate(stage=fct_relevel(stage, 'pg.median','pg.mpc.median','pg.mc.median')) %>%
  ggplot( aes(x=stage, y=round(exp,4),group=1,fill=cluster)) +
  stat_summary(fun.data="mean_cl_boot", geom="ribbon",
               #width=.2, 
               alpha=I(.2)) +
  stat_summary(fun="mean", geom="line") +
  labs(x="Stage(cluster base on scale data)", y="Expression level(log2(TPM+1))") +
  theme_bw() +
  theme(legend.position="none") +
  scale_fill_manual(values = rep('red',50),
                    guide=guide_legend(direction="vertical",
                                       label.position="right",
                                       title=NULL,
                                       ncol=6,
                                       label.hjust=0.8))+
  scale_color_manual(values =  rep('red',50),guide = 'none')+
  # geom_smooth()+
  # ylim(0,1)+
  facet_wrap(~cluster,scales = 'free_y',ncol = 5,drop = F)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(labels = c('PgEpiSC','PgEpiSC.MPC','PgEpiSC.MC'))

library(rio)
export(kmectocluster,file = '../../MSC/kmeans/kmeans.cluster.xlsx')

# kmectocluster <- import('../../MSC/kmeans/kmeans.cluster.xlsx')
pg.pgmc.tpm.clean.cluster <- pg.pgmc.tpm.clean %>% rownames_to_column(var = 'gene') %>% left_join(y = kmectocluster,by = 'gene')

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

pg.pgmc.tpm.clean.cluster <- pg.pgmc.tpm.clean.cluster %>% left_join(y = orth[,c(1:4)],by = c('gene'='genename'))
export(pg.pgmc.tpm.clean.cluster,file = '../../MSC/kmeans/pg.pgmc.tpm.clean.cluster.xlsx')

rownames(pg.pgmc.tpm.clean.cluster) <- pg.pgmc.tpm.clean.cluster$gene
pick.genes.for.cluster <- read.table('../../MSC/kmeans/pick.genes.for.cluster.txt',header = F,sep = '\t')
pick.genes.for.cluster$gene.type <- c(rep('Pluripotency',9), rep('Early myogenesis',8), rep('Late myogenesis',18), rep('Collagen',15))
pick.genes.for.cluster.tpm <- pg.pgmc.tpm.clean.cluster %>% filter(gene %in% pick.genes.for.cluster$V1)

pick.genes.for.cluster.tpm <- separate(pick.genes.for.cluster.tpm, col = cluster, into = c('clusters','genenum'), remove = T)
pick.genes.for.cluster.tpm$cluster <- paste0(pick.genes.for.cluster.tpm$clusters,"_",pick.genes.for.cluster.tpm$genenum)
pick.genes.for.cluster.tpm <- pick.genes.for.cluster.tpm %>% left_join( y = pick.genes.for.cluster, by = c('gene' = 'V1'))
rownames(pick.genes.for.cluster.tpm) <- pick.genes.for.cluster.tpm$gene
pick.genes.for.cluster.tpm$cluster <- factor(pick.genes.for.cluster.tpm$cluster, levels = c(paste0("C_", c(1,2,3,5,7,8,10))))
# pick.genes.for.cluster.tpm <- na.omit(pick.genes.for.cluster.tpm[pick.genes.for.cluster$V1,])
pick.genes.for.cluster.tpm$gene.type <- factor(pick.genes.for.cluster.tpm$gene.type, levels = c('Pluripotency','Early myogenesis','Late myogenesis','Collagen'))

pick.genes.for.cluster.tpm <- pick.genes.for.cluster.tpm[pick.genes.for.cluster$V1,]

ha_left <- list(cluster=pal_npg("nrc")(10)[c(4,2:3,5,7,8,10)], gene.type = npg.col[13:16],Group = col)
names(ha_left$cluster) = levels(pick.genes.for.cluster.tpm$cluster)
names(ha_left$gene.type) = levels(pick.genes.for.cluster.tpm$gene.type)
names(ha_left$Group) = c('pgEpiSC','pgEpiSC.MPC','pgEpiSC.MC')

ComplexHeatmap::pheatmap(t(scale(t(pick.genes.for.cluster.tpm[,c(2:9)]))),
                         border_color = NA,
                         clustering_method = 'ward.D',
                         show_rownames = T,
                         cluster_cols = F,cluster_rows = F,
                         color = scales::alpha(rev(colorRampPalette(RColorBrewer::brewer.pal(5,"RdBu"),alpha=T,bias=.8)(256)),alpha = .8),
                         angle_col = '315',
                         annotation_row = pick.genes.for.cluster.tpm[,c(18,19)],
                         annotation_colors = ha_left,
                         annotation_col = annotation_col_matrix,
                         breaks = unique(c(seq(-2,2, length=256))),fontsize_row = 5
                         )


