files <- list.files("fpkm_reps",pattern="*tab",full.names=T)
fp <- map(files,read.table) %>% reduce(full_join,by="V1")
colnames(fp) <- c("gene",gsub(".fpkm.uniq.tab","",files))
colnames(fp) <- gsub("^.*RNA_","",colnames(fp))
rownames(fp) <- fp$gene
fp_log <- apply(fp[,-1],2,myfun_log2)
pc1 <- prcomp(t(fp_log))
pca <- as.data.frame(pc1$x[,1:2])
samples <- rownames(pca)
pca <- cbind(pca,condition=samples)
pca <- as.data.frame(pca)
pca$condition <- gsub("_[1234]$","",pca$condition)

pca.var <- pc1$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100,1)
require("ggrepel")
pointLabel <- pca$condition
ggplot(pca, aes(PC1, PC2,color=condition)) +
  geom_point(size=3) +theme_bw()+
  xlab(paste0("PC1: ",pca.var.per[1],"% variance")) +
  ylab(paste0("PC2: ",pca.var.per[2],"% variance")) + 
  #scale_color_manual(values = c("darkred","blue","green")) +
  coord_fixed()+
  geom_text_repel(aes(label = pointLabel),size = 3)
dev.copy2pdf(file="pca_fpkm.D0.xEnd.RNA.pdf")

dist_mat <- dist(t(as_tibble(fp[,-1])))
hc_cluster <- hclust(dist_mat,method = "complete")

plot(hc_cluster,cex=0.8,labels=colnames(fp)[-1])
dev.copy2pdf(file="cluster.E4.5_to_E7.5.pdf")



#####eps_es_em
fp <- myfun_f_rbind("fpkm_reps/eps_es_em",".fpkm.uniq.tab")
colnames(fp) <- gsub("RNA.","",colnames(fp))
colnames(fp) <- gsub("-","_",colnames(fp))
fp %<>% dplyr::select(gene,DHK.mES_1,DHK.mES_2,DHK.EPS_1:RC_LCDM_2) %>% select(-ICM_1:-ICM_4,-Morula_1:-Morula_2)
fp_log <- apply(fp[,-1],2,myfun_log2)
pc1 <- prcomp(t(fp_log))
pca <- as.data.frame(pc1$x[,1:3])
samples <- rownames(pca)
pca <- cbind(pca,condition=samples)
pca <- as.data.frame(pca)
pca$condition <- gsub("_[1234]$","",pca$condition)

pca.var <- pc1$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100,1)
require("ggrepel")
pointLabel <- pca$condition
pc12 <- ggplot(pca, aes(PC1, PC2,color=condition)) +
  geom_point(size=3) +theme_bw()+
  xlab(paste0("PC1: ",pca.var.per[1],"% variance")) +
  ylab(paste0("PC2: ",pca.var.per[2],"% variance")) + 
  #scale_color_manual(values = c("darkred","blue","green")) +
  coord_fixed()+
  geom_text_repel(aes(label = pointLabel),size = 3,max.overlaps=20)+theme_classic(12)+theme(legend.position="none")

pc13 <- ggplot(pca, aes(PC1, PC3,color=condition)) +
  geom_point(size=3) +theme_bw()+
  xlab(paste0("PC1: ",pca.var.per[1],"% variance")) +
  ylab(paste0("PC3: ",pca.var.per[3],"% variance")) + 
  #scale_color_manual(values = c("darkred","blue","green")) +
  coord_fixed()+
  geom_text_repel(aes(label = pointLabel),size = 3,max.overlaps=20)+theme_classic(12)+theme(legend.position="none")


pc23 <- ggplot(pca, aes(PC2, PC3,color=condition)) +
  geom_point(size=3) +theme_bw()+
  xlab(paste0("PC2: ",pca.var.per[2],"% variance")) +
  ylab(paste0("PC3: ",pca.var.per[3],"% variance")) + 
  #scale_color_manual(values = c("darkred","blue","green")) +
  coord_fixed()+
  geom_text_repel(aes(label = pointLabel),size = 3)+theme_classic(12)+theme(legend.position="none")
plot_grid(pc12,pc13,pc23, 
          labels = LETTERS[1:3],  #大写字母取6个 当做图标
          ncol = 2)

ggsave(file="pca_fpkm.eps.es.pdf")

        
csif <- data.frame(group_name=colnames(fp_log),batch=c(rep(1,4),rep(2,ncol(fp_log)-4)))
csif$cell <- c(rep("ES",2),rep("EPS",6),rep("ES",4),rep("EPS",2))
fp_log_f <- as.data.frame(fp_log) %>% filter(rowSums(fp_log)!=0,rowVars(fp_log)!=0)
fp_log_f <- as.matrix(fp_log_f)
modcombat = model.matrix(~cell, data = csif)
batch = csif$batch
combat_edata = ComBat(dat=fp_log_f, batch=batch, mod=modcombat)
pc1 <- prcomp(t(combat_edata))
pca <- as.data.frame(pc1$x[,1:3])
# pca$cell <- csif$cell
samples <- rownames(pca)
pca <- cbind(pca,condition=samples)
pca <- as.data.frame(pca)
pca$condition <- gsub("_[1234]$","",pca$condition)
pca.var <- pc1$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100,1)
require("ggrepel")
pointLabel <- pca$condition
pc12 <- ggplot(pca, aes(PC1, PC2,color=condition)) +
  geom_point(size=3) +theme_bw()+
  xlab(paste0("PC1: ",pca.var.per[1],"% variance")) +
  ylab(paste0("PC2: ",pca.var.per[2],"% variance")) + 
  #scale_color_manual(values = c("darkred","blue","green")) +
  coord_fixed()+
  geom_text_repel(aes(label = pointLabel),size = 3,max.overlaps=20)+theme_classic(12)+theme(legend.position="none")
pc13 <- ggplot(pca, aes(PC1, PC3,color=condition)) +
  geom_point(size=3) +theme_bw()+
  xlab(paste0("PC1: ",pca.var.per[1],"% variance")) +
  ylab(paste0("PC3: ",pca.var.per[3],"% variance")) + 
  #scale_color_manual(values = c("darkred","blue","green")) +
  coord_fixed()+
  geom_text_repel(aes(label = pointLabel),size = 3,max.overlaps=20)+theme_classic(12)+theme(legend.position="none")
pc23 <- ggplot(pca, aes(PC2, PC3,color=condition)) +
  geom_point(size=3) +theme_bw()+
  xlab(paste0("PC2: ",pca.var.per[2],"% variance")) +
  ylab(paste0("PC3: ",pca.var.per[3],"% variance")) + 
  #scale_color_manual(values = c("darkred","blue","green")) +
  coord_fixed()+
  geom_text_repel(aes(label = pointLabel),size = 3)+theme_classic(12)+theme(legend.position="none")
plot_grid(pc12,pc13,pc23, 
          labels = LETTERS[1:3],  #大写字母取6个 当做图标
          ncol = 2)
ggsave(file="pca_fpkm.ComBat.eps.es.pdf")
##3D PCA
library("scatterplot3d")
color <- c(rep("#E6550DFF",2),rep("#FD8D3CFF",2),rep("#3182BDFF",2),
  rep("#6BAED6FF",2),rep("#756BB1FF",2),rep("#9E9AC8FF",2),rep("#BCBDDCFF",2))
 pca$col <- color
pca$pch <- c(rep(17,2),rep(19,6),rep(17,4),rep(19,2))
scatterplot3d(pca[,c(2,3,1)],main='PCA fpkm ComBat',type='h',color=pca$col,
              highlight.3d=F,angle=60,grid=T,box=F,scale.y=1,
              cex.symbols=1.2,pch=pca$pch,
              col.grid='lightblue',lty.hplot="dashed",
               xlab=paste0("PC2: ",pca.var.per[2],"% variance"),
              ylab=paste0("PC3: ",pca.var.per[3],"% variance"),
              zlab=paste0("PC1: ",pca.var.per[1],"% variance"))
legend("topleft",pca$condition[c(1,9,11,3,5,7,13)],box.col="white",
  pch=pca$pch[c(1,9,11,3,5,7,13)],col=pca$col[c(1,9,11,3,5,7,13)])
dev.copy2pdf(file="pca_fpkm.ComBat.eps.es.3D.pdf")
