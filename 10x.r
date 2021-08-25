library(Seurat)
library(gridExtra)
library(ggplot2)
library(multtest)
library(metap)
d1 <- Read10X("/home1/xuxc/Blastoid/10x/eps_blastoid/outs/filtered_feature_bc_matrix")
d2 <- Read10X_h5("/home1/xuxc/Blastoid/10x/MAESTRO/MAESTRO_gene_count.smartseq2.blastoid.E3.5.E4.5.h5")
d1.obj <- CreateSeuratObject(counts = d1, project = "10x")
d2.obj <- CreateSeuratObject(counts = d2, project = "smart")
combined <- merge(d1.obj, y = d2.obj, add.cell.ids = c("10x", "smart"), project = "combine")
ifnb.list <- SplitObject(combined,split.by = "ident")
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})


features <- SelectIntegrationFeatures(object.list = ifnb.list)
immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features) ###
# this command creates an 'integrated' data assay
immune.combined <- IntegrateData(anchorset = immune.anchors)

# Perform an integrated analysis

# Now we can run a single integrated analysis on all cells!

# specify that we will perform downstream analysis on the corrected data note that the original
# unmodified data still resides in the 'RNA' assay
DefaultAssay(immune.combined) <- "integrated"
# Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindClusters(immune.combined, resolution = 0.15)

pdf("cluster.raw.pdf")
DimPlot(immune.combined, reduction = "umap", label = TRUE, repel = TRUE)
dev.off()
pdf("integrate.split.pdf")
DimPlot(immune.combined, reduction = "umap", split.by = "orig.ident",label = TRUE, repel = TRUE)
dev.off()


DefaultAssay(immune.combined) <- "RNA"
##rename cluster
immune.combined <- RenameIdents(immune.combined, `0` = "int-2", `1` = "int-1", `2` = "PrE", 
    `3` = "ICM", `4` = "TE",`5` = "ICM")
DimPlot(immune.combined, label = TRUE,split.by = "orig.ident")
ggsave("cluster.name.pdf",width=10,height=5)


##marker gene
PrE <- c("Pdgfra", "Gata6", "Dab2","Sox17" )
TE <- c("Cdx2", "Krt8", "Krt18","Ascl2","Tacstd2","Tfap2c")
ICM <- c("Nanog","Sox2", "Pou5f1", "Sox15",  "Esrrb","Tdgf1" )
selected <- c("Cdx2","Gata6","Sox17")
# VE <- c("Foxa2","Ihh","Cyp26a1","Furin","Afp","Cited1") 
# PE <- c("Plat","Pth1r","Fst","Sparc","Krt19","Lamb1")
# ECT <- c("Nes","Sox1","Fgf5","Gfap","Rest")
# marker <- list(PrE=PrE,TE=TE,ICM=ICM,VE=VE,PE=PE,ECT=ECT)
marker <- list(PrE=PrE,TE=TE,ICM=ICM,selected=selected)

plots <- lapply(X = marker, FUN = function(x) {
FeaturePlot(immune.combined,features = x,split.by = "orig.ident",keep.scale ="all",
	min.cutoff =0,max.cutoff =2,pt.size=0.6,label=TRUE)+theme(legend.position = c(1,0.1)) & 
scale_colour_gradientn(colours=c("grey","#FFBF00FF","#FF8000FF","#FF4000FF","#FF2000FF") )
})
name <- names(plots)
width <- rep(7,4)
height <- c(12,18,18,9)
l <- list(a = name, b = width, c = height)
pmap(l, function(a, b, c) ggsave(filename=paste("FeaturePlot.",a,".pdf",sep=""), plot=plots[[a]],width=b,height=c))



##pie plot
cluster <- as.data.frame(immune.combined@active.ident)
cluster$name <- rownames(cluster)
cluster <- cluster %>% separate(name,c("platform","id"),sep="_")
 colnames(cluster)[1] <- "cl"
pie_data <- cluster %>% dplyr::filter(platform=="10x") %>% group_by(cl) %>% summarize(n=n()) %>% mutate(p=n/sum(n)) 
pdf(file = "blastiod.pie.pdf")
pie(pie_data$p,paste(pie_data$cl,"(",round(100*pie_data$p,2),"%)",sep=""),col=hue_pal()(5))
dev.off()



##add metadata, split by tissue source
CellsMeta = immune.combined@meta.data
source <- c(rep("blastoid",3603),rep("blastoyst",277))
CellsMeta["source"] <- source 
CellsMetaTrim <- subset(CellsMeta, select = c("source"))
immune.combined <- AddMetaData(immune.combined, CellsMetaTrim)
DimPlot(immune.combined, label = TRUE,split.by = "source")
ggsave("cluster.split.by.source.pdf",width=10,height=5)

plots <- lapply(X = marker, FUN = function(x) {
FeaturePlot(immune.combined,features = x,split.by = "source",keep.scale ="all",
	min.cutoff =0,max.cutoff =2,pt.size=0.6,label=TRUE)+theme(legend.position = c(1,0.1)) & 
scale_colour_gradientn(colours=c("grey","#FFBF00FF","#FF8000FF","#FF4000FF","#FF2000FF") )
})
name <- names(plots)
width <- rep(7,4)
height <- c(12,18,18,9)
l <- list(a = name, b = width, c = height)
pmap(l, function(a, b, c) ggsave(filename=paste("FeaturePlot.",a,".V2.pdf",sep=""), plot=plots[[a]],width=b,height=c))



##highlight target cell
umap <- immune.combined@reductions$umap@cell.embeddings %>% as.data.frame()
umap.cl <- cbind(umap,cluster)
cell.id <- paste("c",1:4,sep="")
umap.cl.cource <- cbind(umap.cl,CellsMetaTrim)
umap.cl.cource <- umap.cl.cource %>% mutate(cl.source=paste(cl,source,sep="."))

ggplot(umap.cl, aes(x=UMAP_1, y=UMAP_2, color=cl)) + geom_point(size=0.8)+
labs(title="cdx2 PrE, oct4 int2")+theme_classic()+
geom_point(data=umap.cl[umap.cl$id %in% cell.id,],aes(x=UMAP_1, y=UMAP_2),col="black",size=2)
ggsave("cluster.highlight.pdf",width=7,height=5) #cdx2 PrE, oct4 int2
# ggplot(umap.cl.cource, aes(x=UMAP_1, y=UMAP_2, color=cl)) + geom_point(aes(shape=source),size=0.8)+
# labs(title="cdx2 PrE, oct4 int2")+theme_classic()+scale_shape_manual(values=c(1,4))+
# geom_point(data=umap.cl[umap.cl$id %in% cell.id,],aes(x=UMAP_1, y=UMAP_2),col="black",size=2)
# ggsave("cluster.highlight.V2.pdf",width=7,height=5) #cdx2 PrE, oct4 int2
ggplot(umap.cl.cource, aes(x=UMAP_1, y=UMAP_2, color=cl.source)) + geom_point(size=0.8)+
labs(title="cdx2 PrE, oct4 int2")+theme_classic()+scale_color_brewer(palette = "Paired")+
geom_point(data=umap.cl[umap.cl$id %in% cell.id,],aes(x=UMAP_1, y=UMAP_2),col="black",size=2)
ggsave("cluster.highlight.V3.pdf",width=7,height=5) #cdx2 PrE, oct4 int2
ggplot(umap.cl.cource, aes(x=UMAP_1, y=UMAP_2, color=source)) + geom_point(size=0.8)+
labs(title="cdx2 PrE, oct4 int2")+theme_classic()+
geom_point(data=umap.cl[umap.cl$id %in% cell.id,],aes(x=UMAP_1, y=UMAP_2),col="black",size=2)
ggsave("cluster.highlight.V4.pdf",width=7,height=5) #cdx2 PrE, oct4 int2
# ggplot(umap.cl.cource, aes(x=UMAP_1, y=UMAP_2, color=cl)) + geom_point(aes(alpha=source),size=0.8)+
# labs(title="cdx2 PrE, oct4 int2")+theme_classic()+scale_alpha_manual(values=c(0.2,1))+
# geom_point(data=umap.cl[umap.cl$id %in% cell.id,],aes(x=UMAP_1, y=UMAP_2),col="black",size=2)
# ggsave("cluster.highlight.V5.pdf",width=7,height=5) #cdx2 PrE, oct4 int2
