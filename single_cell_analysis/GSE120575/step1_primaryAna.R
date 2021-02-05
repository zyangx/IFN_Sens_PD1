library(org.Hs.eg.db)

inputData.df <- read.csv("GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt", header=T, row.names=1, sep="\t", check.names=F)
colname.v <- as.vector(colnames(inputData.df))[-1]
inputData.df <- inputData.df[,-16292]
colnames(inputData.df) <- colname.v

data.m <- as.matrix(inputData.df)

metadata.df <- read.csv("metadata.txt", header=T,  sep="\t", check.names=F)
rownames(metadata.df) <- metadata.df$title
metadata.df <- metadata.df[,5:7]
 
sce <- CreateSeuratObject(counts = data.m, assay = "RNA", meta.data=metadata.df, min.cells=0, min.features=0,  project = "Melanoma") 

sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")
#sce <- NormalizeData(sce)
#sce <- ScaleData(sce), display.progress = F) 

sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000)


plot1 <- VariableFeaturePlot(sce)
top10 <- head(VariableFeatures(sce), 10)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2



all.genes <- rownames(sce)
sce <- ScaleData(sce, features = all.genes)


### 
#
sce <- RunPCA(sce, features = VariableFeatures(object = sce))
#
print(sce[["pca"]], dims = 1:5, nfeatures = 5)
#
VizDimLoadings(sce, dims = 1:2, reduction = "pca")
graph2ppt(file="4.PC1nPC2_value.pptx", width=9, aspectr=1.5)
#
DimPlot(sce, reduction = "pca")
graph2ppt(file="5.PCA降维.pptx", width=9, aspectr=1.5)
#
DimHeatmap(sce, dims = 1, cells = 500, balanced = TRUE)
graph2ppt(file="6.单个PC维度_heatmap.pptx", width=9, aspectr=1.5)

DimHeatmap(sce, dims = 1:9, cells = 500, balanced = TRUE)
graph2ppt(file="7-1.1到9个PC维度_heatmap.pptx", width=12, aspectr=1.5)

DimHeatmap(sce, dims = 10:18, cells = 500, balanced = TRUE)
graph2ppt(file="7-2.9到18个PC维度_heatmap.pptx", width=12, aspectr=1.5)

DimHeatmap(sce, dims = 19:27, cells = 500, balanced = TRUE)
graph2ppt(file="7-3.18到27个PC维度_heatmap.pptx", width=12, aspectr=1.5)

###
saveRDS(sce, file = "./PCA.rds")



### 
#
sce <- JackStraw(sce, num.replicate = 100)
sce <- ScoreJackStraw(sce, dims = 1:20)#
JackStrawPlot(sce, dims = 1:20) #
graph2ppt(file="8.PCA维度_JackStraw.pptx", width=12, aspectr=1.5)
#
#
ElbowPlot(sce)
graph2ppt(file="9.PCA维度_ElbowPlot.pptx", width=12, aspectr=1.5)
#


### 
#
sce <- FindNeighbors(sce, dims = 1:20)
#
sce <- FindClusters(sce, resolution = 0.05)
##
sce <- FindClusters(sce, resolution = c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))
##
sce@meta.data %>% View()
#install.packages('clustree')
library(clustree)
clustree(sce@meta.data, prefix = "RNA_snn_res.")
graph2ppt(file="10.clustree.pptx", width=9, aspectr=1.5)

#
sce <- FindClusters(sce, resolution = 0.3)



### 
#
sce <- RunUMAP(sce, dims = 1:20)
DimPlot(sce, reduction = "umap")
graph2ppt(file="11.umap.pptx", width=9, aspectr=1.5)
#
sce <- RunTSNE(sce, dims = 1:20)
pdf("12.tsne.pdf",width=9,height=9)
DimPlot(sce, reduction = "tsne")
dev.off()

graph2ppt(file="12.tsne.pptx", width=9, aspectr=1.5)

###
saveRDS(sce, file = "./ForDoublets.rds")

cluster1.markers <- FindMarkers(sce, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)

pdf("All_cells_CD8_vlnplot.pdf", width=5.5, height=5)
VlnPlot(sce, features = c("CD8A"))
dev.off()

pdf("All_cells_CD8_featureplot.pdf", width=5.5, height=5)
FeaturePlot(sce, features = c("CD8A"))
dev.off()


cluster1.markers <- FindMarkers(sce, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)

sce.markers <- FindAllMarkers(sce, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

top10 <- sce.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

pdf("Top10Gens_heatmap_inCluster.pdf", width=12, height=14)
DoHeatmap(sce, features = top10$gene) 
dev.off()


####################################
library(SingleR)
#ref <- BlueprintEncodeData()
load("SingleRBlueprintEncode.Rdata")

#ref <- HumanPrimaryCellAtlasData()
#load("HumanPrimaryCellAtlasData.Rdata")

pred.data <- SingleR(test=data.m, ref=ref, labels=ref$label.main)
table(pred.data$labels)
sum(is.na(pred.data$pruned.labels))
tab.data <- table(Assigned=pred.data$pruned.labels, Cluster=Idents(sce))


# Adding a pseudo-count of 10 to avoid strong color jumps with just 1 cell.
library(pheatmap)
pdf("All_cell_types_assign_by_singleR.pdf", width=5, height=5)
pheatmap(log2(tab.data+10), color=colorRampPalette(c("white", "blue"))(101))
dev.off()

#AddMetaData(object, metadata, col.name = NULL)
######################################

sce1 <- SubsetData(sce, ident.use= "0")

sce1 <- FindVariableFeatures(sce1, selection.method = "vst", nfeatures = 2000)

sce1 <- RunPCA(sce1, features = VariableFeatures(object = sce1))

DimPlot(sce1, reduction = "pca")

pdf("CD8cell_PC1_DimHeatmap.pdf", width=7, height=5)
DimHeatmap(sce1, dims = 1, cells = 500, balanced = TRUE)
dev.off()

pdf("CD8cell_PC1-9_DimHeatmap.pdf", width=10, height=10)
DimHeatmap(sce1, dims = 1:9, cells = 500, balanced = TRUE)
dev.off()

sce1 <- FindClusters(sce1, resolution = 0.3)

sce1 <- RunUMAP(sce1, dims = 1:20)
DimPlot(sce1, reduction = "umap")

VlnPlot(sce1, features = c("CD27"))

pdf("CD8cells_feature_by_CD27_tsne.pdf", width=5.5, height=5)
FeaturePlot(sce1, features = c("CD27"), reduction="tsne")
dev.off()

pdf("CD8cells_feature_by_CD27_umap.pdf", width=5.5, height=5)
FeaturePlot(sce1, features = c("CD27"))
dev.off()







