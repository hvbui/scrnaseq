library(Seurat)
library(readr)
library(dplyr)
library(patchwork)
library(data.table)
library(ggplot2)
library(cowplot)
library(pheatmap)
library(sctransform)
library(SeuratDisk)
library(scater)
library(loomR)
library(devtools)
library(funcutils)
library(rlist)
library(tidyverse)
library(SeuratWrappers)
library(gplots)
#read in data
wat.data = Read10X(data.dir = "//research.drive.wisc.edu/galmozzi/Galmozzi lab/Hoang/Projects/Adipocyte Subtypes/new_analysis/0WAT/filtered_feature_bc_matrix")
wat.data = Read10X(data.dir = "F:/new_analysis/0WAT/outs/filtered_feature_bc_matrix_CAP")
#convert into seurat object
wat = CreateSeuratObject(counts = wat.data, project = "0watpreads", min.cells = 3, min.features = 200)
wat$data.set = rep("wat", length(wat$orig.ident))

#calculate percent mitochondrial RNA for watpreads
wat[["percent.mt"]] = PercentageFeatureSet(wat, pattern = "^mt-")
wat[["percent.mt"]] = PercentageFeatureSet(wat, pattern = "^MT-")
#generate violin plot #genes/cell, mtRNA/cell, reads/cell
VlnPlot(wat, features= c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size=0.5, ncol = 3)
#rerun without the dots
VlnPlot(wat, features= c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size=0, ncol = 3)

#Feature Scatter visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(wat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(wat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#apply filter-- UMI 6000-60000, features 3900-7700, mito umi 1.7-5
wat <- subset(wat, subset = nFeature_RNA > 3900 & nFeature_RNA < 7700 & percent.mt < 5& percent.mt > 1.7 & nCount_RNA > 6000 & nCount_RNA < 60000)

#normalizing the data -normalizes the feature expression measurements for each cell by the total expression, multiplies this by a scale factor (10,000 by default), and log-transforms the result
wat <- NormalizeData(wat, normalization.method = "LogNormalize", scale.factor = 10000)

#identification of highly variable features
wat <- FindVariableFeatures(wat, selection.method = "vst", nfeatures = 2000)

#Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(wat), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(wat)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#scaling the data -pre-processing step prior to dimensional reduction techniques like PCA. The ScaleData() function:
#Shifts the expression of each gene, so that the mean expression across cells is 0
#Scales the expression of each gene, so that the variance across cells is 1
#This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
#The results of this are stored in pbmc[["RNA"]]@scale.data 
all.genes <- rownames(wat)
wat <- ScaleData(wat, features = all.genes)

#PCA
wat <- RunPCA(wat, features = VariableFeatures(object = wat))

#determine dimensionality
ElbowPlot(wat, ndims = 50, reduction = "pca")
wat <- JackStraw(wat, num.replicate=100)
wat <- ScoreJackStraw(wat, dims = 1:20)
JackStrawPlot(wat, dims = 1:15)

#clustering
wat <- FindNeighbors(wat, dims = 1:10)
wat <- FindClusters(wat, resolution = 0.5)

#umap
wat <- RunUMAP(wat, dims = 1:10)
DimPlot(wat, reduction = "umap")

#save object so it doesnt need to rerun all previous codes
#saveRDS(wat, file ="//research.drive.wisc.edu/galmozzi/Galmozzi lab/Hoang/Projects/Bioinformatics/new_analysis/0wat/R/23-1-10.rds")

# find all markers of cluster 2
cluster2.markers <- FindMarkers(wat, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)

#find all markers
wat.markers<- FindAllMarkers(wat)
write.csv(wat.markers, "F:/wat.markers.csv")

#remove cluster 6 7 and recluster
wat.exclude <- subset(x=wat, idents = c(6,7), invert = TRUE)
wat.exclude <- NormalizeData(wat.exclude, normalization.method = "LogNormalize", scale.factor = 10000)
wat.exclude <- FindVariableFeatures(wat.exclude, selection.method = "vst")
wat.exclude <- RunPCA(wat.exclude, features = VariableFeatures(object = wat.exclude))
wat.exclude <- FindNeighbors(wat.exclude, dims = 1:10)
wat.exclude <- FindClusters(wat.exclude, resolution = 0.5)
wat.exclude <- RunUMAP(wat.exclude, dims = 1:10)
DimPlot(wat.exclude, reduction = "umap")

#Find all markers from new cells - removing cluster 7 8 9
cluster0.markers <- FindMarkers(wat.exclude, ident.1 = 0, min.pct = 0.25)
cluster1.markers <- FindMarkers(wat.exclude, ident.1 = 1, min.pct = 0.25)
cluster2.markers <- FindMarkers(wat.exclude, ident.1 = 2, min.pct = 0.25)
cluster3.markers <- FindMarkers(wat.exclude, ident.1 = 3, min.pct = 0.25)
cluster4.markers <- FindMarkers(wat.exclude, ident.1 = 4, min.pct = 0.25)

#export cluster genes
write.csv(cluster0.markers, "F:/cluster0.csv")
write.csv(cluster1.markers, "F:/cluster1.csv")
write.csv(cluster2.markers, "F:/cluster2.csv")
write.csv(cluster3.markers, "F:/cluster3.csv")
write.csv(cluster4.markers, "F:/cluster4.csv")
write.csv(cluster5.markers, "F:/cluster5.csv")

#export uMAP projections
umap = cbind("Barcode" = rownames(Embeddings(object = wat.exclude, reduction = "umap")), Embeddings(object = wat.exclude, reduction = "umap"))
write.table(umap, file="F:/uMAP_seurat.csv", sep = ",", quote = F, row.names = F, col.names = T)

#export metadata to get clusters export
write.csv(wat.exclude@meta.data,"F:/seurat_metadata.csv")

#CELL CYCLE SCORING
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
wat.exclude <- CellCycleScoring(wat.exclude, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

#for uncapitalized genes
m.s.genes <- str_to_title(cc.genes$s.genes)
m.g2m.genes <- str_to_title(cc.genes$g2m.genes)
wat.exclude <- CellCycleScoring(wat.exclude, s.features = m.s.genes, g2m.features = m.g2m.genes, set.ident = TRUE)

# view cell cycle scores and phase assignments
head(wat.exclude[[]])
DimPlot(wat.exclude, reduction="umap", group.by="Phase")

#run PCA on cell cycle 
bat <- RunPCA(bat, features = c(s.genes, g2m.genes))
DimPlot(bat)

#ALTERNATE WORKFLOW- CELL CYCLE DIFFERENCE REGRESSION
wat.exclude$CC.Difference <- wat.exclude$S.Score - wat.exclude$G2M.Score
wat.exclude_ccdiff <- ScaleData(wat.exclude, vars.to.regress = "CC.Difference", features = rownames(wat.exclude))
# cell cycle effects strongly mitigated in PCA
wat.exclude_cc <- RunPCA(exclude, features = VariableFeatures(exclude), nfeatures.print = 10)
# when running a PCA on cell cycle genes, actively proliferating cells remain distinct from G1
# cells however, within actively proliferating cells, G2M and S phase cells group together
#group.by recognizes all metadata column names
wat.exclude <- RunPCA(wat.exclude, features = c(s.genes, g2m.genes))
wat.exclude_ccdiff <- RunPCA(wat.exclude_ccdiff, features = c(s.genes, g2m.genes))
DimPlot(wat.exclude_cc, reduction="pca", group.by="Phase")

#REGRESS OUT CELL CYCLE SCORES DURING DATA SCALING
wat.exclude_cc <- ScaleData(wat.exclude, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(wat.exclude))
# Now, a PCA on the variable genes no longer returns components associated with cell cycle
wat.exclude_cc <- RunPCA(wat.exclude_cc, features = VariableFeatures(wat.exclude_cc), nfeatures.print = 10)

# When running a PCA on only cell cycle genes, cells no longer separate by cell-cycle phase
bat <- RunPCA(bat, features = c(s.genes, g2m.genes))
DimPlot(bat, reduction= "pca")

#RUN UMAP ON CC 
wat.exclude_cc <- FindNeighbors(wat.exclude_cc, dims = 1:10)
wat.exclude_cc <- FindClusters(wat.exclude_cc, resolution = 0.5)
wat.exclude_cc <- RunUMAP(wat.exclude_cc, dims = 1:10)
DimPlot(wat.exclude_cc, reduction = "umap",group.by="seurat_clusters")
DimPlot(wat.exclude_cc, reduction="pca", group.by="Phase")
#Find all markers for CC regression
cluster0cc.markers <- FindMarkers(wat.exclude_cc, ident.1 = 0, min.pct = 0.25)
cluster1cc.markers <- FindMarkers(wat.exclude_cc, ident.1 = 1, min.pct = 0.25)
cluster2cc.markers <- FindMarkers(wat.exclude_cc, ident.1 = 2, min.pct = 0.25)
cluster3cc.markers <- FindMarkers(wat.exclude_cc, ident.1 = 3, min.pct = 0.25)
cluster4cc.markers <- FindMarkers(wat.exclude_cc, ident.1 = 4, min.pct = 0.25)
cluster5cc.markers <- FindMarkers(wat.exclude_cc, ident.1 = 5, min.pct = 0.25)
cluster6cc.markers <- FindMarkers(wat.exclude_cc, ident.1 = 6, min.pct = 0.25)
#export markers to csv
write.csv(cluster0cc.markers, "F:/cluster0.csv")
write.csv(cluster1cc.markers, "F:/cluster1.csv")
write.csv(cluster2cc.markers, "F:/cluster2.csv")
write.csv(cluster3cc.markers, "F:/cluster3.csv")
write.csv(cluster4cc.markers, "F:/cluster4.csv")
write.csv(cluster5cc.markers, "F:/cluster5.csv")
write.csv(cluster6cc.markers, "F:/cluster6.csv")

#Heatmap for each cluster
# Remove variable feature selection
wat.exclude_cc@assays[["RNA"]]$counts <- wat.data
wat.exclude_cc <- NormalizeData(wat.exclude_cc)
wat.exclude_cc <- ScaleData(wat.exclude_cc)
#cluster 0:
cluster0.watcc <- readLines("//research.drive.wisc.edu/galmozzi/Galmozzi lab/Hoang/Projects/Bioinformatics/new_analysis/0WAT/R/ex_67/cc_regression/gene lists/watcc_cluster0_surface_highconfidence.txt")
watcc.cluster0.markers.heatmap <- DoHeatmap(wat.exclude_cc, features = cluster0.watcc, angle = 0) + scale_fill_distiller(palette = "PRGn") + theme(axis.text.y = element_text(size = 12))
save_plot(filename = "F:/watccsurfacecluster0.png", plot = watcc.cluster0.markers.heatmap,limitsize = FALSE,base_height = 20, base_width = 20)
#cluster 1
cluster1.watcc <- readLines("//research.drive.wisc.edu/galmozzi/Galmozzi lab/Hoang/Projects/Bioinformatics/new_analysis/0WAT/R/ex_67/cc_regression/gene lists/watcc_cluster1_surface_highconfidence.txt")
watcc.cluster1.markers.heatmap <- DoHeatmap(wat.exclude_cc, features = cluster1.watcc, angle = 0) + scale_fill_distiller(palette = "PRGn") + theme(axis.text.y = element_text(size = 12))
save_plot(filename = "F:/watccsurfacecluster1.png", plot = watcc.cluster1.markers.heatmap,limitsize = FALSE,base_height = 20, base_width = 20)
#cluster 2
cluster2.watcc <- readLines("//research.drive.wisc.edu/galmozzi/Galmozzi lab/Hoang/Projects/Bioinformatics/new_analysis/0WAT/R/ex_67/cc_regression/gene lists/watcc_cluster2_surface_highconfidence.txt")
watcc.cluster2.markers.heatmap <- DoHeatmap(wat.exclude_cc, features = cluster2.watcc, angle = 0) + scale_fill_distiller(palette = "PRGn") + theme(axis.text.y = element_text(size = 12))
save_plot(filename = "F:/watccsurfacecluster2.png", plot = watcc.cluster2.markers.heatmap,limitsize = FALSE,base_height = 20, base_width = 20)

#cluster 3
cluster3.watcc <- readLines("//research.drive.wisc.edu/galmozzi/Galmozzi lab/Hoang/Projects/Bioinformatics/new_analysis/0WAT/R/ex_67/cc_regression/gene lists/watcc_cluster3_surface_highconfidence.txt")
watcc.cluster3.markers.heatmap <- DoHeatmap(wat.exclude_cc, features = cluster3.watcc, angle = 0) + scale_fill_distiller(palette = "PRGn") + theme(axis.text.y = element_text(size = 12))
save_plot(filename = "F:/watccsurfacecluster3.png", plot = watcc.cluster3.markers.heatmap,limitsize = FALSE,base_height = 20, base_width = 20)
#cluster 4
cluster4.watcc <- readLines("//research.drive.wisc.edu/galmozzi/Galmozzi lab/Hoang/Projects/Bioinformatics/new_analysis/0WAT/R/ex_67/cc_regression/gene lists/watcc_cluster4_surface_highconfidence.txt")
watcc.cluster4.markers.heatmap <- DoHeatmap(wat.exclude_cc, features = cluster4.watcc, angle = 0) + scale_fill_distiller(palette = "PRGn") + theme(axis.text.y = element_text(size = 12))
save_plot(filename = "F:/watccsurfacecluster4.png", plot = watcc.cluster4.markers.heatmap,limitsize = FALSE,base_height = 20, base_width = 20)
#cluster 5
cluster5.watcc <- readLines("//research.drive.wisc.edu/galmozzi/Galmozzi lab/Hoang/Projects/Bioinformatics/new_analysis/0WAT/R/ex_67/cc_regression/gene lists/watcc_cluster5_surface_highconfidence.txt")
watcc.cluster5.markers.heatmap <- DoHeatmap(wat.exclude_cc, features = cluster5.watcc, angle = 0) + scale_fill_distiller(palette = "PRGn") + theme(axis.text.y = element_text(size = 12))
save_plot(filename = "F:/watccsurfacecluster5.png", plot = watcc.cluster5.markers.heatmap,limitsize = FALSE,base_height = 20, base_width = 20)

#cluster 6
cluster6.watcc <- readLines("//research.drive.wisc.edu/galmozzi/Galmozzi lab/Hoang/Projects/Bioinformatics/new_analysis/0WAT/R/ex_67/cc_regression/gene lists/watcc_cluster6_surface_highconfidence.txt")
watcc.cluster6.markers.heatmap <- DoHeatmap(wat.exclude_cc, features = cluster6.watcc, angle = 0) + scale_fill_distiller(palette = "PRGn") + theme(axis.text.y = element_text(size = 12))
save_plot(filename = "F:/watccsurfacecluster6.png", plot = watcc.cluster6.markers.heatmap,limitsize = FALSE,base_height = 50, base_width = 20)

#laying all clusters together for heatmap cell surface markers
watcc.markers <- readLines("//research.drive.wisc.edu/galmozzi/Galmozzi lab/Hoang/Projects/Bioinformatics/new_analysis/0WAT/R/ex_67/cc_regression/gene lists/surface_combined_for_heatmap.txt")
watcc.markers.heatmap <- DoHeatmap(wat.exclude_cc, features = watcc.markers, angle = 0) + scale_fill_distiller(palette = "PRGn") + theme(axis.text.y = element_text(size = 4))
save_plot(filename = "F:/watccsurfacemap.png", plot = watcc.markers.heatmap,limitsize = FALSE,base_height = 50, base_width = 30)
#heatmap cell surface markers without duplicates
watcc.markers_duperm <- readLines("//research.drive.wisc.edu/galmozzi/Galmozzi lab/Hoang/Projects/Bioinformatics/new_analysis/0WAT/R/ex_67/cc_regression/gene lists/surface_combined_for_heatmap_DUPE_removed.txt")
watcc.markers.heatmap_duperm <- DoHeatmap(wat.exclude_cc, features = watcc.markers_duperm, angle = 0) + scale_fill_distiller(palette = "PRGn") + theme(axis.text.y = element_text(size = 8))
save_plot(filename = "F:/watccsurfacemap_DUPErm.png", plot = watcc.markers.heatmap_duperm,limitsize = FALSE,base_height = 50, base_width = 30)


#export uMAP projections
#1. for unexcluded data:
umap.unexcluded = cbind("Barcode" = rownames(Embeddings(object = wat, reduction = "umap")), Embeddings(object = wat, reduction = "umap"))
write.table(umap.unexcluded, file="F:/uMAP_seurat.csv", sep = ",", quote = F, row.names = F, col.names = T)
#2. For excluded data
umap = cbind("Barcode" = rownames(Embeddings(object = exclude, reduction = "umap")), Embeddings(object = exclude, reduction = "umap"))
write.table(umap, file="F:/uMAP_seurat.csv", sep = ",", quote = F, row.names = F, col.names = T)
#3. for CCreg data
umap = cbind("Barcode" = rownames(Embeddings(object = wat.exclude_cc, reduction = "umap")), Embeddings(object = wat.exclude_cc, reduction = "umap"))
write.table(umap, file="F:/wat.exclude_cc_UMAP.csv", sep = ",", quote = F, row.names = F, col.names = T)

#export metadata to get clusters export
#1. for unexcluded data
write.csv(wat@meta.data,"F:/seurat_metadata.csv")
write.csv(exclude@meta.data,"F:/seurat_metadata.csv")
#3. for ccreg data
write.csv(wat.exclude_cc@meta.data,"F:/seurat_metadata.csv")


# To subset and remove single cluster and keep the remaining clusters for new analysis
#sub_obj <- subset(wat, idents = 1, invert = TRUE)

# Can also provide multiple idents using the typical R syntax "c()"
#sub_obj <- subset(object = obj_name, idents = c(1, 2, 3))

# If you have renamed the clusters be sure to provide their names in quotes in the function
#sub_obj <- subset(object = obj_name, idents = "Monocytes")

#create .loom file from wat.exclude object
#watclusters.loom <- as.loom(wat.exclude, filename = "F:/watclusters.loom", verbose = FALSE)
#watclusters.loom$close_all()

#create H5ad object (loom doesnt work)
#SaveH5Seurat(wat, filename="watclusters.h5Seurat")
