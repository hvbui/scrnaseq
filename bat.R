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
library(stringr)
require("Seurat.utils")
#read in data
bat.data = Read10X(data.dir = "//research.drive.wisc.edu/galmozzi/Galmozzi lab/Hoang/Projects/Adipocyte Subtypes/new_analysis/0BAT/filtered_feature_bc_matrix")
#use this read in for capitalized gene names (for regression)
bat.data = Read10X(data.dir = "F:/new_analysis/0BAT/outs/filtered_feature_bc_matrix_CAP")

#convert into seurat object
bat = CreateSeuratObject(counts = bat.data, project = "0BATpreads", min.cells = 3, min.features = 200)
bat$data.set = rep("bat", length(bat$orig.ident))

#convert new names too for regression
batcc = CreateSeuratObject(counts = bat.data, project = "0BATpreads", min.cells = 3, min.features = 200)
batcc$data.set = rep("batcc", length(batcc$orig.ident))

#calculate percent mitochondrial RNA for batpreads
bat[["percent.mt"]] = PercentageFeatureSet(bat, pattern = "^mt-")
bat[["percent.mt"]] = PercentageFeatureSet(bat, pattern = "^MT-")
#for cc use ^MT because genes are capitalized. idk why but it works like that
batcc[["percent.mt"]] = PercentageFeatureSet(batcc, pattern = "^MT-")

#generate violin plot #genes/cell, mtRNA/cell, reads/cell
VlnPlot(bat, features= c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size=0.5, ncol = 3)
#rerun without the dots
VlnPlot(bat, features= c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size=0, ncol = 3)

#Feature Scatter visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(bat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(bat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#apply filter-- UMI 6000-60000, features 3900-7700, mito umi 1.7-5
bat <- subset(bat, subset = nFeature_RNA > 3900 & nFeature_RNA < 7700 & percent.mt < 5& percent.mt > 1.7 & nCount_RNA > 6000 & nCount_RNA < 60000)
batcc <- subset(batcc, subset = nFeature_RNA > 3900 & nFeature_RNA < 7700 & percent.mt < 5& percent.mt > 1.7 & nCount_RNA > 6000 & nCount_RNA < 60000)

#normalizing the data -normalizes the feature expression measurements for each cell by the total expression, multiplies this by a scale factor (10,000 by default), and log-transforms the result
bat <- NormalizeData(bat, normalization.method = "LogNormalize", scale.factor = 10000)
batcc <- NormalizeData(batcc, normalization.method = "LogNormalize", scale.factor = 10000)
#identification of highly variable features
bat <- FindVariableFeatures(bat, selection.method = "vst", nfeatures = 2000)
batcc <- FindVariableFeatures(batcc, selection.method = "vst", nfeatures = 2000)

#Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(bat), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(bat)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#scaling the data -pre-processing step prior to dimensional reduction techniques like PCA. The ScaleData() function:
#Shifts the expression of each gene, so that the mean expression across cells is 0
#Scales the expression of each gene, so that the variance across cells is 1
  #This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
  #The results of this are stored in pbmc[["RNA"]]@scale.data 
all.genes <- rownames(bat)
bat <- ScaleData(bat, features = all.genes)

all.genes <- rownames(batcc)
batcc <- ScaleData(batcc, features = all.genes)


#PCA
bat <- RunPCA(bat, features = VariableFeatures(object = bat))
batcc <- RunPCA(batcc, features = VariableFeatures(object = batcc))


## CELL CYCLE REGRESSION
# Running a PCA on cell cycle genes reveals, unsurprisingly, that cells separate entirely by
# phase
cc.genes$s.genes <- str_to_title(cc.genes$s.genes)
cc.genes$g2m.genes <- str_to_title(cc.genes$g2m.genes)

#cc.genes <- as.character(cc.genes)
#cc.genes <- str_to_title(cc.genes)
#s.genes <- cc.genes$s.genes
s.genes <- cc.genes["s.genes"]
#g2m.genes <- cc.genes$g2m.genes
g2m.genes <- cc.genes["g2m.genes"]


#DIRECTLY FROM SATJIA LAB
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
bat <- CellCycleScoring(bat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
#uncapitalizing cc genes
m.s.genes <- str_to_title(cc.genes$s.genes)
m.g2m.genes <- str_to_title(cc.genes$g2m.genes)

bat <- CellCycleScoring(bat, s.features = m.s.genes, g2m.features = m.g2m.genes, set.ident = TRUE)
# view cell cycle scores and phase assignments
head(bat[[]])

#run PCA on cell cycle 
bat <- RunPCA(bat, features = c(s.genes, g2m.genes))
#for uncapitalized cc genes:
bat <- RunPCA(bat, features = c(m.s.genes, m.g2m.genes))
DimPlot(bat)
# Visualize the distribution of cell cycle markers across
RidgePlot(bat, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)

#REGRESS OUT CELL CYCLE SCORES DURING DATA SCALING
bat <- ScaleData(bat, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(bat))
# Now, a PCA on the variable genes no longer returns components associated with cell cycle
bat <- RunPCA(bat, features = VariableFeatures(bat), nfeatures.print = 10)

# When running a PCA on only cell cycle genes, cells no longer separate by cell-cycle phase
bat <- RunPCA(bat, features = c(s.genes, g2m.genes))
DimPlot(bat, reduction= "pca")


#Find all markers from new cells (run clustering line 124 and uMAP line 128 first)
cluster0.markers <- FindMarkers(bat, ident.1 = 0, min.pct = 0.25)
cluster1.markers <- FindMarkers(bat, ident.1 = 1, min.pct = 0.25)
cluster2.markers <- FindMarkers(bat, ident.1 = 2, min.pct = 0.25)
cluster3.markers <- FindMarkers(bat, ident.1 = 3, min.pct = 0.25)
cluster4.markers <- FindMarkers(bat, ident.1 = 4, min.pct = 0.25)
cluster5.markers <- FindMarkers(bat, ident.1 = 5, min.pct = 0.25)


#####################################################
#determine dimensionality
ElbowPlot(bat, ndims = 50, reduction = "pca")
bat <- JackStraw(bat, num.replicate=100)
bat <- ScoreJackStraw(bat, dims = 1:20)
JackStrawPlot(bat, dims = 1:15)

#clustering
bat <- FindNeighbors(bat, dims = 1:10)
bat <- FindClusters(bat, resolution = 0.5)

batcc <- FindNeighbors(batcc, dims = 1:10)
batcc <- FindClusters(batcc, resolution = 0.5)
#umap
bat <- RunUMAP(bat, dims = 1:10)
DimPlot(bat, reduction = "umap")
DimPlot(bat, reduction = "umap", group.by="Phase")
batcc <- RunUMAP(batcc, dims = 1:10)
DimPlot(batcc, reduction = "umap")

#save object so it doesnt need to rerun all previous codes
#saveRDS(bat, file ="//research.drive.wisc.edu/galmozzi/Galmozzi lab/Hoang/Projects/Bioinformatics/new_analysis/0BAT/R/23-1-10.rds")

#find all markers
#1. For unexcluded data
markers.noexclude <- FindAllMarkers(bat,features = VariableFeatures(object = bat))
write.csv(markers.noexclude, "F:/markers_noexclude.csv")
#2. For excluded data
markers.exclude <- FindAllMarkers(exclude,features = VariableFeatures(object = exclude))
write.csv(markers.exclude, "F:/markers_exclude.csv")

# find all markers of cluster 2
cluster2.markers <- FindMarkers(bat, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)

#get expression data for ALL genes
#1. For unexcluded data
expression_data <- AverageExpression(bat)
write.csv(expression_data, "F:/expression_data.csv")
#2. For excluded data
expression_data_exclude <- AverageExpression(exclude)
write.csv(expression_data_exclude, "F:/expression_data.csv")


#remove cluster 7 8 9 and recluster
bat.exclude <- subset(x=bat, idents = c(7,8,9), invert = TRUE)
bat.exclude <- NormalizeData(bat.exclude, normalization.method = "LogNormalize", scale.factor = 10000)
bat.exclude <- FindVariableFeatures(bat.exclude, selection.method = "vst", nfeatures = 2000)
bat.exclude <- RunPCA(bat.exclude, features = VariableFeatures(object = bat.exclude))
bat.exclude <- FindNeighbors(bat.exclude, dims = 1:10)
bat.exclude <- FindClusters(bat.exclude, resolution = 0.5)
bat.exclude <- RunUMAP(bat.exclude, dims = 1:10)
DimPlot(bat.exclude, reduction = "umap")

saveRDS(bat.exclude, file= "F:/bat.exclude_uncap.rds")

#generating heatmap for cell surface markers
#cluster 0:
cluster0.bat <- readLines("//research.drive.wisc.edu/galmozzi/Galmozzi lab/Hoang/Projects/Bioinformatics/new_analysis/0BAT/R/ex_789/gene list/bat_cluster0_surface.txt")
bat.cluster0.markers.heatmap <- DoHeatmap(bat.exclude, features = cluster0.bat, angle = 0) + scale_fill_distiller(palette = "PRGn") + theme(axis.text.y = element_text(size = 12))
save_plot(filename = "F:/batsurfacecluster0.png", plot = bat.cluster0.markers.heatmap,limitsize = FALSE,base_height = 20, base_width = 20)
#cluster 1
cluster1.bat <- readLines("//research.drive.wisc.edu/galmozzi/Galmozzi lab/Hoang/Projects/Bioinformatics/new_analysis/0BAT/R/ex_789/gene list/bat_cluster1_surface.txt")
bat.cluster1.markers.heatmap <- DoHeatmap(bat.exclude, features = cluster1.bat, angle = 0) + scale_fill_distiller(palette = "PRGn") + theme(axis.text.y = element_text(size = 12))
save_plot(filename = "F:/batsurfacecluster1.png", plot = bat.cluster1.markers.heatmap,limitsize = FALSE,base_height = 20, base_width = 20)
#cluster 2
cluster2.bat <- readLines("//research.drive.wisc.edu/galmozzi/Galmozzi lab/Hoang/Projects/Bioinformatics/new_analysis/0BAT/R/ex_789/gene list/bat_cluster2_surface.txt")
bat.cluster2.markers.heatmap <- DoHeatmap(bat.exclude, features = cluster2.bat, angle = 0) + scale_fill_distiller(palette = "PRGn") + theme(axis.text.y = element_text(size = 12))
save_plot(filename = "F:/batsurfacecluster2.png", plot = bat.cluster2.markers.heatmap,limitsize = FALSE,base_height = 20, base_width = 20)

#cluster 3
cluster3.bat <- readLines("//research.drive.wisc.edu/galmozzi/Galmozzi lab/Hoang/Projects/Bioinformatics/new_analysis/0BAT/R/ex_789/gene list/bat_cluster3_surface.txt")
bat.cluster3.markers.heatmap <- DoHeatmap(bat.exclude, features = cluster3.bat, angle = 0) + scale_fill_distiller(palette = "PRGn") + theme(axis.text.y = element_text(size = 12))
save_plot(filename = "F:/batsurfacecluster3.png", plot = bat.cluster3.markers.heatmap,limitsize = FALSE,base_height = 20, base_width = 20)
#cluster 4
cluster4.bat <- readLines("//research.drive.wisc.edu/galmozzi/Galmozzi lab/Hoang/Projects/Bioinformatics/new_analysis/0BAT/R/ex_789/gene list/bat_cluster4_surface.txt")
bat.cluster4.markers.heatmap <- DoHeatmap(bat.exclude, features = cluster4.bat, angle = 0) + scale_fill_distiller(palette = "PRGn") + theme(axis.text.y = element_text(size = 12))
save_plot(filename = "F:/batsurfacecluster4.png", plot = bat.cluster4.markers.heatmap,limitsize = FALSE,base_height = 20, base_width = 20)
#cluster 5
cluster5.bat <- readLines("//research.drive.wisc.edu/galmozzi/Galmozzi lab/Hoang/Projects/Bioinformatics/new_analysis/0BAT/R/ex_789/gene list/bat_cluster5_surface.txt")
bat.cluster5.markers.heatmap <- DoHeatmap(bat.exclude, features = cluster5.bat, angle = 0) + scale_fill_distiller(palette = "PRGn") + theme(axis.text.y = element_text(size = 12))
save_plot(filename = "F:/batsurfacecluster5.png", plot = bat.cluster5.markers.heatmap,limitsize = FALSE,base_height = 20, base_width = 20)

#laying all clusters together for heatmap cell surface markers
bat.surface.markers <- readLines("//research.drive.wisc.edu/galmozzi/Galmozzi lab/Hoang/Projects/Bioinformatics/new_analysis/0BAT/R/ex_789/gene list/bat_combined_surface.txt")
bat.markers.surface.heatmap <- DoHeatmap(bat.exclude, features = bat.surface.markers, angle = 0) + scale_fill_distiller(palette = "PRGn") + theme(axis.text.y = element_text(size = 12))
save_plot(filename = "F:/batsurfacemap.png", plot = bat.markers.surface.heatmap,limitsize = FALSE,base_height = 50, base_width = 30)
#heatmap cell surface markers without duplicates
bat.markers_duperm <- readLines("//research.drive.wisc.edu/galmozzi/Galmozzi lab/Hoang/Projects/Bioinformatics/new_analysis/0BAT/R/ex_789/gene list/bat_combined_surface_dupe_removed.txt")
bat.markers.heatmap_duperm <- DoHeatmap(bat.exclude, features = bat.markers_duperm, angle = 0) + scale_fill_distiller(palette = "PRGn") + theme(axis.text.y = element_text(size = 12))
save_plot(filename = "F:/batsurfacemap_DUPErm.png", plot = bat.markers.heatmap_duperm,limitsize = FALSE,base_height = 20, base_width = 20)


#ALTERNATE WORKFLOW- CELL CYCLE DIFFERENCE REGRESSION
exclude$CC.Difference <- exclude$S.Score - exclude$G2M.Score
exclude_ccdiff <- ScaleData(exclude, vars.to.regress = "CC.Difference", features = rownames(exclude))
# cell cycle effects strongly mitigated in PCA
exclude <- RunPCA(exclude, features = VariableFeatures(exclude), nfeatures.print = 10)

# when running a PCA on cell cycle genes, actively proliferating cells remain distinct from G1
# cells however, within actively proliferating cells, G2M and S phase cells group together
#group.by recognizes all metadata column names
exclude <- RunPCA(exclude, features = c(s.genes, g2m.genes))
exclude_ccdiff <- RunPCA(exclude_ccdiff, features = c(s.genes, g2m.genes))
DimPlot(exclude, reduction="pca", group.by="Phase")
#then run UMAP
exclude <- FindNeighbors(exclude, dims = 1:10)
exclude <- FindClusters(exclude, resolution = 0.5)
exclude <- RunUMAP(exclude, dims = 1:10)

exclude_ccdiff <- FindNeighbors(exclude_ccdiff, dims = 1:10)

exclude_ccdiff <- FindClusters(exclude_ccdiff, resolution = 0.5)

exclude_ccdiff <- RunUMAP(exclude_ccdiff, dims = 1:10)
DimPlot(exclude, reduction = "umap")
DimPlot(exclude_ccdiff, reduction = "umap")

#Find all markers from new cells - removing cluster 7 8 9
cluster0.markers <- FindMarkers(exclude, ident.1 = 0, min.pct = 0.25)
cluster1.markers <- FindMarkers(exclude, ident.1 = 1, min.pct = 0.25)
cluster2.markers <- FindMarkers(exclude, ident.1 = 2, min.pct = 0.25)
cluster3.markers <- FindMarkers(exclude, ident.1 = 3, min.pct = 0.25)
cluster4.markers <- FindMarkers(exclude, ident.1 = 4, min.pct = 0.25)
cluster5.markers <- FindMarkers(exclude, ident.1 = 5, min.pct = 0.25)
#export cluster genes
write.csv(cluster0.markers, "F:/cluster0.csv")
write.csv(cluster1.markers, "F:/cluster1.csv")
write.csv(cluster2.markers, "F:/cluster2.csv")
write.csv(cluster3.markers, "F:/cluster3.csv")
write.csv(cluster4.markers, "F:/cluster4.csv")
write.csv(cluster5.markers, "F:/cluster5.csv")

#export gene names and normalized expression values

# Extract the gene names
bat_gene_names <- rownames(bat@assays$RNA@counts)

# Extract the normalized expression values
normalized_expression_values <- bat@meta.data

# Combine gene names and normalized expression values into a data frame
df <- data.frame(normalized_expression_values, row.names = bat_gene_names)
# Save the data frame to a .csv file
write.csv(df, file = "gene_expression_data.csv")


#export uMAP projections
#1. for unexcluded data:
umap.unexcluded = cbind("Barcode" = rownames(Embeddings(object = bat, reduction = "umap")), Embeddings(object = bat, reduction = "umap"))
write.table(umap.unexcluded, file="F:/uMAP_seurat.csv", sep = ",", quote = F, row.names = F, col.names = T)
#2. For excluded data
umap = cbind("Barcode" = rownames(Embeddings(object = exclude, reduction = "umap")), Embeddings(object = exclude, reduction = "umap"))
write.table(umap, file="F:/uMAP_seurat.csv", sep = ",", quote = F, row.names = F, col.names = T)

#export metadata to get clusters export
#1. for unexcluded data
write.csv(bat@meta.data,"F:/seurat_metadata.csv")
write.csv(exclude@meta.data,"F:/seurat_metadata.csv")

# To subset and remove single cluster and keep the remaining clusters for new analysis
#sub_obj <- subset(bat, idents = 1, invert = TRUE)

# Can also provide multiple idents using the typical R syntax "c()"
#sub_obj <- subset(object = obj_name, idents = c(1, 2, 3))

# If you have renamed the clusters be sure to provide their names in quotes in the function
#sub_obj <- subset(object = obj_name, idents = "Monocytes")

#create .loom file from exclude object
#batclusters.loom <- as.loom(exclude, filename = "F:/batclusters.loom", verbose = FALSE)
#batclusters.loom$close_all()

#create H5ad object (loom doesnt work)
#SaveH5Seurat(bat, filename="batclusters.h5Seurat")

meta.data.groups <- unique(x = bat@meta.data$seurat_clusters)

#export expression matrix
expression_matrix <- as.data.frame(GetAssayData(object = bat, slot = "scale.data", assay = "RNA"))
cluster_names <- bat@meta.data[group.cells, "seurat_clusters"]
expression_matrix_with_cluster_names <- cbind(cluster_names, expression_matrix)
write.csv(expression_matrix_with_cluster_names, file = "F:/expression_matrix_with_cluster_names.csv")

#set idents
Idents(object=exclude) <- "seurat_clusters"
Idents(object=exclude_ccdiff) <- "seurat_clusters"

#find DE markers for ccdiff 3 and 5
ccdiff.markers <-FindMarkers(exclude_ccdiff, ident.1 = "3", ident.2 = "5")
write.csv(ccdiff.markers, "F:/ccdiff.csv")
#find DE markers for exclude 2 and 3
cc23.exclude.markers <-FindMarkers(exclude, ident.1 = "2", ident.2 = "3")
write.csv(cc23.exclude.markers, "F:/2and3.exclude.markers.csv")

#find DE markers for cluster 0 vs 1 exclude
de.0v1.markers <-FindMarkers(exclude, ident.1 = "0", ident.2 = "1")
write.csv(de.0v1.markers, "F:/de.0v1.markers.csv")

#find DE markers for 0 vs 5 exclude
de.0v5.markers <-FindMarkers(exclude, ident.1 = "0", ident.2 = "5")
write.csv(de.0v5.markers, "F:/de.0v5.markers.csv")

#find DE markers for 1 vs 4 exclude
de.1v4.markers <-FindMarkers(exclude, ident.1 = "1", ident.2 = "4")
write.csv(de.1v4.markers, "F:/de.1v4.markers.csv")

#find DE markers for 4v5  exclude
de.4v5.markers <-FindMarkers(exclude, ident.1 = "4", ident.2 = "5")
write.csv(de.4v5.markers, "F:/de.4v5.markers.csv")

#FIND MARKERS WITH CLUSTERING(chatGPT)
de_genes_exclude <- FindMarkers(exclude, ident.1 = 2, ident.2 = 3)

cluster_assignments_exclude <- colData(exclude)$seurat_clusters

df_de_genes <- data.frame(gene = rownames(de_genes),
                          logFC = de_genes$logFC,
                          p_val = de_genes$p_val,
                          cluster_1 = rep("cluster 2", nrow(de_genes)),
                          cluster_2 = rep("cluster 3", nrow(de_genes)))

#find conserved markers for ccdiff 3 and 5
library(multtest)
library(metap)
ccdiff.markers_conserved <-FindConservedMarkers(exclude_ccdiff, ident.1 = 3, ident.2 = 5, grouping.var="seurat_clusters")
write.csv(ccdiffconserverd.markers, "F:/ccdiffconserved.csv")

#just do it manually

#save expression value, and then add in cluster assignment from full_seurat_cluster.csv
#define genes of interest
genes_bat_full <- c("MYOD1","PLP1","ITGAM", "PTPRC", "PDGFRA", "ALDH1A1", "PPARG", "FABP4", "DLK1","FLT1" )
genes_wat_full <-c("PLP1","ITGAM", "PTPRC", "PDGFRA","PPARG", "FABP4", "LY6C1","FLT1")
# Get the expression values for the specified genes
expression_values_batfull <- GetAssayData(bat, assay = "RNA")[genes_bat_full, ]
expression_values_watfull <- GetAssayData(wat, assay = "RNA")[genes_wat_full, ]
expression_values_batex <- GetAssayData(bat.exclude, assay = "RNA")[genes_bat_full, ]
expression_values_watex <- GetAssayData(wat.exclude, assay = "RNA")[genes_wat_full, ]
#write csv
write.csv(expression_values_batex, "F:/expression_values_batex")
write.csv(expression_values_watex, "F:/expression_values_watex")
