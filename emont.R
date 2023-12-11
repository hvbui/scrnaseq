# integrated analysis using RPCA
library(Seurat)
library(SeuratData)
library(glmGampoi)
library(Seurat.utils)
AG_WAT <- readRDS("//research.drive.wisc.edu/galmozzi/Galmozzi lab/Hoang/Projects/Adipocyte Subtypes/new_analysis/0WAT/R/ex_67/cc_regression/seurat objects/wat.exclude_cc.rds")
Emont_WAT <- readRDS("//research.drive.wisc.edu/galmozzi/Galmozzi lab/Hoang/Projects/Adipocyte Subtypes/integrated_analysis/Emont_2022/mouse_ASPCs.rds")

AG_WAT <- wat.exclude_cc
#seems like Emont has capitalized gene names. lets correct it.

# Extract gene names from the Dimnames slot and capitalize them
capitalized_gene_names <- toupper(rownames(Emont_WAT$RNA$counts))
# Update the row names in the Seurat object
rownames(Emont_WAT$RNA$counts) <- capitalized_gene_names

# Extract gene names from the Dimnames slot and capitalize them
capitalized_gene_names <- toupper(rownames(Emont_WAT@assays$RNA@counts))

# Update the row names in the Seurat object
Emont_WAT@assays$RNA@counts@Dimnames[[1]] <- capitalized_gene_names
Emont_WAT@assays$RNA@data@Dimnames[[1]] <- capitalized_gene_names
Emont_WAT@assays$RNA@meta.features <- data.frame(genes = capitalized_gene_names)

#attempting to fix error by SCTransform
AG_WAT <-SCTransform(AG_WAT)
Emont_WAT <-SCTransform(Emont_WAT)
AG_SCT <- AG_WAT
Emont_SCT <- Emont_WAT

#repeat capitalize gene names for SCT assay
# Extract gene names from the Dimnames slot and capitalize them
capitalized_gene_names <- toupper(rownames(Emont_WAT@assays$SCT@counts))

# Update the row names in the Seurat object
Emont_WAT@assays$SCT@counts@Dimnames[[1]] <- capitalized_gene_names
Emont_WAT@assays$SCT@data@Dimnames[[1]] <- capitalized_gene_names
Emont_WAT@assays$SCT@meta.features <- data.frame(genes = capitalized_gene_names)
#same goes for scale.data, meta.features, var.features slot
capitalized_gene_names <- toupper(rownames(Emont_WAT@assays$SCT@scale.data))# this capitalize scale.data genes
rownames(Emont_WAT@assays$SCT@scale.data) <- capitalized_gene_names 

capitalized_gene_names <- toupper(Emont_WAT@assays$SCT@var.features)
Emont_WAT@assays$SCT@var.features <-capitalized_gene_names

capitalized_gene_names <- toupper(rownames(Emont_WAT@assays$SCT@data)
rownames(Emont_WAT@assays$SCT@data)<- capitalized_gene_names

#creating a list of 2 seurat object
AG_Emont.list <-list(AG_WAT, Emont_WAT)

AG_Emont.list <- SplitObject(AG_Emont.list, split.by = "cell_type") #this code doesnt work

# normalize and identify variable features for each dataset independently
AG_Emont.list <- lapply(X = AG_Emont.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- SelectIntegrationFeatures(object.list = AG_Emont.list, nfeatures=2000)

AG_Emont.list <- lapply(X = AG_Emont.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

#Perform integration
#find anchor
AG_Emont.anchors <- FindIntegrationAnchors(object.list = AG_Emont.list, anchor.features = features, reduction = "rpca", dims=1:50)
AG_Emont.AG_Emont.combined <- IntegrateData(anchorset = AG_Emont.anchors, dims=1:50)

###transfer label
# Create metadata columns
AG_WAT_clusters <- Idents(AG_WAT) 
Emont_WAT_clusters <- Idents(Emont_WAT)

# Transfer labels
AG_Emont.combined[['AG_WAT_clusters']] <- "no match"
AG_Emont.combined[['Emont_WAT_clusters']] <- "no match"

for(i in 1:ncol(AG_Emont.combined)){
  AG_Emont.combined <- SetIdent(AG_Emont.combined, cells = colnames(AG_Emont.combined)[i], ident = AG_WAT_clusters[match(colnames(AG_Emont.combined)[i], colnames(AG_WAT))] )
  
  AG_Emont.combined <- SetIdent(AG_Emont.combined, cells = colnames(AG_Emont.combined)[i], ident = Emont_WAT_clusters[match(colnames(AG_Emont.combined)[i], colnames(Emont_WAT))] )
}

DimPlot(AG_Emont.combined, group.by = 'Emont_WAT_clusters')
###
# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(AG_Emont.AG_Emont.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
AG_Emont.AG_Emont.combined <- ScaleData(AG_Emont.AG_Emont.combined, verbose = FALSE)
AG_Emont.AG_Emont.combined <- RunPCA(AG_Emont.AG_Emont.combined, npcs = 10, verbose = FALSE)
AG_Emont.AG_Emont.combined <- RunUMAP(AG_Emont.AG_Emont.combined, reduction = "pca", dims = 1:10)
AG_Emont.AG_Emont.combined <- FindNeighbors(AG_Emont.AG_Emont.combined, reduction = "pca", dims = 1:10)
AG_Emont.AG_Emont.combined <- FindClusters(AG_Emont.AG_Emont.combined, resolution = 0.5)

wat.exclude_ccUMAP <- DimPlot(wat.exclude_cc)
AG_Emont_combined_AGlabel <- DimPlot(AG_Emont.AG_Emont.combined, group.by = "RNA_snn_res.0.5"
AG_Emont_combined_Emontlabel <- DimPlot(AG_Emont.AG_Emont.combined, group.by = "cell_type")
AG_Emont_combined_Seuratlabel <- DimPlot(AG_Emont.AG_Emont.combined, group.by= "seurat_clusters")
#save da plot
save_plot(filename = "F:/AG_Emont_combined_AGlabel.png", plot = AG_Emont_combined_AGlabel,limitsize = FALSE,base_height = 5, base_width = 5)
save_plot(filename = "F:/AG_Emont_combined_Emontlabel.png", plot = AG_Emont_combined_Emontlabel,limitsize = FALSE,base_height = 5, base_width = 5.5)
save_plot(filename = "F:/AG_Emont_combined_Seuratlabel.png", plot = AG_Emont_combined_Seuratlabel,limitsize = FALSE,base_height = 5, base_width = 5)
save_plot(filename = "F:/wat.exclude_cc_UPDATE.png", plot = wat.exclude_ccUMAP,limitsize = FALSE,base_height = 5, base_width = 5)
