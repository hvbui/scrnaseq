# integrated analysis using RPCA
library(Seurat)
library(SeuratData)
library(glmGampoi)
library(Seurat.utils)
AG_WAT <- readRDS("//research.drive.wisc.edu/galmozzi/Galmozzi lab/Hoang/Projects/Adipocyte Subtypes/new_analysis/0WAT/R/ex_67/cc_regression/seurat objects/wat.exclude_cc.rds")
Emont_WAT <- readRDS("//research.drive.wisc.edu/galmozzi/Galmozzi lab/Hoang/Projects/Adipocyte Subtypes/integrated_analysis/Sarvari/eWAT_FAP.Rds")

#update seurat object
Sarvari_WAT = UpdateSeuratObject(Emont_WAT)

#attempting to fix error by SCTransform
AG_WAT <-SCTransform(AG_WAT)
Sarvari_WAT <-SCTransform(Sarvari_WAT)
AG_SCT <- AG_WAT
Sarvari_SCT <- Sarvari_WAT

#creating a list of 2 seurat object
AG_Sarvari.list <-list(AG_WAT, Sarvari_WAT)


# normalize and identify variable features for each dataset independently
AG_Sarvari.list <- lapply(X = AG_Sarvari.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- SelectIntegrationFeatures(object.list = AG_Sarvari.list, nfeatures=2000)

AG_Sarvari.list <- lapply(X = AG_Sarvari.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

#Perform integration
#find anchor
AG_Sarvari.anchors <- FindIntegrationAnchors(object.list = AG_Sarvari.list, anchor.features = features, reduction = "rpca", dims=1:50)
AG_Sarvari.AG_Sarvari.combined <- IntegrateData(anchorset = AG_Sarvari.anchors, dims=1:50)

###transfer label
# Create metadata columns
AG_WAT_clusters <- Idents(AG_WAT) 
Sarvari_WAT_clusters <- Idents(Sarvari_WAT)

# Transfer labels
AG_Sarvari.combined[['AG_WAT_clusters']] <- "no match"
AG_Sarvari.combined[['Sarvari_WAT_clusters']] <- "no match"

for(i in 1:ncol(AG_Sarvari.combined)){
  AG_Sarvari.combined <- SetIdent(AG_Sarvari.combined, cells = colnames(AG_Sarvari.combined)[i], ident = AG_WAT_clusters[match(colnames(AG_Sarvari.combined)[i], colnames(AG_WAT))] )
  
  AG_Sarvari.combined <- SetIdent(AG_Sarvari.combined, cells = colnames(AG_Sarvari.combined)[i], ident = Sarvari_WAT_clusters[match(colnames(AG_Sarvari.combined)[i], colnames(Sarvari_WAT))] )
}

DimPlot(AG_Sarvari.AG_Sarvari.combined, group.by = 'Sarvari_WAT_clusters')
###
# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(AG_Sarvari.AG_Sarvari.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
AG_Sarvari.AG_Sarvari.combined <- ScaleData(AG_Sarvari.AG_Sarvari.combined, verbose = FALSE)
AG_Sarvari.AG_Sarvari.combined <- RunPCA(AG_Sarvari.AG_Sarvari.combined, npcs = 10, verbose = FALSE)
AG_Sarvari.AG_Sarvari.combined <- RunUMAP(AG_Sarvari.AG_Sarvari.combined, reduction = "pca", dims = 1:10)
AG_Sarvari.AG_Sarvari.combined <- FindNeighbors(AG_Sarvari.AG_Sarvari.combined, reduction = "pca", dims = 1:10)
AG_Sarvari.AG_Sarvari.combined <- FindClusters(AG_Sarvari.AG_Sarvari.combined, resolution = 0.5)

wat.exclude_ccUMAP <- DimPlot(wat.exclude_cc)
AG_Sarvari_combined_AGlabel <- DimPlot(AG_Sarvari.AG_Sarvari.combined, group.by = "RNA_snn_res.0.5")
AG_Sarvari_combined_Sarvarilabel <- DimPlot(AG_Sarvari.AG_Sarvari.combined, group.by = "Subtype")
AG_Sarvari_combined_Seuratlabel <- DimPlot(AG_Sarvari.AG_Sarvari.combined, group.by= "seurat_clusters")
#save da plot
save_plot(filename = "F:/AG_Sarvari_combined_AGlabel.png", plot = AG_Sarvari_combined_AGlabel,limitsize = FALSE,base_height = 5, base_width = 5)
save_plot(filename = "F:/AG_Sarvari_combined_Sarvarilabel.png", plot = AG_Sarvari_combined_Sarvarilabel,limitsize = FALSE,base_height = 5, base_width = 5.5)
save_plot(filename = "F:/AG_Sarvari_combined_Seuratlabel.png", plot = AG_Sarvari_combined_Seuratlabel,limitsize = FALSE,base_height = 5, base_width = 5)
save_plot(filename = "F:/wat.exclude_cc_UPDATE.png", plot = wat.exclude_ccUMAP,limitsize = FALSE,base_height = 5, base_width = 5)
                                     