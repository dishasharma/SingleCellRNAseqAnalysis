IntegrationSeurat <- function(QCedRDS,Dims,saveDir, Assay = "RNA", samplename, condition, splitby, resolution)
{
message(paste("\nLoading Libraries\n"))
library(Seurat)
library(magrittr)
library(dplyr)
library(patchwork)
library(tibble)
library(AnnotationDbi)
library(ggplot2)
library(stringr)
library(dplyr)
library(clustree)

message(paste("\nSpliting Object\n"))
QCedRDS.list <- SplitObject(QCedRDS, split.by = splitby)

message(paste("\nFinding Anchors\n"))
# normalize and identify variable features for each dataset independently
QCedRDS.list <- lapply(X = QCedRDS.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = QCedRDS.list)
QCedRDS.anchors <- FindIntegrationAnchors(object.list = QCedRDS.list, anchor.features = features)
QCedRDS.combined <- IntegrateData(anchorset = QCedRDS.anchors)
DefaultAssay(QCedRDS.combined) <- "integrated"

message(paste("\nScaling and Processing\n"))
# Run the standard workflow for visualization and clustering
QCedRDS.combined <- ScaleData(QCedRDS.combined, verbose = FALSE)
QCedRDS.combined <- RunPCA(QCedRDS.combined, npcs = Dims, verbose = FALSE)
QCedRDS.combined <- RunUMAP(QCedRDS.combined, reduction = "pca", dims = 1:Dims)
QCedRDS.combined <- FindNeighbors(QCedRDS.combined, reduction = "pca", dims = 1:Dims)
QCedRDS.combined <- FindClusters(QCedRDS.combined, resolution = resolution)

dir.create(paste(saveDir,"UMAP/",sep = ""), showWarnings = FALSE)
pdf(paste(saveDir,"UMAP/umap_splitby",samplename,"_",Assay,"_",splitby,"_",condition,".pdf",sep = ""), width = 7, height = 12)
print(DimPlot(QCedRDS.combined, reduction = "umap", group.by = splitby))
dev.off()

dir.create(paste(saveDir,"UMAP/",sep = ""), showWarnings = FALSE)
pdf(paste(saveDir,"UMAP/umap_",samplename,"_",Assay,"_",condition,".pdf",sep = ""), width = 7, height = 12)
print(DimPlot(QCedRDS.combined, reduction = "umap", label = TRUE, repel = TRUE))
dev.off()
message("Saving RDS file")
dir.create(paste(saveDir,"saveRDS_obj/",sep = ""), showWarnings = FALSE)
saveRDS(QCedRDS.combined,paste(saveDir,"saveRDS_obj/Integrated_",samplename,"_",splitby,".RDS",sep = ""))
return(QCedRDS.combined)
}