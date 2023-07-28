scRNASeurat <- function(QCedRDS,Dims,saveDir, Assay = "RNA", samplename, condition)
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

message(paste("\nNormalizing\n"))
QCedRDS <- NormalizeData(QCedRDS, normalization.method="LogNormalize", scale.factor = 10000)

message(paste("\nFinding Variable Genes\n"))
QCedRDS <- FindVariableFeatures(QCedRDS, selection.method="vst", nfeatures = 2000)
all_genes_QCedRDS <- rownames(QCedRDS)

message(paste("\nScaling\n"))
QCedRDS <- ScaleData(Projec1_subset, features = all_genes_QCedRDS)

    
message(paste("\nRun PCA\n"))
QCedRDS <- RunPCA(QCedRDS, features = VariableFeatures(object=QCedRDS))
print(QCedRDS[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(QCedRDS, dims = 1:2, reduction = "pca")
DimPlot(QCedRDS, reduction = "pca")

message(paste("\nRun UMAP\n"))
QCedRDS <- RunUMAP(QCedRDS, dims = 1:Dims)

message(paste("\nShowing UMAP\n"))
print(DimPlot(QCedRDS, reduction = "umap"))

dir.create(paste(saveDir,"UMAP/",sep = ""), showWarnings = FALSE)
pdf(paste(saveDir,"UMAP/umap_",samplename,"_",Assay,"_",condition,".pdf",sep = ""), width = 7, height = 12)
print(DimPlot(QCedRDS, reduction = "umap"))
dev.off()
    
message(paste("\nFinding Neighbours for Clustering\n"))
QCedRDS <- FindNeighbors(QCedRDS, dims = 1:Dims)

resolution <- c(0.001,0.01,0.03,0.05,0.07,0.09,0.1,0.125,0.15,0.2,0.3,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2)
for(j in 1:length(resolution))
{
    QCedRDS <- FindClusters(QCedRDS, resolution = resolution[j], verbose = FALSE)
}
dir.create(paste(saveDir,"clustertree/",sep = ""), showWarnings = FALSE)
pdf(paste(saveDir,"clustertree/clustering_",samplename,"_",Assay,"_",condition,".pdf",sep = ""), width = 7, height = 12)
print(clustree(QCedRDS))
dev.off()
message("\n\nplease check for the resolution using clustertree is at this location: ",saveDir,"clustertree/\n\n")
return(QCedRDS)
}