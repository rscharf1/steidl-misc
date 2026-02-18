library(data.table)
library(Seurat)
library(ggplot2)
library(Matrix)
library(dplyr)

seu_H1 <- readRDS("outputs/H1_out_no_whitelist/alevin/processed.rds")
seu_H2 <- readRDS("outputs/H2_out_no_whitelist/alevin/processed.rds")

seu_H1$condition <- "preleukemic"
seu_H2$condition <- "leukemic"

objs <- list(seu_H1, seu_H2)

features <- SelectIntegrationFeatures(object.list = objs)

objs <- lapply(objs, function(x) {
  x <- ScaleData(x, features = features)
  x <- RunPCA(x, features = features)
})

anchors <- FindIntegrationAnchors(object.list = objs, anchor.features = features)
combined <- IntegrateData(anchorset = anchors)

DefaultAssay(combined) <- "integrated"

combined <- ScaleData(combined)
combined <- RunPCA(combined)
combined <- RunUMAP(combined, dims = 1:20)
combined <- FindNeighbors(combined, dims = 1:20)
combined <- FindClusters(combined, resolution = 0.5)

FeaturePlot(combined, features = c("Kit","Ly6a","Procr","Hlf","Gata2"))
VlnPlot(combined, features = c("Kit","Ly6a","Hlf"), group.by="seurat_clusters")
