library(data.table)
library(Seurat)
library(ggplot2)
library(Matrix)
library(dplyr)

# Save data as sparse matrix 
counts <- fread("inputs/GSE210963_counts.csv.gz", data.table = FALSE)

rownames(counts) <- counts$V1
counts$V1 <- NULL

counts_sparse <- as(as.matrix(counts), "dgCMatrix")

saveRDS(counts_sparse, file = "inputs/counts_sparse.rds")

cells <- fread("inputs/GSE210963_metadata.csv.gz")

# Create seurat object 
	# nCount_RNA: total reads per cell
	# nFeature_RNA: total genes detected per cell 
seu <- CreateSeuratObject(
  counts = counts_sparse,
  project = "MDSC"
)

seu@meta.data$patient <- cells$patient

seu <- subset(seu, subset = nFeature_RNA >= 500)

seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

pdf("outputs/qc1.pdf")

	VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))

dev.off()

seu.list <- SplitObject(seu, split.by = "patient")

seu.list <- lapply(seu.list, function(x) {
  x <- NormalizeData(x, normalization.method = "LogNormalize")
  x <- FindVariableFeatures(
    x,
    selection.method = "vst",
    nfeatures = 2000
  )
  x
})

anchors <- FindIntegrationAnchors(
  object.list = seu.list,
  dims = 1:30
)

seu.integrated <- IntegrateData(
  anchorset = anchors,
  dims = 1:30
)

# saveRDS(seu.integrated, file = "inputs/integrated.rds")

seu.integrated <- readRDS("inputs/integrated.rds")

DefaultAssay(seu.integrated) <- "integrated"

seu.integrated <- ScaleData(seu.integrated)
seu.integrated <- RunPCA(seu.integrated, npcs = 30)

pdf("outputs/qc2.pdf")

	ElbowPlot(seu.integrated)

dev.off()

seu.integrated <- RunUMAP(seu.integrated, dims = 1:30)

seu.integrated <- FindNeighbors(seu.integrated, dims = 1:30)
seu.integrated <- FindClusters(
  seu.integrated,
  resolution = 0.5
)

pdf("outputs/umap.pdf")
	DimPlot(seu.integrated, reduction = "umap", label = TRUE)
	DimPlot(seu.integrated, reduction = "umap", group.by = "patient")
dev.off()

library(SingleR)
library(celldex)

BiocManager::install("celldex")






