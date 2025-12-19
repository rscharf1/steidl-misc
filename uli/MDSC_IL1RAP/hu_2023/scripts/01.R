library(dplyr)
library(Seurat)
library(data.table)
library(ggplot2)
library(readr)
library(tidyverse)
library(Matrix)

# Load and create Seurat object 
matrix_data <- fread("inputs/AML1_expression_counts.csv.gz")

colnames(matrix_data)[1] <- "gene"

matrix_data <- as.data.frame(matrix_data)  # ensure data.frame
rownames(matrix_data) <- matrix_data$gene
matrix_data$gene <- NULL

matrix_data <- as.matrix(matrix_data)
storage.mode(matrix_data) <- "double"

sparse_mat <- Matrix(matrix_data, sparse = TRUE)

obj <- CreateSeuratObject(counts = sparse_mat)

# Create UMAP 
obj <- NormalizeData(obj)

obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)

obj <- ScaleData(obj)

obj <- RunPCA(obj, features = VariableFeatures(obj))

obj <- FindNeighbors(obj, dims = 1:20)

obj <- FindClusters(obj, resolution = 0.5)

obj <- RunUMAP(obj, dims = 1:20)

genes <- Features(obj[["RNA"]])
genes[grep("NCAM", genes)]

obj$celltype <- "Unknown"

obj$celltype[obj$seurat_clusters == 8] <- "B cells"
obj$celltype[obj$seurat_clusters == 6] <- "CD8 T cells"
obj$celltype[obj$seurat_clusters %in% c(2,3)] <- "CD4 T cells"
obj$celltype[obj$seurat_clusters == 5] <- "Monocytes"

pdf("out.pdf")
	DimPlot(obj, reduction = "umap", label = TRUE, label.size = 5)
	VlnPlot(obj, features = "CD19", group.by = "seurat_clusters", pt.size = 0) # B cells: 8
	VlnPlot(obj, features = "CD3D", group.by = "seurat_clusters", pt.size = 0) # T cells: 2, 3, 6
	VlnPlot(obj, features = "CD3E", group.by = "seurat_clusters", pt.size = 0) # T cells: 2, 3, 6
	VlnPlot(obj, features = "CD3G", group.by = "seurat_clusters", pt.size = 0) # T cells: 2, 3, 6
	VlnPlot(obj, features = "CD4", group.by = "seurat_clusters", pt.size = 0) # CD4: 2, 3
	VlnPlot(obj, features = "CD8A", group.by = "seurat_clusters", pt.size = 0) # CD8: 6
	VlnPlot(obj, features = "CD8B", group.by = "seurat_clusters", pt.size = 0) # CD8: 3, 6
	VlnPlot(obj, features = "CD14", group.by = "seurat_clusters", pt.size = 0) # Monocytes: 5
	VlnPlot(obj, features = "NCAM1", group.by = "seurat_clusters", pt.size = 0) # NK
	VlnPlot(obj, features = "NCAM2", group.by = "seurat_clusters", pt.size = 0) # NK
	DimPlot(obj, group.by = "celltype", label = TRUE)
dev.off()





