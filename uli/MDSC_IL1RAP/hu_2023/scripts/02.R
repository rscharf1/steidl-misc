library(dplyr)
library(Seurat)
library(data.table)
library(ggplot2)
library(readr)
library(tidyverse)
library(Matrix)

files <- list.files("inputs", pattern = "gz$", full.names = TRUE)

aml <- files[grep("AML", files)]

load_seurat_object <- function(path) {
	message("Loading: ", path)

	matrix_data <- fread(path)

	colnames(matrix_data)[1] <- "gene"

	matrix_data <- as.data.frame(matrix_data)  # ensure data.frame
	rownames(matrix_data) <- matrix_data$gene
	matrix_data$gene <- NULL

	matrix_data <- as.matrix(matrix_data)
	storage.mode(matrix_data) <- "double"

	sparse_mat <- Matrix(matrix_data, sparse = TRUE)

	obj <- CreateSeuratObject(counts = sparse_mat)

	obj
}

obj_list <- lapply(aml, load_seurat_object)

merged <- merge(obj_list[[1]], y = obj_list[-1], add.cell.ids = basename(aml[1:2]))

merged <- NormalizeData(merged)
merged <- FindVariableFeatures(merged)
merged <- ScaleData(merged)
merged <- RunPCA(merged)
merged <- FindNeighbors(merged, dims = 1:20)
merged <- FindClusters(merged)
merged <- RunUMAP(merged, dims = 1:20)

merged$celltype <- "Unknown"
merged$celltype[merged$seurat_clusters == 18] <- "B cells"
merged$celltype[merged$seurat_clusters %in% c(4, 6, 11, 14)] <- "CD4 T cells"
merged$celltype[merged$seurat_clusters %in% c(7, 8)] <- "CD8 T cells"
merged$celltype[merged$seurat_clusters %in% c(0, 1, 2)] <- "Monocytes"

pdf("out2.pdf")
	DimPlot(merged, reduction = "umap", label = TRUE, label.size = 5)
	VlnPlot(merged, features = "CD19", group.by = "seurat_clusters", pt.size = 0) # B cells: 16
	VlnPlot(merged, features = "CD3D", group.by = "seurat_clusters", pt.size = 0) # T cells: 4, 6, 7, 8, 12, 14
	VlnPlot(merged, features = "CD3E", group.by = "seurat_clusters", pt.size = 0) # T cells: 4, 6, 7, 8, 12, 14
	VlnPlot(merged, features = "CD3G", group.by = "seurat_clusters", pt.size = 0) # T cells: 4, 6, 7, 8, 12, 14
	VlnPlot(merged, features = "CD4", group.by = "seurat_clusters", pt.size = 0) # CD4: 4, 6, 11, 14
	VlnPlot(merged, features = "CD8A", group.by = "seurat_clusters", pt.size = 0) # CD8: 7, 8
	VlnPlot(merged, features = "CD8B", group.by = "seurat_clusters", pt.size = 0) # CD8: 7, 8
	VlnPlot(merged, features = "CD14", group.by = "seurat_clusters", pt.size = 0) # Monocytes: 0, 1, 2
	VlnPlot(merged, features = "NCAM1", group.by = "seurat_clusters", pt.size = 0) # NK
	VlnPlot(merged, features = "NCAM2", group.by = "seurat_clusters", pt.size = 0) # NK
	DimPlot(merged, group.by = "celltype", label = TRUE)
dev.off()

# library(SeuratWrappers)

# obj_list <- PrepSCTIntegration(obj_list)
# features <- SelectIntegrationFeatures(obj_list, nfeatures = 3000)
# obj_list <- SCTNormalize(obj_list)
# obj_list <- RunPCA(obj_list)

# anchors <- FindIntegrationAnchors(obj_list, dims = 1:30, features = features)
# merged <- IntegrateData(anchorset = anchors, dims = 1:30)

# merged <- RunPCA(merged)
# merged <- RunUMAP(merged, dims = 1:30)
# merged <- FindClusters(merged, resolution = 0.5)