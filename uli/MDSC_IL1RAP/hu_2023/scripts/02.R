library(dplyr)
library(Seurat)
library(data.table)
library(ggplot2)
library(readr)
library(tidyverse)
library(Matrix)

# Merge AML samples into single Seurat object
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

cell_ids <- sub("_expression_counts\\.csv\\.gz$", "", basename(aml))

merged <- merge(obj_list[[1]], y = obj_list[-1], add.cell.ids = cell_ids)

merged <- NormalizeData(merged)
merged <- FindVariableFeatures(merged)
merged <- ScaleData(merged)
merged <- RunPCA(merged)
merged <- FindNeighbors(merged, dims = 1:20)
merged <- FindClusters(merged)
merged <- RunUMAP(merged, dims = 1:20)

# saveRDS(merged, "outputs/merged_AML.rds")

# Analysis
merged <- readRDS("outputs/merged_AML.rds")

merged$celltype <- "Unknown"
merged$celltype[merged$seurat_clusters %in% c(4, 25, 27)] <- "B cells"
merged$celltype[merged$seurat_clusters %in% c(0, 1, 8, 15, 20)] <- "CD4 T cells"
merged$celltype[merged$seurat_clusters %in% c(2, 6, 7, 10, 17, 19, 27)] <- "CD8 T cells"
merged$celltype[merged$seurat_clusters %in% c()] <- "Monocytes"


pdf("out2.pdf")
	DimPlot(merged, reduction = "umap", label = TRUE, label.size = 5)
	VlnPlot(merged, features = "CD19", group.by = "seurat_clusters", pt.size = 0) # B cells: 4, 25, 27
	VlnPlot(merged, features = "CD3D", group.by = "seurat_clusters", pt.size = 0) # T cells: 0, 1, 2, 6, 7, 8, 9, 10, 15, 17, 19, 20, 21, 24, 27
	VlnPlot(merged, features = "CD3E", group.by = "seurat_clusters", pt.size = 0) # T cells: 0, 1, 2, 6, 7, 8, 9, 10, 15, 17, 19, 20, 21, 24, 27
	VlnPlot(merged, features = "CD3G", group.by = "seurat_clusters", pt.size = 0) # T cells: 
	VlnPlot(merged, features = "CD4", group.by = "seurat_clusters", pt.size = 0) # CD4: 0, 1, 8, 10, 11, 15, 16, 19, 20, 22, 26, 28
	VlnPlot(merged, features = "CD8A", group.by = "seurat_clusters", pt.size = 0) # CD8: 2, 6, 7, 10, 17, 19, 27
	VlnPlot(merged, features = "CD8B", group.by = "seurat_clusters", pt.size = 0) # CD8: 2, 6, 7, 10, 17, 19, 27
	VlnPlot(merged, features = "CD14", group.by = "seurat_clusters", pt.size = 0) # Monocytes: 3, 5, 11, 13, 16, 21, 23, 28
	VlnPlot(merged, features = "NCAM1", group.by = "seurat_clusters", pt.size = 0) # NK: 9
	VlnPlot(merged, features = "NCAM2", group.by = "seurat_clusters", pt.size = 0) # NK
	DimPlot(merged, group.by = "celltype", label = TRUE)
dev.off()

# In the other paper, MDSCs had ITGAM (CD11b) and CD33



intersect(
c(0, 1, 2, 6, 7, 8, 9, 10, 15, 17, 19, 20, 21, 24, 27), 
c(0, 1, 8, 10, 11, 15, 16, 19, 20, 22, 26, 28)
) %>% dput()





