library(data.table)
library(Seurat)
library(ggplot2)
library(Matrix)
library(dplyr)
library(stringr)

# Save data as sparse matrix 
counts <- fread("inputs/GSE210963_counts.csv.gz", data.table = FALSE)

rownames(counts) <- counts$V1
counts$V1 <- NULL

counts_sparse <- as(as.matrix(counts), "dgCMatrix")

# saveRDS(counts_sparse, file = "inputs/counts_sparse.rds")
counts_sparse <- readRDS("inputs/counts_sparse.rds")

cells <- fread("inputs/GSE210963_metadata.csv.gz")

# Create seurat object 
	# nCount_RNA: total reads per cell
	# nFeature_RNA: total genes detected per cell 
seu <- CreateSeuratObject(
  counts = counts_sparse,
  project = "MDSC"
)

seu@meta.data$patient <- cells$patient
seu@meta.data$cell_type <- cells$cell_type

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

# DefaultAssay(seu.integrated) <- "integrated"

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

seu.integrated@meta.data$patient <- rownames(seu.integrated@meta.data) %>% 
  str_sub(., start = -1, end = -1) %>% 
  as.numeric()

DefaultAssay(seu.integrated) <- "RNA"

pdf("outputs/umap.pdf")
	DimPlot(seu.integrated, reduction = "umap", label = TRUE)

	DimPlot(seu.integrated, reduction = "umap", group.by = "patient")

  FeaturePlot(
    seu.integrated,
    features = "ITGAM"
  ) + ggtitle("ITGAM (CD11b)")

  FeaturePlot(
    seu.integrated,
    features = "CD33"
  ) + ggtitle("CD33")

  FeaturePlot(
    seu.integrated,
    features = "HLA-DRA"
  ) + ggtitle("HLA-DRA")

  FeaturePlot(
    seu.integrated,
    features = "HLA-DRB1"
  ) + ggtitle("HLA-DRB1")

  FeaturePlot(
    seu.integrated,
    features = "IL1RAP"
  ) + ggtitle("IL1RAP")

dev.off()

cells <- fread("inputs/GSE210963_metadata.csv.gz")
seu.integrated@meta.data$cell_type <- cells$cell_type

cancer_type <- c(
  ME = "Melanoma",
  BR = "Breast cancer",
  HN = "HNSCC"
)

seu.integrated@meta.data$cancer_type <- cancer_type[seu.integrated@meta.data$cell_type]

seu.integrated@meta.data$label <- paste0(seu.integrated@meta.data$patient, "_", seu.integrated@meta.data$cancer_type)

pdf("outputs/il1rap.pdf", width = 25, height = 15)

  p <- FeaturePlot(
    seu.integrated,
    features = "IL1RAP",
    split.by = "label",
    combine = FALSE
  ) # + ggtitle("IL1RAP")

  wrap_plots(p, ncol = 2)

dev.off()


pdf("outputs/violin.pdf")

  VlnPlot(
    seu.integrated,
    features = "IL1RAP"
  )

dev.off()




rownames(seu.integrated[["RNA"]])[
  grep("HLA", rownames(seu.integrated[["RNA"]]), ignore.case=TRUE)
]






# install.packages(
#   "https://cran.r-project.org/src/contrib/Archive/Seurat/Seurat_4.3.0.1.tar.gz",
#   repos = NULL,
#   type = "source",
#   lib = "/gs/gsfs0/home/rscharf/R/x86_64-pc-linux-gnu-library/4.4"
# )

# library(Seurat, lib.loc="/gs/gsfs0/home/rscharf/R/x86_64-pc-linux-gnu-library/4.4")



# packageVersion("Seurat")

# install.packages(
#   "https://cran.r-project.org/src/contrib/Archive/SeuratObject/SeuratObject_4.1.4.tar.gz",
#   repos = NULL,
#   type = "source",
#   lib = "/gs/gsfs0/home/rscharf/R/x86_64-pc-linux-gnu-library/4.4"
# )





