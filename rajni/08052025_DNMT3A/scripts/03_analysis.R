library(Seurat)
library(data.table)
library(dplyr)
library(BoneMarrowMap)
library(Matrix)

# Now I have a RDS object that is mapped to cell types 
# Only really care about the HSCs and how the expression of these 5 genes looks in each of the samples 

obj <- readRDS("outputs/combined_mapped.rds")


obj@meta.data$predicted_CellType_No_Abbreviations %>% unique()

hsc <- subset(
  obj, 
  subset = predicted_CellType_No_Abbreviations == "Hematopoietic Stem Cell"
)

# Normalize + UMAP
hsc <- NormalizeData(hsc)                                   # creates "data" layer in v5
hsc <- FindVariableFeatures(hsc, nfeatures = 3000)

# PCA → neighbors → UMAP → clusters
hsc <- RunPCA(hsc, features = VariableFeatures(hsc))
hsc <- FindNeighbors(hsc, dims = 1:20)
hsc <- FindClusters(hsc, resolution = 0.2)
hsc <- RunUMAP(hsc, dims = 1:20)

pdf("outputs/HSC_umap.pdf")
	DimPlot(hsc, group.by = "batch", reduction = "umap")
	FeaturePlot(hsc, features = c("IRF8","RUNX1"))
dev.off()

# One sample
sample_id <- "GSM5171988"

# 1) subset to just that sample (assumes metadata column 'batch')
s <- subset(hsc, subset = batch == sample_id)

# 2) make sure we have a normalized layer for plotting
if (!"data" %in% Layers(s[["RNA"]])) s <- NormalizeData(s)

# 3) side-by-side UMAP feature plots for Irf8 and Runx1
p <- FeaturePlot(
  s,
  features = c("Irf8","Runx1"),
  reduction = "umap",
  layer = "data",        # Seurat v5: plot normalized data layer
  order = TRUE,
  ncol = 2
)
p
