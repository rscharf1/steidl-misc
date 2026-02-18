library(Seurat)
library(patchwork)
library(dplyr)
library(ggplot2)

# PY-Hi

data <- Read10X_h5("../cellranger_count/PY-Hi/filtered_feature_bc_matrix.h5")

seu <- CreateSeuratObject(
  counts = data,
  project = "scRNA",
  min.cells = 3,
  min.features = 200
)

seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

pdf("PY-Hi-QC.pdf")

	VlnPlot(
	  seu,
	  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
	  ncol = 3
	)

	FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "percent.mt")

	seu <- subset(
	  seu,
	  subset =
	    nFeature_RNA > 500 &
	    nFeature_RNA < 7000 &
	    percent.mt < 20
	)

	VlnPlot(
	  seu,
	  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
	  ncol = 3
	)

	FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "percent.mt")

dev.off()


#########################
# PY-Lo
data <- Read10X_h5("../cellranger_count/PY-Lo/filtered_feature_bc_matrix.h5")

seu <- CreateSeuratObject(
  counts = data,
  project = "scRNA",
  min.cells = 3,
  min.features = 200
)

seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

pdf("PY-Lo-QC.pdf")

	VlnPlot(
	  seu,
	  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
	  ncol = 3
	)

	FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "percent.mt")

	seu <- subset(
	  seu,
	  subset =
	    nFeature_RNA > 500 &
	    nFeature_RNA < 7000 &
	    percent.mt < 20
	)

	VlnPlot(
	  seu,
	  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
	  ncol = 3
	)

	FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "percent.mt")

dev.off()

#########################
# Combine 

pyhi <- Read10X_h5("../cellranger_count/PY-Hi/filtered_feature_bc_matrix.h5")
pylo <- Read10X_h5("../cellranger_count/PY-Lo/filtered_feature_bc_matrix.h5")

pyhi <- CreateSeuratObject(
  counts = pyhi,
  project = "scRNA",
  min.cells = 3,
  min.features = 200
)

pylo <- CreateSeuratObject(
  counts = pylo,
  project = "scRNA",
  min.cells = 3,
  min.features = 200
)

pyhi[["percent.mt"]] <- PercentageFeatureSet(pyhi, pattern = "^MT-")
pylo[["percent.mt"]] <- PercentageFeatureSet(pylo, pattern = "^MT-")

pyhi$PY_status <- "PY_hi"
pylo$PY_status <- "PY_lo"

pyhi <- subset(
  pyhi,
  subset =
    nFeature_RNA > 500 &
    nFeature_RNA < 7000 &
    percent.mt < 20
)

pylo <- subset(
  pylo,
  subset =
    nFeature_RNA > 500 &
    nFeature_RNA < 7000 &
    percent.mt < 20
)

seu <- merge(pyhi, y = pylo)

table(seu$PY_status)

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)
seu <- ScaleData(seu)
seu <- RunPCA(seu)
seu <- RunUMAP(seu, dims = 1:30)

pdf("UMAP.pdf")

	DimPlot(seu, group.by = "PY_status")

dev.off()

tfs <- c("GATA1", "SPI1", "RUNX1", "CEBPA", "FOXO4", "MYC")

tfs %in% rownames(seu)

pdf("Compare_TFs_UMAP.pdf", width = 12, height = 20)

	FeaturePlot(
	  seu,
	  features = tfs,
	  split.by = "PY_status",
	  max.cutoff = "q95"
	)

dev.off()

pdf("Compare_TFs_Violin.pdf")

	VlnPlot(
	  seu,
	  features = tfs,
	  group.by = "PY_status",
	  pt.size = 0
	)

dev.off()

#####################
# Unique states 

seu <- FindNeighbors(seu, dims = 1:30)
seu <- FindClusters(seu, resolution = 0.4)

table(seu$seurat_clusters, seu$PY_status)

prop.table(table(seu$seurat_clusters, seu$PY_status), 1)

#####################

data_hi <- Read10X_h5("../cellranger_count/PY-Hi/filtered_feature_bc_matrix.h5")

seu_hi <- CreateSeuratObject(
  counts = data_hi,
  project = "PY_Hi",
  min.cells = 3,
  min.features = 200
)

seu_hi[["percent.mt"]] <- PercentageFeatureSet(seu_hi, pattern = "^MT-")

seu_hi <- subset(
  seu_hi,
  subset =
    nFeature_RNA > 500 &
    nFeature_RNA < 7000 &
    percent.mt < 20
)

seu_hi$PY_status <- "PY_hi"


data_lo <- Read10X_h5("../cellranger_count/PY-Lo/filtered_feature_bc_matrix.h5")

seu_lo <- CreateSeuratObject(
  counts = data_lo,
  project = "PY_Lo",
  min.cells = 3,
  min.features = 200
)

seu_lo[["percent.mt"]] <- PercentageFeatureSet(seu_lo, pattern = "^MT-")

seu_lo <- subset(
  seu_lo,
  subset =
    nFeature_RNA > 500 &
    nFeature_RNA < 7000 &
    percent.mt < 20
)

seu_lo$PY_status <- "PY_lo"


df_hi <- data.frame(
  nFeature_RNA = seu_hi$nFeature_RNA,
  nCount_RNA   = seu_hi$nCount_RNA,
  PY_status    = seu_hi$PY_status
)

df_lo <- data.frame(
  nFeature_RNA = seu_lo$nFeature_RNA,
  nCount_RNA   = seu_lo$nCount_RNA,
  PY_status    = seu_lo$PY_status
)

qc_df <- rbind(df_hi, df_lo)

pdf("Pyronin_comparison.pdf")

	ggplot(qc_df, aes(x = PY_status, y = nFeature_RNA, fill = PY_status)) +
	  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
	  coord_cartesian(ylim = c(500, 7000)) +
	  theme_classic() +
	  labs(y = "Detected genes per cell", x = NULL)

	ggplot(qc_df, aes(x = PY_status, y = nCount_RNA, fill = PY_status)) +
	  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
	  theme_classic() +
	  labs(y = "UMI counts per cell", x = NULL)

	ggplot(qc_df, aes(x = PY_status, y = nCount_RNA, fill = PY_status)) +
	  geom_violin(trim = FALSE, alpha = 0.6) +
	  geom_jitter(
	    width = 0.15,
	    size = 0.3,
	    alpha = 0.3,
	    color = "black"
	  ) +
	  # scale_y_log10() +
	  theme_classic() +
	  labs(
	    x = NULL,
	    y = "UMI counts per cell"
	  ) +
	  theme(
	    legend.position = "none"
	  )

dev.off()











