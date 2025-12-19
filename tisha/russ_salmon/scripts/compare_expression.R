library(data.table)
library(Seurat)
library(ggplot2)
library(Matrix)
library(dplyr)

seu_H1 <- readRDS("outputs/H1_out_no_whitelist/alevin/processed.rds")
seu_H2 <- readRDS("outputs/H2_out_no_whitelist/alevin/processed.rds")

seu_merged <- merge(
  seu_H1,
  y = seu_H2,
  add.cell.ids = c("H1", "H2"),
  project = "H1_vs_H2"
)

seu_merged <- NormalizeData(seu_merged)
seu_merged <- FindVariableFeatures(seu_merged)
seu_merged <- ScaleData(seu_merged)
seu_merged <- RunPCA(seu_merged)
seu_merged <- FindNeighbors(seu_merged, dims = 1:10)
seu_merged <- FindClusters(seu_merged, resolution = 0.5)
seu_merged <- RunUMAP(seu_merged, dims = 1:10)

# saveRDS(seu_merged, file = "outputs/no_whitelist_merged.rds")

seu_merged <- readRDS("outputs/no_whitelist_merged.rds")

seu_merged$sample <- ifelse(
	grepl("H1", rownames(seu_merged@meta.data)),
	"H1", "H2"
)

pdf("outputs/compare1.pdf")

	DimPlot(seu_merged, reduction = "umap")

	DimPlot(
	  seu_merged,
	  reduction = "umap",
	  group.by = "sample"
	)

	FeaturePlot(
	  seu_merged,
	  features = "Cnot6l",
	  split.by = "sample"
	)

	FeaturePlot(
	  subset(seu_merged, sample == "H1"),
	  features = "Cnot6l"
	) + ggtitle("Cnot6l – H1")

	FeaturePlot(
	  subset(seu_merged, sample == "H2"),
	  features = "Cnot6l"
	) + ggtitle("Cnot6l – H2")

	VlnPlot(
	  seu_merged,
	  features = "Cnot6l",
	  group.by = "seurat_clusters",
	  split.by = "sample",
	  pt.size = 0
	)

dev.off()

# Positive control
pb <- AggregateExpression(
  seu_merged,
  group.by = "sample",
  assays = "RNA",
  slot = "counts"
)$RNA

expr_filter <- rowSums(pb > 10) == 2

pb <- pb[expr_filter, ]

log_fc <- log1p(pb[, "H2"]) - log1p(pb[, "H1"])

log_fc <- sort(log_fc, decreasing = TRUE)

top_genes <- c(
	names(head(log_fc, 2)),
	names(tail(log_fc, 2))
)

pdf("outputs/compare2.pdf", width = 15, height = 15)

	FeaturePlot(
	  seu_merged,
	  features = top_genes,
	  split.by = "sample",
	  ncol = length(top_genes)
	)

	VlnPlot(
	  seu_merged,
	  features = top_genes,
	  group.by = "seurat_clusters",
	  split.by = "sample",
	  pt.size = 0
	)

dev.off()

# Negative control

pb <- AggregateExpression(
  seu_merged,
  group.by = "sample",
  assays = "RNA",
  slot = "counts"
)$RNA

expr_filter <- rowSums(pb > 2000) == 2

pb <- pb[expr_filter, ]

log_fc <- log1p(pb[, "H2"]) - log1p(pb[, "H1"])

top_genes <- log_fc %>% abs() %>% sort() %>% head(2) %>% names()

pb[log_fc %>% abs() %>% sort() %>% head(20) %>% names(), ]

pdf("outputs/compare3.pdf", width = 10, height = 7)

	FeaturePlot(
	  seu_merged,
	  features = top_genes,
	  split.by = "sample",
	  ncol = length(top_genes)
	)

	VlnPlot(
	  seu_merged,
	  features = top_genes,
	  group.by = "seurat_clusters",
	  split.by = "sample",
	  pt.size = 0
	)

dev.off()

# similar_genes <- log_fc[log_fc > -0.001 & log_fc < 0.001]

# expr_filter <- rowSums(pb > 1000) == 2

# pb[names(similar_genes), ][]

















