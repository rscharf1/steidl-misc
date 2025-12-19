library(data.table)
library(Seurat)
library(ggplot2)
library(Matrix)
library(dplyr)

# PROCESS SAMPLE 1

# No whitelist first 
counts <- readMM("outputs/H1_out_no_whitelist/alevin/quants_mat.mtx.gz")
rows <- fread("outputs/H1_out_no_whitelist/alevin/quants_mat_rows.txt", header = FALSE)

cols <- fread("outputs/H1_out_no_whitelist/alevin/quants_mat_cols.txt", header = FALSE)
gene_names <- fread("inputs/txp2gene_vM25.tsv", header = FALSE)

gene_table <- merge(cols, unique(gene_names[, c(2,3)]), by.x="V1", by.y="V2", all.y=FALSE)

gene_table <- gene_table[cols, ]

gene_table[, V3 := make.unique(V3)]

rownames(counts) <- rows$V1
colnames(counts) <- gene_table$V3

counts <- t(counts)

# saveRDS(counts, file = "outputs/H1_out_no_whitelist/alevin/counts_sparse.rds")

counts <- readRDS("outputs/H1_out_no_whitelist/alevin/counts_sparse.rds")

seu <- CreateSeuratObject(
  counts = counts
)

seu <- subset(seu, subset = nFeature_RNA >= 500)

pdf("outputs/qc1.pdf")
	VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA"))
dev.off()

seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000)

seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)

seu <- ScaleData(seu, features = rownames(seu))

seu <- RunPCA(seu, features = VariableFeatures(object = seu))

seu <- FindNeighbors(seu, dims = 1:10)
seu <- FindClusters(seu, resolution = 0.5)

seu <- RunUMAP(seu, dims = 1:10)

# saveRDS(seu, file = "outputs/H1_out_no_whitelist/alevin/processed.rds")

seu <- readRDS("outputs/H1_out_no_whitelist/alevin/processed.rds")

pdf("outputs/umap.pdf")
	DimPlot(seu, reduction = "umap")
	VlnPlot(seu, features = c("Cnot6l"))
	FeaturePlot(seu, features = c("Cnot6l"))
dev.off()

########################################
# PROCESS SAMPLE 2

# No whitelist first 
counts <- readMM("outputs/H2_out_no_whitelist/alevin/quants_mat.mtx.gz")
rows <- fread("outputs/H2_out_no_whitelist/alevin/quants_mat_rows.txt", header = FALSE)

cols <- fread("outputs/H2_out_no_whitelist/alevin/quants_mat_cols.txt", header = FALSE)
gene_names <- fread("inputs/txp2gene_vM25.tsv", header = FALSE)

gene_table <- merge(cols, unique(gene_names[, c(2,3)]), by.x="V1", by.y="V2", all.y=FALSE)

gene_table <- gene_table[cols, ]

gene_table[, V3 := make.unique(V3)]

rownames(counts) <- rows$V1
colnames(counts) <- gene_table$V3

counts <- t(counts)

# saveRDS(counts, file = "outputs/H2_out_no_whitelist/alevin/counts_sparse.rds")

# counts <- readRDS("outputs/H2_out_no_whitelist/alevin/counts_sparse.rds")

seu <- CreateSeuratObject(
  counts = counts
)

seu <- subset(seu, subset = nFeature_RNA >= 500)

pdf("outputs/H2_qc1.pdf")
	VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA"))
dev.off()

seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000)

seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)

seu <- ScaleData(seu, features = rownames(seu))

seu <- RunPCA(seu, features = VariableFeatures(object = seu))

seu <- FindNeighbors(seu, dims = 1:10)
seu <- FindClusters(seu, resolution = 0.5)

seu <- RunUMAP(seu, dims = 1:10)

# saveRDS(seu, file = "outputs/H2_out_no_whitelist/alevin/processed.rds")

seu <- readRDS("outputs/H2_out_no_whitelist/alevin/processed.rds")

pdf("outputs/H2_umap.pdf")
	DimPlot(seu, reduction = "umap")
	VlnPlot(seu, features = c("Cnot6l"))
	FeaturePlot(seu, features = c("Cnot6l"))
dev.off()




