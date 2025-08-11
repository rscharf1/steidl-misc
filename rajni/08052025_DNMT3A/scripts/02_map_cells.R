library(Seurat)
library(data.table)
library(dplyr)
library(BoneMarrowMap)
library(Matrix)

# Plan
	# Send each sample through the BMM 
	# Label cells that are bonified LT-HSCs -> Save intermediate file
	# Subset to include only LT-HSCs 
	# Look at expression of genes of interest 


# Set up map 
projection_path = '~/Tools/references/bonemarrowmap/'

ref <- readRDS(paste0(projection_path, 'BoneMarrowMap_SymphonyReference.rds'))
ref$save_uwot_path <- paste0(projection_path, 'BoneMarrowMap_uwot_model.uwot')

ReferenceSeuratObj <- create_ReferenceObject(ref)

ReferenceSeuratObj$CellType_Annotation_formatted <- ReferenceSeuratObj$CellType_Annotation_formatted %>% 
	gsub("\t\n|\t|\n", " ", .)

pdf("outputs/ref.pdf", width = 25, height = 15)
	DimPlot(ReferenceSeuratObj, reduction = 'umap', group.by = 'CellType_Annotation_formatted', 
	        raster=FALSE, label=TRUE, label.size = 4) + NoAxes()
	DimPlot(ReferenceSeuratObj, reduction = 'umap', group.by = 'CellType_Broad', 
	        raster=FALSE, label=TRUE, label.size = 4)
	FeaturePlot(ReferenceSeuratObj, reduction = 'umap', features = 'Pseudotime', raster=FALSE) + NoAxes()
dev.off()

# Single sample through BMM 
d <- "inputs/GSE168807_scRNA_mtx/GSM5171988"

m <- ReadMtx(
  mtx = file.path(d, "GSM5171988_PBS_VM_matrix.mtx.gz"),
  features = file.path(d, "GSM5171988_PBS_VM_features.tsv.gz"),
  cells = file.path(d, "GSM5171988_PBS_VM_barcodes.tsv.gz"),
  feature.column = 2   # use gene symbols (col 2 of features.tsv.gz)
)

rownames(m) <- toupper(rownames(m))

obj <- CreateSeuratObject(
  counts = m,
  project = "GSM5171988",
  min.cells = 3,
  min.features = 200
)

obj$batch <- "GSM5171988"

query <- map_Query(
  exp_query = GetAssayData(obj, assay = "RNA", slot = "counts"),
  metadata_query = obj@meta.data,
  ref_obj = ref,             # should already be loaded from bonemarrowmap
  vars = "batch"             # or whatever metadata column you're using
)

# Run QC based on mapping error score, flag cells with mapping error >= 2.5 MADs above median
query <- query %>% calculate_MappingError(., reference = ref, MAD_threshold = 2.5) 

# Get QC Plots
QC_plots <- plot_MappingErrorQC(query)

pdf("outputs/sample_qc.pdf", width = 15)
	patchwork::wrap_plots(QC_plots, ncol = 4, widths = c(0.8, 0.3, 0.8, 0.3))
dev.off()

query <- subset(query, mapping_error_QC == 'Pass')

query <- predict_CellTypes(
  query_obj = query, 
  ref_obj = ref, 
  initial_label = 'initial_CellType', # celltype assignments before filtering on mapping QC
  final_label = 'predicted_CellType'  # celltype assignments with map QC failing cells assigned as NA
) 

pdf("outputs/sample_umap.pdf", width = 15)
	DimPlot(subset(query, mapping_error_QC == 'Pass'), reduction = 'umap_projected', group.by = c('predicted_CellType'), 
	        raster=FALSE, label=TRUE, label.size = 4)
	DimPlot(subset(query, mapping_error_QC == 'Pass'), reduction = 'umap_projected', group.by = c('predicted_CellType_Broad'), 
	        raster=FALSE, label=TRUE, label.size = 4)
dev.off()

# Compile all samples into single seurat object 
samples <- list.dirs("inputs/GSE168807_scRNA_mtx", full.names = TRUE, recursive = FALSE)

obj_list <- lapply(seq_along(samples), function(sample) {
	message(paste0(sample, " / ", length(samples)))

	files <- list.files(samples[sample], full.names = TRUE, recursive = FALSE)

	m <- ReadMtx(
		cells = files[grepl("barcode", files, ignore.case=TRUE)],
		features = files[grepl("feature", files, ignore.case=TRUE)],
		mtx = files[grepl("matrix", files, ignore.case=TRUE)],
		feature.column = 2   # use gene symbols (col 2 of features.tsv.gz)
	)

	rownames(m) <- toupper(rownames(m))

	obj <- CreateSeuratObject(
	  counts = m,
	  project = basename(samples[sample]),
	  min.cells = 3,
	  min.features = 200
	)

	obj$batch <- basename(samples[sample])

	obj
})

obj_merged <- Reduce(function(x, y) merge(x, y), obj_list)

rna <- obj_merged[["RNA"]]
lays <- Layers(rna)

# Union of all genes in the merged assay
union_genes <- unique(rownames(rna))

get_layer_aligned <- function(layer) {
  m <- LayerData(rna, layer = layer)
  if (!inherits(m, "dgCMatrix")) {
    m <- as(Matrix(m, sparse = TRUE), "dgCMatrix")
  }
  # Add zero rows for genes missing from this layer
  missing <- setdiff(union_genes, rownames(m))
  if (length(missing) > 0) {
    zeros <- Matrix(0, nrow = length(missing), ncol = ncol(m), sparse = TRUE)
    rownames(zeros) <- missing
    m <- rbind(m, zeros)
  }
  # Reorder to the union order
  m <- m[union_genes, , drop = FALSE]
  m
}

# Align all layers, then column-bind
mats <- lapply(lays, get_layer_aligned)
cts  <- do.call(cbind, mats)

# Replace RNA assay with a single counts layer
obj_merged[["RNA"]] <- CreateAssay5Object(counts = cts)
DefaultAssay(obj_merged) <- "RNA"

# sanity check: should now show a single "counts" layer
Layers(obj_merged[["RNA"]])

counts_q <- LayerData(obj_merged[["RNA"]], layer = "counts")

query <- map_Query(
  exp_query = counts_q,
  metadata_query = obj_merged@meta.data,
  ref_obj = ref,
  vars = "batch"
)

# Run QC based on mapping error score, flag cells with mapping error >= 2.5 MADs above median
query <- query %>% calculate_MappingError(., reference = ref, MAD_threshold = 2.5) 

# Get QC Plots
QC_plots <- plot_MappingErrorQC(query)

pdf("outputs/combined_qc.pdf", width = 15)
	patchwork::wrap_plots(QC_plots, ncol = 4, widths = c(0.8, 0.3, 0.8, 0.3))
dev.off()

query <- subset(query, mapping_error_QC == 'Pass')

query <- predict_CellTypes(
  query_obj = query, 
  ref_obj = ref, 
  initial_label = 'initial_CellType', # celltype assignments before filtering on mapping QC
  final_label = 'predicted_CellType'  # celltype assignments with map QC failing cells assigned as NA
) 

pdf("outputs/combined_umap.pdf", width = 15)
	DimPlot(subset(query, mapping_error_QC == 'Pass'), reduction = 'umap_projected', group.by = c('predicted_CellType'), 
	        raster=FALSE, label=TRUE, label.size = 4)
	DimPlot(subset(query, mapping_error_QC == 'Pass'), reduction = 'umap_projected', group.by = c('predicted_CellType_Broad'), 
	        raster=FALSE, label=TRUE, label.size = 4)
dev.off()

saveRDS(query, "outputs/combined_mapped.rds")









