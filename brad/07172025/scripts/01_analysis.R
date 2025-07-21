library(dplyr)
library(Seurat)
library(patchwork)
library(data.table)
library(ggplot2)
library(stringr)
library(BoneMarrowMap)

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

counts <- Read10X(data.dir = "inputs/GBA13/")

# Create Seurat object
seurat_obj <- CreateSeuratObject(counts = counts)

seurat_obj$batch <- "sample1"

query <- map_Query(
  exp_query = GetAssayData(seurat_obj, assay = "RNA", slot = "counts"),
  metadata_query = seurat_obj@meta.data,
  ref_obj = ref,             # should already be loaded from bonemarrowmap
  vars = "batch"             # or whatever metadata column you're using
)

query <- query %>% calculate_MappingError(., reference = ref, MAD_threshold = 2.5) 

# Get QC Plots
QC_plots <- plot_MappingErrorQC(query)

pdf("outputs/qc.pdf", width = 15)
	patchwork::wrap_plots(QC_plots, ncol = 4, widths = c(0.8, 0.3, 0.8, 0.3))
dev.off()

query <- subset(query, mapping_error_QC == 'Pass')

query <- predict_CellTypes(
  query_obj = query, 
  ref_obj = ref, 
  initial_label = 'initial_CellType', # celltype assignments before filtering on mapping QC
  final_label = 'predicted_CellType'  # celltype assignments with map QC failing cells assigned as NA
) 

pdf("outputs/proj_umap.pdf", width = 15)
	DimPlot(subset(query, mapping_error_QC == 'Pass'), reduction = 'umap_projected', group.by = c('predicted_CellType'), 
	        raster=FALSE, label=TRUE, label.size = 4)
	DimPlot(subset(query, mapping_error_QC == 'Pass'), reduction = 'umap_projected', group.by = c('predicted_CellType_Broad'), 
	        raster=FALSE, label=TRUE, label.size = 4)
	DimPlot(subset(query, mapping_error_QC == 'Pass'), reduction = 'umap_projected', group.by = c('predicted_CellType_Broad'), 
	        raster=FALSE)
dev.off()

