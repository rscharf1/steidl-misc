suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(data.table)
  library(apeglm)
  library(ggrepel)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(pheatmap)
  library(msigdbr) 
  library(patchwork)
  library(matrixStats)
})

counts_file <- "../04_counts/featureCounts_output/counts.txt"
meta_file   <- "sample_metadata.csv"

counts <- read.delim(
  counts_file,
  comment.char = "#",
  check.names = FALSE
)

sample_names <- fread(meta_file)

coldata <- as.data.frame(sample_names)
rownames(coldata) <- coldata$sample_id
coldata$sample_id <- NULL

coldata$protocol <- ifelse(
  grepl("^25444R-04", rownames(coldata)),
  "smartseq",
  "standard"
)

coldata$protocol <- factor(coldata$protocol)

coldata$condition <- factor(
  coldata$condition,
  levels = c("control", "AraC_1_day", "AraC_3_day")
)

# featureCounts columns:
# Geneid Chr Start End Strand Length sample1 sample2 ...
count_matrix <- counts %>%
  dplyr::select(-Chr, -Start, -End, -Strand, -Length) %>%
  column_to_rownames("Geneid")

colnames(count_matrix) <- colnames(count_matrix) %>% 
  basename() %>% 
  gsub(".Aligned.sortedByCoord.out.bam", "", .)

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData   = coldata,
  design    = ~ 1
)

dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- estimateSizeFactors(dds)

vsd <- vst(dds, blind = TRUE)
mat <- assay(vsd)

gene_var <- rowVars(mat)
summary(gene_var)

top_n <- 1000
top_genes <- names(sort(gene_var, decreasing = TRUE))[1:top_n]

mat_top <- mat[top_genes, ]
dim(mat_top)

mat_z <- t(scale(t(mat_top)))

annotation_col <- coldata[, c("condition", "protocol")]

annotation_col$condition <- factor(
  annotation_col$condition,
  levels = c("control", "AraC_1_day", "AraC_3_day")
)

annotation_col$protocol <- factor(
  annotation_col$protocol,
  levels = c("standard", "smartseq")
)

pdf("v4/heatmap.pdf")

  pheatmap(
    mat_z,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = FALSE,
    annotation_col = annotation_col,
    fontsize_col = 10,
    border_color = NA,
    main = "Global gene expression heatmap (VST, Z-score)"
  )

dev.off()









