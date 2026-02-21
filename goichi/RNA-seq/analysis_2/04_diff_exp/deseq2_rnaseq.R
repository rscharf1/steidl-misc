suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(data.table)
  library(apeglm)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(ggrepel)
  library(matrixStats)
  library(pheatmap)
})

counts_file <- "../03_align/counts/counts.txt"
meta_file   <- "sample_metadata.csv"

counts <- read.delim(
  counts_file,
  comment.char = "#",
  check.names = FALSE
)

# featureCounts columns:
# Geneid Chr Start End Strand Length sample1 sample2 ...
count_matrix <- counts %>%
  select(-Chr, -Start, -End, -Strand, -Length) %>%
  column_to_rownames("Geneid")

colnames(count_matrix) <- colnames(count_matrix) %>% 
  basename() %>% 
  gsub(".Aligned.sortedByCoord.out.bam", "", .)

sample_names <- fread(meta_file)

coldata <- as.data.frame(sample_names)
rownames(coldata) <- coldata$sample_id
coldata$sample_id <- NULL

coldata$condition <- factor(coldata$condition)

coldata$day <- factor(c(1,1,2,2,3,3))
coldata$treatment <- factor(c("control","treated",
                              "control","treated",
                              "control","treated"))

count_matrix <- count_matrix[, rownames(coldata)]

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = coldata,
  design = ~ 1
)

dds <- estimateSizeFactors(dds)
norm_counts <- counts(dds, normalized = TRUE)

ensg_ids <- gsub("\\..*$", "", rownames(norm_counts))

symbol_map <- mapIds(
  org.Hs.eg.db,
  keys = ensg_ids,
  keytype = "ENSEMBL",
  column = "SYMBOL",
  multiVals = "first"
)

TFs <- c("SPI1", "GATA1", "RUNX1", "CEBPA")

make_day_volcano <- function(treated, control, day_label) {

  log2fc <- log2((norm_counts[, treated] + 1) /
                 (norm_counts[, control] + 1))

  mean_expr <- rowMeans(norm_counts[, c(treated, control)])

  dt <- data.table(
    gene = rownames(norm_counts),
    log2FoldChange = log2fc,
    meanExpr = mean_expr
  )

  dt$gene_name <- symbol_map

  dt <- dt[meanExpr > 10 | gene_name %in% TFs]

  dt[, signal := log10(meanExpr + 1)]

  dt$sig <- FALSE
  dt[abs(log2FoldChange) > 2]$sig <- TRUE
  dt$sig %>% table

  dt[gene_name %in% TFs]

  dt$plot_group <- "NS"
  dt$plot_group[dt$sig == TRUE] <- "Significant"
  dt$plot_group[dt$gene_name %in% TFs] <- "TF"

  ggplot(dt, aes(x = log2FoldChange, y = signal, color = plot_group)) +
    geom_point(
      data = dt[plot_group != "TF"],
      alpha = 0.6,
      size = 1.2
    ) +
    geom_point(
      data = dt[plot_group == "TF"],
      size = 2.5
    ) +
    geom_vline(xintercept = c(-2, 2), linetype = "dashed") +
    geom_text_repel(
      data = dt[gene_name %in% TFs],
      aes(label = gene_name),
      color = "black",
      size = 3.5,
      box.padding = 0.3,
      point.padding = 0.2,
      segment.color = "grey50",
      max.overlaps = Inf
    ) +
    scale_color_manual(
      values = c(
        "NS" = "grey80",
        "Significant" = "#d73027",   # nice red
        "TF" = "#4575b4"             # blue
      )
    ) +
    theme_classic(base_size = 14) +
    labs(
      title = paste("Day", day_label),
      x = "log2 Fold Change",
      y = "log10 Mean Expression",
      color = NULL
    )
}

p1 <- make_day_volcano("C1T", "C1U", "1")
p2 <- make_day_volcano("C2T", "C2U", "2")
p3 <- make_day_volcano("C3T", "C3U", "3")

treated <- "C1T"
control <- "C1U"
day_label <- "1"

pdf("Volcano.pdf")

  print(p1)
  print(p2)
  print(p3)

dev.off()

#############
# Heatmap 

log_mat <- log2(norm_counts + 1)
keep <- rowMeans(log_mat) > 1
log_mat <- log_mat[keep, ]

vars <- rowVars(log_mat)
top_genes <- order(vars, decreasing = TRUE)[1:2000]
mat_top <- log_mat[top_genes, ]

annotation_col <- data.frame(
  Day = factor(c("1","1","2","2","3","3")),
  Treatment = factor(c("Control","Treated",
                       "Control","Treated",
                       "Control","Treated"))
)

rownames(annotation_col) <- colnames(mat_top)

pdf("Heat.pdf")

  pheatmap(
    mat_top,
    scale = "row",                 # Z-score per gene
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    clustering_method = "complete",
    annotation_col = annotation_col,
    show_rownames = FALSE,
    fontsize_col = 12,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
    border_color = NA
  )

dev.off()













