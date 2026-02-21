library(dplyr)
library(data.table)
library(DESeq2)
library(ggplot2)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(GenomicRanges)
library(ggrepel)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(msigdbr)
library(pheatmap)

# Assemble counts matrix 
counts_dir <- "../03_count_mat/counts"
peaks_bed  <- "../03_count_mat/merged_peaks_500bp.filtered.sorted.bed"

count_files <- list.files(
  counts_dir,
  full.names = TRUE
)

peaks <- read.table(peaks_bed, header = FALSE, stringsAsFactors = FALSE)
colnames(peaks) <- c("chr", "start", "end")
peak_ids <- paste(peaks$chr, peaks$start, peaks$end, sep = ":")

counts_list <- lapply(count_files, function(f) {
  df <- read.table(f, header = FALSE)
  
  # Expect: chr start end count
  stopifnot(ncol(df) >= 4)
  
  df[[4]]
})

count_matrix <- do.call(cbind, counts_list)
rownames(count_matrix) <- peak_ids

colnames(count_matrix) <- gsub(
  "\\.counts\\.txt$",
  "",
  basename(count_files)
)

sample_info <- tibble(
  sample = colnames(count_matrix),
  condition = case_when(
    grepl("01_|02_", sample) ~ "control",
    grepl("03_|04_", sample) ~ "day_1",
    grepl("05_|06_", sample) ~ "day_14",
    TRUE ~ NA_character_
  )
)

# Create DESeq2 object 

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData   = sample_info,
  design    = ~ condition
)

dds <- dds[rowSums(counts(dds)) >= 10, ]

vsd <- vst(dds, blind = TRUE)

dds <- DESeq(dds)

res_day1_vs_ctrl <- results(dds, contrast = c("condition", "day_1", "control"))
res_day14_vs_ctrl <- results(dds, contrast = c("condition", "day_14", "control"))
res_day14_vs_day1 <- results(dds, contrast = c("condition", "day_14", "day_1"))

########
link_peaks_to_genes <- function(dt) {
  res_dt <- as.data.table(as.data.frame(dt), keep.rownames = "peak_id")
  
  coords <- tstrsplit(res_dt$peak_id, ":", fixed = TRUE)

  res_dt[, `:=`(
    chr   = coords[[1]],
    start = as.integer(coords[[2]]),
    end   = as.integer(coords[[3]])
  )]

  peaks <- readPeakFile("../03_count_mat/merged_peaks_500bp.filtered.sorted.bed")

  peak_anno <- annotatePeak(
    peaks,
    TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
    annoDb = "org.Hs.eg.db",
    tssRegion = c(-2000, 2000)
  )

  anno_df <- as.data.frame(peak_anno)

  anno_df$peak_id <- paste(
    anno_df$seqnames,
    anno_df$start,
    anno_df$end,
    sep = ":"
  )

  peak_dt <- as.data.table(anno_df)

  setnames(peak_dt,
           old = c("seqnames"),
           new = c("chr"))

  # Ensure integer coordinates
  peak_dt[, `:=`(
    start = as.integer(start),
    end   = as.integer(end)
  )]

  setkey(res_dt, chr, start, end)
  setkey(peak_dt, chr, start, end)

  res_annotated <- foverlaps(
    res_dt,
    peak_dt,
    type = "any",
    nomatch = NA
  )

  res_annotated <- res_annotated[, .(
    peak_id,
    chr,
    start,
    end,
    log2FoldChange,
    padj,
    SYMBOL,
    GENENAME,
    annotation,
    distanceToTSS
  )]

  res_annotated
}

res <- link_peaks_to_genes(res_day1_vs_ctrl)

gene_df_unique <- res %>%
  filter(!is.na(SYMBOL)) %>%
  arrange(desc(log2FoldChange)) %>%
  distinct(SYMBOL, .keep_all = TRUE)

gene_df_unique$ENSEMBL <- mapIds(
  org.Hs.eg.db,
  keys     = gene_df_unique$SYMBOL,
  keytype  = "SYMBOL",
  column   = "ENSEMBL",
  multiVals = "first"
)

gene_df_unique <- gene_df_unique %>%
  filter(!is.na(ENSEMBL))

gene_df_unique <- gene_df_unique %>%
  group_by(ENSEMBL) %>%
  slice_max(order_by = abs(log2FoldChange), n = 1, with_ties = FALSE) %>%
  ungroup()

symbol_map <- mapIds(
  org.Hs.eg.db,
  keys = gene_df_unique$ENSEMBL,
  keytype = "ENSEMBL",
  column = "SYMBOL",
  multiVals = "first"
)

TF_sets <- msigdbr(
  species = "Homo sapiens",
  category = "C3"
) %>%
  filter(
    grepl(
      "PU1|GATA1|CEBPA|FOXO",
      gs_name
    )
  ) %>%
  split(x = .$ensembl_gene, f = .$gs_name)

pu1_symbols <- symbol_map[TF_sets[["PU1_Q6"]]] %>%
  unname() %>%        # drop ENSG names
  na.omit() %>%       # remove ENSGs not in mapping
  unique()

gata1_symbols <- symbol_map[TF_sets[["GATA1_01"]]] %>%
  unname() %>%        # drop ENSG names
  na.omit() %>%       # remove ENSGs not in mapping
  unique()

objs <- list(res_day1_vs_ctrl, res_day14_vs_ctrl, res_day14_vs_day1)
names(objs) <- c("1_ctrl", "14_ctrl", "14_1")

dt <- lapply(seq_len(length(objs)), function(x) {
	print(x)

	dt <- objs[[x]]

	res <- link_peaks_to_genes(dt)

	df_plot <- as.data.frame(res)

	# define significance
	df_plot <- df_plot %>%
	  mutate(sig = padj < 0.05)

	df_plot$condition <- names(objs)[x]

	df_plot

}) %>% rbindlist()

dt <- as.data.table(dt)

pu1 <- dt %>% 
	filter(SYMBOL %in% pu1_symbols) %>%
  arrange(desc(log2FoldChange)) %>%
  distinct(SYMBOL, .keep_all = TRUE)

heat_df <- pu1 %>%
  group_by(SYMBOL, condition) %>%
  summarise(
    log2FC = log2FoldChange[which.max(abs(log2FoldChange))],
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(
    names_from = condition,
    values_from = log2FC
  )

mat <- heat_df %>%
  tibble::column_to_rownames("SYMBOL") %>%
  as.matrix()

# Replace NA with 0 (important for clustering)
mat[is.na(mat)] <- 0

pdf("v4/pu1_heat.pdf")

pheatmap(
  mat,
  # scale = "row",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  color = colorRampPalette(c("blue", "white", "red"))(100),
  border_color = NA,
  fontsize_row = 2   # try 4–8 depending on how many genes
)

dev.off()

gata1 <- dt %>% 
	filter(SYMBOL %in% gata1_symbols) %>%
  arrange(desc(log2FoldChange)) %>%
  distinct(SYMBOL, .keep_all = TRUE)

heat_df <- gata1 %>%
  group_by(SYMBOL, condition) %>%
  summarise(
    log2FC = log2FoldChange[which.max(abs(log2FoldChange))],
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(
    names_from = condition,
    values_from = log2FC
  )

mat <- heat_df %>%
  tibble::column_to_rownames("SYMBOL") %>%
  as.matrix()

# Replace NA with 0 (important for clustering)
mat[is.na(mat)] <- 0

pdf("v4/gata1_heat.pdf")

pheatmap(
  mat,
  # scale = "row",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  color = colorRampPalette(c("blue", "white", "red"))(100),
  border_color = NA,
  fontsize_row = 2   # try 4–8 depending on how many genes
)

dev.off()

###############################
link_peaks_to_genes <- function(dt) {
  res_dt <- as.data.table(as.data.frame(dt), keep.rownames = "peak_id")
  
  coords <- tstrsplit(res_dt$peak_id, ":", fixed = TRUE)

  res_dt[, `:=`(
    chr   = coords[[1]],
    start = as.integer(coords[[2]]),
    end   = as.integer(coords[[3]])
  )]

  peaks <- readPeakFile("../03_count_mat/merged_peaks_500bp.filtered.sorted.bed")

  peak_anno <- annotatePeak(
    peaks,
    TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
    annoDb = "org.Hs.eg.db",
    tssRegion = c(-2000, 2000)
  )

  anno_df <- as.data.frame(peak_anno)

  anno_df$peak_id <- paste(
    anno_df$seqnames,
    anno_df$start,
    anno_df$end,
    sep = ":"
  )

  peak_dt <- as.data.table(anno_df)

  setnames(peak_dt,
           old = c("seqnames"),
           new = c("chr"))

  # Ensure integer coordinates
  peak_dt[, `:=`(
    start = as.integer(start),
    end   = as.integer(end)
  )]

  setkey(res_dt, chr, start, end)
  setkey(peak_dt, chr, start, end)

  res_annotated <- foverlaps(
    res_dt,
    peak_dt,
    type = "any",
    nomatch = NA
  )

  res_annotated <- res_annotated[, .(
    peak_id,
    chr,
    start,
    end,
    log2FoldChange,
    padj,
    SYMBOL,
    GENENAME,
    annotation,
    distanceToTSS
  )]

  res_annotated
}

res <- link_peaks_to_genes(res_day1_vs_ctrl)

gene_df_unique <- res %>%
  filter(!is.na(SYMBOL)) %>%
  arrange(desc(log2FoldChange)) %>%
  distinct(SYMBOL, .keep_all = TRUE)

gene_df_unique$ENSEMBL <- mapIds(
  org.Hs.eg.db,
  keys     = gene_df_unique$SYMBOL,
  keytype  = "SYMBOL",
  column   = "ENSEMBL",
  multiVals = "first"
)

gene_df_unique <- gene_df_unique %>%
  filter(!is.na(ENSEMBL))

gene_df_unique <- gene_df_unique %>%
  group_by(ENSEMBL) %>%
  slice_max(order_by = abs(log2FoldChange), n = 1, with_ties = FALSE) %>%
  ungroup()

symbol_map <- mapIds(
  org.Hs.eg.db,
  keys = gene_df_unique$ENSEMBL,
  keytype = "ENSEMBL",
  column = "SYMBOL",
  multiVals = "first"
)

symbol_to_ensg <- setNames(names(symbol_map), symbol_map)

targets <- fread("TF_gene_pairs.tsv", header = FALSE)
targets <- targets[nchar(V2) > 0]
targets$Ensembl <- symbol_to_ensg[targets$V2]
colnames(targets) <- c("TF", "gene_name", "gene_id")
# targets <- targets[!is.na(gene_id)]

pu1_symbols <- targets[TF == "PU.1"]
gata1_symbols <- targets[TF == "GATA1"]
runx1_symbols <- targets[TF == "RUNX1"]
cebpa_symbols <- targets[TF == "CEBPA"]

objs <- list(res_day1_vs_ctrl, res_day14_vs_ctrl, res_day14_vs_day1)
names(objs) <- c("1_ctrl", "14_ctrl", "14_1")

dt <- lapply(seq_len(length(objs)), function(x) {
  print(x)

  dt <- objs[[x]]

  res <- link_peaks_to_genes(dt)

  df_plot <- as.data.frame(res)

  # define significance
  df_plot <- df_plot %>%
    mutate(sig = padj < 0.05)

  df_plot$condition <- names(objs)[x]

  df_plot

}) %>% rbindlist()

dt <- as.data.table(dt)

target_list <- list(pu1_symbols, gata1_symbols, runx1_symbols, cebpa_symbols)
names(target_list) <- c("PU1", "GATA1", "RUNX1", "CEBPA")

plots <- lapply(seq_len(length(target_list)), function(x) {
  tmp_targets <- target_list[[x]]$gene_name

  tmp <- dt %>% 
  filter(SYMBOL %in% tmp_targets) %>%
  arrange(desc(log2FoldChange)) %>%
  distinct(SYMBOL, .keep_all = TRUE)

  heat_df <- tmp %>%
    group_by(SYMBOL, condition) %>%
    summarise(
      log2FC = log2FoldChange[which.max(abs(log2FoldChange))],
      .groups = "drop"
    ) %>%
    tidyr::pivot_wider(
      names_from = condition,
      values_from = log2FC
    )

  mat <- heat_df %>%
    tibble::column_to_rownames("SYMBOL") %>%
    as.matrix()

  # Replace NA with 0 (important for clustering)
  mat[is.na(mat)] <- 0

  pheatmap(
    mat,
    main = paste0(names(target_list)[x], " Targets"),
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    clustering_method = "complete",
    color = colorRampPalette(c("blue", "white", "red"))(100),
    border_color = NA,
    show_rownames = FALSE,
    silent = TRUE
  )
})

pdf("v5/test.pdf")

for (p in plots) {
  grid::grid.newpage()     # force new page
  grid::grid.draw(p$gtable)
}

dev.off()





