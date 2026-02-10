library(dplyr)
library(data.table)
library(DESeq2)
library(ggplot2)

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

pdf("01_PCA.pdf")

	plotPCA(vsd, intgroup = "condition")

dev.off()

dds <- DESeq(dds)

res_day1_vs_ctrl <- results(dds, contrast = c("condition", "day_1", "control"))
res_day14_vs_ctrl <- results(dds, contrast = c("condition", "day_14", "control"))
res_day14_vs_day1 <- results(dds, contrast = c("condition", "day_14", "day_1"))

summary(res_day1_vs_ctrl)
summary(res_day14_vs_ctrl)
summary(res_day14_vs_day1)

df_day_1 <- as.data.frame(res_day1_vs_ctrl)
df_day_14 <- as.data.frame(res_day14_vs_ctrl)

pdf("02_Changes.pdf")

	ggplot() +
	  geom_density(data = df_day_1, aes(log2FoldChange), color = "blue") +
	  geom_density(data = df_day_14, aes(log2FoldChange), color = "green") +
	  geom_vline(xintercept = 0, linetype = "dashed") +
	  labs(
	    x = "log2 fold-change (accessibility)",
	    y = "Density",
	    title = "Global accessibility shifts vs control"
	  )

	df_compare <- data.frame(
	  day_1 = res_day1_vs_ctrl$log2FoldChange,
	  day_14 = res_day14_vs_ctrl$log2FoldChange
	)

	ggplot(df_compare, aes(x = day_1, y = day_14)) +
	  geom_point(size = 0.5, alpha = 0.5) +
	  geom_abline(intercept = 0, slope = 1, color = "red") +
	  labs(
	    x = "Day 1 vs control log2FC",
	    y = "Day 14 vs control log2FC"
	  ) +
	  theme_classic()

dev.off()


sum(res_day1_vs_ctrl$padj < 0.05, na.rm = TRUE)
sum(res_day14_vs_ctrl$padj < 0.05, na.rm = TRUE)

sum(res_day1_vs_ctrl$padj < 0.05 & res_day1_vs_ctrl$log2FoldChange > 0, na.rm = TRUE)
sum(res_day1_vs_ctrl$padj < 0.05 & res_day1_vs_ctrl$log2FoldChange < 0, na.rm = TRUE)

sum(res_day14_vs_ctrl$padj < 0.05 & res_day14_vs_ctrl$log2FoldChange > 0, na.rm = TRUE)
sum(res_day14_vs_ctrl$padj < 0.05 & res_day14_vs_ctrl$log2FoldChange < 0, na.rm = TRUE)

pdf("03_Day1_v_Ctrl.pdf")

	df_day1 <- as.data.frame(res_day1_vs_ctrl)

	ggplot(df_day1, aes(x = log2FoldChange)) +
	  geom_density(fill = "darkgreen", alpha = 0.4) +
	  geom_vline(xintercept = 0, linetype = "dashed") +
	  labs(
	    x = "log2 fold-change (Day 1 vs control)",
	    y = "Density",
	    title = "Global chromatin accessibility shifts in Day 1"
	  ) +
	  theme_classic()

	df_day1 %>%
	  filter(!is.na(log2FoldChange)) %>%
	  arrange(log2FoldChange) %>%
	  mutate(rank = row_number() / n()) %>%
	  ggplot(aes(x = log2FoldChange, y = rank)) +
	  geom_line(color = "darkgreen") +
	  geom_vline(xintercept = 0, linetype = "dashed") +
	  labs(
	    x = "log2 fold-change",
	    y = "Cumulative fraction of peaks",
	    title = "Cumulative accessibility changes (Day 1 vs control)"
	  ) +
	  theme_classic()

	ggplot(df_day1, aes(x = log2FoldChange, y = -log10(padj))) +
	  geom_point(size = 0.4, alpha = 0.4) +
	  geom_vline(xintercept = 0, linetype = "dashed") +
	  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
	  labs(
	    x = "log2 fold-change",
	    y = "-log10(FDR)",
	    title = "Differential accessibility: Day 1 vs control"
	  ) +
	  theme_classic()

	data.frame(
	  direction = c("Opening", "Closing"),
	  n = c(
	    sum(df_day1$padj < 0.05 & df_day1$log2FoldChange > 0, na.rm = TRUE),
	    sum(df_day1$padj < 0.05 & df_day1$log2FoldChange < 0, na.rm = TRUE)
	  )
	) %>%
	  ggplot(aes(x = direction, y = n, fill = direction)) +
	  geom_col() +
	  theme_classic() +
	  labs(
	    y = "Number of peaks",
	    title = "Direction of accessibility changes in Day 1"
	  )

dev.off()

pdf("04_Day14_v_Ctrl.pdf")

	df_day14 <- as.data.frame(res_day14_vs_ctrl)

	ggplot(df_day14, aes(x = log2FoldChange)) +
	  geom_density(fill = "darkgreen", alpha = 0.4) +
	  geom_vline(xintercept = 0, linetype = "dashed") +
	  labs(
	    x = "log2 fold-change (Day 14 vs control)",
	    y = "Density",
	    title = "Global chromatin accessibility shifts in Day 14"
	  ) +
	  theme_classic()

	df_day14 %>%
	  filter(!is.na(log2FoldChange)) %>%
	  arrange(log2FoldChange) %>%
	  mutate(rank = row_number() / n()) %>%
	  ggplot(aes(x = log2FoldChange, y = rank)) +
	  geom_line(color = "darkgreen") +
	  geom_vline(xintercept = 0, linetype = "dashed") +
	  labs(
	    x = "log2 fold-change",
	    y = "Cumulative fraction of peaks",
	    title = "Cumulative accessibility changes (Day 14 vs control)"
	  ) +
	  theme_classic()

	ggplot(df_day14, aes(x = log2FoldChange, y = -log10(padj))) +
	  geom_point(size = 0.4, alpha = 0.4) +
	  geom_vline(xintercept = 0, linetype = "dashed") +
	  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
	  labs(
	    x = "log2 fold-change",
	    y = "-log10(FDR)",
	    title = "Differential accessibility: Day 14 vs control"
	  ) +
	  theme_classic()

	data.frame(
	  direction = c("Opening", "Closing"),
	  n = c(
	    sum(df_day14$padj < 0.05 & df_day14$log2FoldChange > 0, na.rm = TRUE),
	    sum(df_day14$padj < 0.05 & df_day14$log2FoldChange < 0, na.rm = TRUE)
	  )
	) %>%
	  ggplot(aes(x = direction, y = n, fill = direction)) +
	  geom_col() +
	  theme_classic() +
	  labs(
	    y = "Number of peaks",
	    title = "Direction of accessibility changes in Day 14"
	  )

dev.off()

# Venn

alpha <- 0.05

open_day1 <- rownames(res_day1_vs_ctrl)[
  res_day1_vs_ctrl$padj < alpha &
  res_day1_vs_ctrl$log2FoldChange > 0
]

open_day14 <- rownames(res_day14_vs_ctrl)[
  res_day14_vs_ctrl$padj < alpha &
  res_day14_vs_ctrl$log2FoldChange > 0
]

close_day1 <- rownames(res_day1_vs_ctrl)[
  res_day1_vs_ctrl$padj < alpha &
  res_day1_vs_ctrl$log2FoldChange < 0
]

close_day14 <- rownames(res_day14_vs_ctrl)[
  res_day14_vs_ctrl$padj < alpha &
  res_day14_vs_ctrl$log2FoldChange < 0
]

library(VennDiagram)
library(grid)

pdf("05_venn.pdf", width = 8, height = 8)

	venn.plot <- venn.diagram(
	  x = list(
	    "Day 1 opening"  = open_day1,
	    "Day 14 opening" = open_day14
	  ),
	  filename = NULL,
	  fill = c("lightblue", "darkgreen"),
	  alpha = 0.5,
	  cex = 1.5,
	  cat.cex = 1.2,
	  cat.pos = c(0, 0),     # angles (left, right)
	  cat.dist = c(-0.05, -0.05),
	  main = "Opening chromatin regions"
	)

	grid.draw(venn.plot)
	grid.newpage()

	venn.plot <- venn.diagram(
	  x = list(
	    "Day 1 closing"  = close_day1,
	    "Day 14 closing" = close_day14
	  ),
	  filename = NULL,          # <-- critical
	  fill = c("lightblue", "darkgreen"),
	  alpha = 0.5,
	  cex = 1.5,
	  cat.cex = 1.2,
  	  cat.pos = c(0, 0),     # angles (left, right)
	  cat.dist = c(-0.05, -0.05),
	  main = "Closing chromatin regions"
	)

	grid.draw(venn.plot)

dev.off()

library(UpSetR)

all_peaks <- rownames(res_day1_vs_ctrl)

upset_mat <- data.frame(
  open_day1  = as.integer(all_peaks %in% open_day1),
  open_day14  = as.integer(all_peaks %in% open_day14),
  close_day1 = as.integer(all_peaks %in% close_day1),
  close_day14 = as.integer(all_peaks %in% close_day14)
)

# Force base data.frame (not tibble)
upset_mat <- as.data.frame(upset_mat)

pdf("06_upset_peaks.pdf", width = 8, height = 5)

upset(
  upset_mat,
  nsets = 4,
  nintersects = 15,
  order.by = "freq",
  main.bar.color = "black",
  sets.bar.color = "grey40",
  text.scale = 1.4
)

dev.off()

####
# Gene Level Analysis

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(GenomicRanges)
library(ggrepel)

res_day1_vs_ctrl
res_day14_vs_ctrl
res_day14_vs_day1

dt <- res_day14_vs_day1
comparison <- "Day 14 v. Day 1"
file_names <- "14_1"

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

res <- link_peaks_to_genes(dt)

df_plot <- as.data.frame(res)

# define significance
df_plot <- df_plot %>%
  mutate(sig = padj < 0.05)

pdf(paste0("v2/volcano_", file_names, ".pdf"))

	label_candidates <- df_plot %>%
	filter(
	  sig,
	  !is.na(SYMBOL),
	  !grepl("^LOC", SYMBOL),
	  !grepl("uncharacterized", GENENAME, ignore.case = TRUE),
	  !grepl("pseudogene", GENENAME, ignore.case = TRUE)
	)

	label_up <- label_candidates %>%
	  filter(log2FoldChange > 0) %>%
	  arrange(desc(log2FoldChange)) %>%
	  slice_head(n = 10)

	label_down <- label_candidates %>%
	  filter(log2FoldChange < 0) %>%
	  arrange(log2FoldChange) %>%
	  slice_head(n = 10)

	label_df <- bind_rows(label_up, label_down)

	ggplot(df_plot, aes(x = log2FoldChange, y = -log10(padj))) +
	  geom_point(size = 0.4, alpha = 0.4, color = "grey50") +
	  geom_point(
	    data = df_plot %>% filter(sig),
	    size = 0.5,
	    alpha = 0.6,
	    color = "firebrick"
	  ) +
	  geom_text_repel(
	    data = label_df,
	    aes(label = SYMBOL),
	    size = 3,
	    max.overlaps = Inf,
	    box.padding = 0.4,
	    point.padding = 0.3
	  ) +
	  geom_vline(xintercept = 0, linetype = "dashed") +
	  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
	  labs(
	    x = paste0("log2 fold-change ", comparison),
	    y = "-log10(FDR)",
	    title = paste0("Differential chromatin accessibility: ", comparison)
	  ) +
	  theme_classic()

	label_df <- df_plot %>%
		filter(
			grepl("RUNX1$|CEBPA$|PU1|GATA1|MYC$|FOXO4|SPI1", SYMBOL, ignore.case=TRUE)
		) %>%
	  arrange(desc(log2FoldChange)) %>%
	  distinct(SYMBOL, .keep_all = TRUE)

	ggplot(df_plot, aes(x = log2FoldChange, y = -log10(padj))) +
	  geom_point(size = 0.4, alpha = 0.4, color = "grey50") +
	  geom_point(
	    data = df_plot %>% filter(sig),
	    size = 0.5,
	    alpha = 0.6,
	    color = "firebrick"
	  ) +
	  geom_text_repel(
	    data = label_df,
	    aes(label = SYMBOL),
	    size = 3,
	    max.overlaps = Inf,
	    box.padding = 0.4,
	    point.padding = 0.3
	  ) +
	  geom_vline(xintercept = 0, linetype = "dashed") +
	  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
	  labs(
	    x = paste0("log2 fold-change ", comparison),
	    y = "-log10(FDR)",
	    title = paste0("Differential chromatin accessibility: ", comparison)
	  ) +
	  theme_classic()

dev.off()

# GSEA Analysis 
library(org.Hs.eg.db)
library(AnnotationDbi)
library(msigdbr)

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

ranks <- gene_df_unique %>%
  arrange(desc(log2FoldChange)) %>%
  dplyr::select(ENSEMBL, log2FoldChange)

ranked_genes <- ranks$log2FoldChange
names(ranked_genes) <- ranks$ENSEMBL

hallmark_sets <- msigdbr(
  species  = "Homo sapiens",
  category = "H"
) %>%
  split(x = .$ensembl_gene, f = .$gs_name)

  # TF target genes 
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

symbol_map <- mapIds(
  org.Hs.eg.db,
  keys = gene_df_unique$ENSEMBL,
  keytype = "ENSEMBL",
  column = "SYMBOL",
  multiVals = "first"
)

runx_targets <- fread("RUNX1_targets_JASPAR.csv")

symbol_df <- tibble::enframe(symbol_map, name = "ENSG", value = "Symbol") %>%
  filter(!is.na(Symbol))

runx_targets_ensg <- runx_targets %>%
  left_join(symbol_df, by = "Symbol")

library(fgsea)

fgsea_res <- fgsea(
  pathways = hallmark_sets,
  stats    = ranked_genes
)

fgsea_res[order(fgsea_res$padj), ][, c("pathway", "NES", "padj")][1:10, ]

top_fgsea <- fgsea_res %>%
  filter(padj < 0.1)

pdf(paste0("v2/GSEA_", file_names, ".pdf"))

  ggplot(top_fgsea,
         aes(x = reorder(pathway, NES), y = NES, fill = NES > 0)) +
    geom_col() +
    coord_flip() +
    scale_fill_manual(values = c("steelblue", "firebrick")) +
    labs(
      x = "Hallmark pathway",
      y = "Normalized Enrichment Score (NES)",
      title = paste0("GSEA of ", comparison)
    ) +
    theme_classic()

  plotEnrichment(
    hallmark_sets[["HALLMARK_TGF_BETA_SIGNALING"]],
    ranked_genes
  ) +
    labs(title = "GSEA: TGF Beta")

  plotEnrichment(
    hallmark_sets[["HALLMARK_IL2_STAT5_SIGNALING"]],
    ranked_genes
  ) +
    labs(title = "GSEA: IL2 STAT5")

  plotEnrichment(
    hallmark_sets[["HALLMARK_MITOTIC_SPINDLE"]],
    ranked_genes
  ) +
    labs(title = "GSEA: Mitotic Spindle")

  plotEnrichment(
    TF_sets[["PU1_Q6"]],
    ranked_genes
  ) +
    labs(title = "GSEA: PU1 Targets")

  plotEnrichment(
    TF_sets[["GATA1_01"]],
    ranked_genes
  ) +
    labs(title = "GSEA: GATA1 Targets")

  plotEnrichment(
    TF_sets[["CEBPA_01"]],
    ranked_genes
  ) +
    labs(title = "GSEA: CEBPA Targets")

  plotEnrichment(
    TF_sets[["FOXO4_TARGET_GENES"]],
    ranked_genes
  ) +
    labs(title = "GSEA: FOXO4 Targets")

  plotEnrichment(
    as.vector(na.omit(runx_targets_ensg$ENSG)),
    ranked_genes
  ) +
    labs(title = "JASPAR: RUNX1 Targets")

dev.off()

# Label PU1 and GATA1 targets in volcano plots 
# PU1 targets 
pu1_symbols <- symbol_map[TF_sets[["PU1_Q6"]]] %>%
  unname() %>%        # drop ENSG names
  na.omit() %>%       # remove ENSGs not in mapping
  unique()

gata1_symbols <- symbol_map[TF_sets[["GATA1_01"]]] %>%
  unname() %>%        # drop ENSG names
  na.omit() %>%       # remove ENSGs not in mapping
  unique()

df_plot %>% filter(!is.na(SYMBOL)) %>% head()

pdf(paste0("v2/volcano_TF_", file_names, ".pdf"))

	label_df <- df_plot %>% 
		filter(SYMBOL %in% pu1_symbols) %>%
	  arrange(desc(log2FoldChange)) %>%
	  distinct(SYMBOL, .keep_all = TRUE)

	ggplot(df_plot, aes(x = log2FoldChange, y = -log10(padj))) +
	  geom_point(size = 0.4, alpha = 0.4, color = "grey50") +
	  geom_point(
	    data = label_df,
	    size = 0.8,
	    alpha = 0.9,
	    color = "dodgerblue3"
	  ) +
	  geom_vline(xintercept = 0, linetype = "dashed") +
	  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
	  labs(
	    x = paste0("log2 fold-change ", comparison),
	    y = "-log10(FDR)",
	    title = paste0("PU1 Targets: ", comparison)
	  ) +
	  theme_classic()

	label_df <- df_plot %>% 
		filter(SYMBOL %in% gata1_symbols) %>%
	  arrange(desc(log2FoldChange)) %>%
	  distinct(SYMBOL, .keep_all = TRUE)

	ggplot(df_plot, aes(x = log2FoldChange, y = -log10(padj))) +
	  geom_point(size = 0.4, alpha = 0.4, color = "grey50") +
	  geom_point(
	    data = label_df,
	    size = 0.8,
	    alpha = 0.9,
	    color = "dodgerblue3"
	  ) +
	  geom_vline(xintercept = 0, linetype = "dashed") +
	  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
	  labs(
	    x = paste0("log2 fold-change ", comparison),
	    y = "-log10(FDR)",
	    title = paste0("GATA1 Targets: ", comparison)
	  ) +
	  theme_classic()

dev.off()





