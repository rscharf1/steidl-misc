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
    grepl("07_|08_", sample) ~ "py_lo",
    grepl("09_|10_", sample) ~ "py_hi",
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

res_pylo_vs_ctrl <- results(dds, contrast = c("condition", "py_lo", "control"))
res_pyhi_vs_ctrl <- results(dds, contrast = c("condition", "py_hi", "control"))
res_pyhi_vs_pylo <- results(dds, contrast = c("condition", "py_hi", "py_lo"))

summary(res_pylo_vs_ctrl)
summary(res_pyhi_vs_ctrl)
summary(res_pyhi_vs_pylo)

df_lo <- as.data.frame(res_pylo_vs_ctrl)
df_hi <- as.data.frame(res_pyhi_vs_ctrl)

pdf("02_Changes.pdf")

	ggplot() +
	  geom_density(data = df_lo, aes(log2FoldChange), color = "blue") +
	  geom_density(data = df_hi, aes(log2FoldChange), color = "green") +
	  geom_vline(xintercept = 0, linetype = "dashed") +
	  labs(
	    x = "log2 fold-change (accessibility)",
	    y = "Density",
	    title = "Global accessibility shifts vs control"
	  )

	df_compare <- data.frame(
	  pylo = res_pylo_vs_ctrl$log2FoldChange,
	  pyhi = res_pyhi_vs_ctrl$log2FoldChange
	)

	ggplot(df_compare, aes(x = pylo, y = pyhi)) +
	  geom_point(size = 0.5, alpha = 0.5) +
	  geom_abline(intercept = 0, slope = 1, color = "red") +
	  labs(
	    x = "py-lo vs control log2FC",
	    y = "py-hi vs control log2FC"
	  ) +
	  theme_classic()

dev.off()


sum(res_pylo_vs_ctrl$padj < 0.05, na.rm = TRUE)
sum(res_pyhi_vs_ctrl$padj < 0.05, na.rm = TRUE)

sum(res_pyhi_vs_ctrl$padj < 0.05 & res_pyhi_vs_ctrl$log2FoldChange > 0, na.rm = TRUE)
sum(res_pyhi_vs_ctrl$padj < 0.05 & res_pyhi_vs_ctrl$log2FoldChange < 0, na.rm = TRUE)

pdf("03_High_v_Ctrl.pdf")

	df_hi <- as.data.frame(res_pyhi_vs_ctrl)

	ggplot(df_hi, aes(x = log2FoldChange)) +
	  geom_density(fill = "darkgreen", alpha = 0.4) +
	  geom_vline(xintercept = 0, linetype = "dashed") +
	  labs(
	    x = "log2 fold-change (py-hi vs control)",
	    y = "Density",
	    title = "Global chromatin accessibility shifts in py-hi"
	  ) +
	  theme_classic()

	df_hi %>%
	  filter(!is.na(log2FoldChange)) %>%
	  arrange(log2FoldChange) %>%
	  mutate(rank = row_number() / n()) %>%
	  ggplot(aes(x = log2FoldChange, y = rank)) +
	  geom_line(color = "darkgreen") +
	  geom_vline(xintercept = 0, linetype = "dashed") +
	  labs(
	    x = "log2 fold-change",
	    y = "Cumulative fraction of peaks",
	    title = "Cumulative accessibility changes (py-hi vs control)"
	  ) +
	  theme_classic()

	ggplot(df_hi, aes(x = log2FoldChange, y = -log10(padj))) +
	  geom_point(size = 0.4, alpha = 0.4) +
	  geom_vline(xintercept = 0, linetype = "dashed") +
	  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
	  labs(
	    x = "log2 fold-change",
	    y = "-log10(FDR)",
	    title = "Differential accessibility: py-hi vs control"
	  ) +
	  theme_classic()

	data.frame(
	  direction = c("Opening", "Closing"),
	  n = c(
	    sum(df_hi$padj < 0.05 & df_hi$log2FoldChange > 0, na.rm = TRUE),
	    sum(df_hi$padj < 0.05 & df_hi$log2FoldChange < 0, na.rm = TRUE)
	  )
	) %>%
	  ggplot(aes(x = direction, y = n, fill = direction)) +
	  geom_col() +
	  theme_classic() +
	  labs(
	    y = "Number of peaks",
	    title = "Direction of accessibility changes in py-hi"
	  )

dev.off()


pdf("04_Low_v_Ctrl.pdf")

	df_low <- as.data.frame(res_pylo_vs_ctrl)

	ggplot(df_low, aes(x = log2FoldChange)) +
	  geom_density(fill = "darkgreen", alpha = 0.4) +
	  geom_vline(xintercept = 0, linetype = "dashed") +
	  labs(
	    x = "log2 fold-change (py-lo vs control)",
	    y = "Density",
	    title = "Global chromatin accessibility shifts in py-lo"
	  ) +
	  theme_classic()

	df_low %>%
	  filter(!is.na(log2FoldChange)) %>%
	  arrange(log2FoldChange) %>%
	  mutate(rank = row_number() / n()) %>%
	  ggplot(aes(x = log2FoldChange, y = rank)) +
	  geom_line(color = "darkgreen") +
	  geom_vline(xintercept = 0, linetype = "dashed") +
	  labs(
	    x = "log2 fold-change",
	    y = "Cumulative fraction of peaks",
	    title = "Cumulative accessibility changes (py-lo vs control)"
	  ) +
	  theme_classic()

	ggplot(df_low, aes(x = log2FoldChange, y = -log10(padj))) +
	  geom_point(size = 0.4, alpha = 0.4) +
	  geom_vline(xintercept = 0, linetype = "dashed") +
	  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
	  labs(
	    x = "log2 fold-change",
	    y = "-log10(FDR)",
	    title = "Differential accessibility: py-lo vs control"
	  ) +
	  theme_classic()

	data.frame(
	  direction = c("Opening", "Closing"),
	  n = c(
	    sum(df_low$padj < 0.05 & df_low$log2FoldChange > 0, na.rm = TRUE),
	    sum(df_low$padj < 0.05 & df_low$log2FoldChange < 0, na.rm = TRUE)
	  )
	) %>%
	  ggplot(aes(x = direction, y = n, fill = direction)) +
	  geom_col() +
	  theme_classic() +
	  labs(
	    y = "Number of peaks",
	    title = "Direction of accessibility changes in py-lo"
	  )

dev.off()

# Venn

alpha <- 0.05

open_lo <- rownames(res_pylo_vs_ctrl)[
  res_pylo_vs_ctrl$padj < alpha &
  res_pylo_vs_ctrl$log2FoldChange > 0
]

open_hi <- rownames(res_pyhi_vs_ctrl)[
  res_pyhi_vs_ctrl$padj < alpha &
  res_pyhi_vs_ctrl$log2FoldChange > 0
]

close_lo <- rownames(res_pylo_vs_ctrl)[
  res_pylo_vs_ctrl$padj < alpha &
  res_pylo_vs_ctrl$log2FoldChange < 0
]

close_hi <- rownames(res_pyhi_vs_ctrl)[
  res_pyhi_vs_ctrl$padj < alpha &
  res_pyhi_vs_ctrl$log2FoldChange < 0
]

library(VennDiagram)
library(grid)

pdf("05_venn.pdf", width = 8, height = 8)

	venn.plot <- venn.diagram(
	  x = list(
	    "py-low opening"  = open_lo,
	    "py-high opening" = open_hi
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
	    "py-low closing"  = close_lo,
	    "py-high closing" = close_hi
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

all_peaks <- rownames(res_pyhi_vs_ctrl)

upset_df <- data.frame(
  peak = all_peaks,
  open_lo  = all_peaks %in% open_lo,
  open_hi  = all_peaks %in% open_hi,
  close_lo = all_peaks %in% close_lo,
  close_hi = all_peaks %in% close_hi
)

pdf("06_upset_peaks.pdf", width = 8, height = 5)

	upset(
	  upset_df[, c("open_lo", "open_hi", "close_lo", "close_hi")],
	  sets = c("open_lo", "open_hi", "close_lo", "close_hi"),
	  order.by = "freq",
	  empty.intersections = "on",
	  main.bar.color = "black",
	  sets.bar.color = "grey40",
	  text.scale = 1.4
	)

dev.off()

all_peaks <- rownames(res_pyhi_vs_ctrl)

upset_mat <- data.frame(
  open_lo  = as.integer(all_peaks %in% open_lo),
  open_hi  = as.integer(all_peaks %in% open_hi),
  close_lo = as.integer(all_peaks %in% close_lo),
  close_hi = as.integer(all_peaks %in% close_hi)
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









