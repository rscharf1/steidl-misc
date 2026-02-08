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

alpha <- 0.05

open_day1 <- rownames(res_day1_vs_ctrl)[
  res_day1_vs_ctrl$padj < alpha &
  res_day1_vs_ctrl$log2FoldChange > 0
]

close_day1 <- rownames(res_day1_vs_ctrl)[
  res_day1_vs_ctrl$padj < alpha &
  res_day1_vs_ctrl$log2FoldChange < 0
]

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

# open_day1_genes <- anno_df %>%
#   filter(peak_id %in% open_day1) %>%
#   filter(!is.na(SYMBOL))


# close_day1_genes <- anno_df %>%
#   filter(peak_id %in% close_day1) %>%
#   filter(!is.na(SYMBOL))

peaks_to_granges <- function(peaks) {
  parts <- do.call(rbind, strsplit(peaks, ":"))
  GRanges(
    seqnames = parts[,1],
    ranges = IRanges(
      start = as.integer(parts[,2]),
      end   = as.integer(parts[,3])
    )
  )
}

anno_gr <- GRanges(
  seqnames = anno_df$seqnames,
  ranges = IRanges(
    start = anno_df$start,
    end   = anno_df$end
  )
)

open_day1_gr  <- peaks_to_granges(open_day1)
close_day1_gr <- peaks_to_granges(close_day1)

open_hits  <- findOverlaps(anno_gr, open_day1_gr)
close_hits <- findOverlaps(anno_gr, close_day1_gr)

open_day1_genes <- anno_df[queryHits(open_hits), ] %>%
  filter(!is.na(SYMBOL))







