suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(data.table)
  library(apeglm)
})

counts_file <- "../04_counts/featureCounts_output/counts.txt"
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

count_matrix <- count_matrix[, rownames(coldata)]

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData   = coldata,
  design    = ~ condition
)

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

dds <- DESeq(dds)

res_1d <- results(
  dds,
  contrast = c("condition", "AraC_1_day", "control")
)

res_3d <- results(
  dds,
  contrast = c("condition", "AraC_3_day", "control")
)

res_1d_shrunk <- lfcShrink(
  dds,
  coef = "condition_AraC_1_day_vs_control",
  type = "apeglm"
)

res_3d_shrunk <- lfcShrink(
  dds,
  coef = "condition_AraC_3_day_vs_control",
  type = "apeglm"
)

pdf("01.pdf")

  plotMA(res_1d_shrunk, ylim = c(-4, 4),
         main = "AraC 1 day vs Control")

  plotMA(res_3d_shrunk, ylim = c(-4, 4),
         main = "AraC 3 day vs Control (exploratory)")

dev.off()

plot_volcano <- function(res, title) {
  res_df <- as.data.frame(res) %>%
    rownames_to_column("gene") %>%
    mutate(sig = padj < 0.05 & abs(log2FoldChange) > 1)

  ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color = sig), alpha = 0.6, size = 1) +
    scale_color_manual(values = c("grey70", "red")) +
    theme_classic() +
    labs(title = title, x = "log2 fold change", y = "-log10 adjusted p")
}

pdf("02.pdf")

plot_volcano(res_1d_shrunk, "AraC 1 day vs Control")
plot_volcano(res_3d_shrunk, "AraC 3 day vs Control (exploratory)")

vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup = "condition")

dev.off()

library(pheatmap)

top_genes <- res_1d_shrunk %>%
  as.data.frame() %>%
  filter(padj < 0.05) %>%
  arrange(padj) %>%
  head(50) %>%
  rownames()

mat <- assay(vsd)[top_genes, ]

pdf("03.pdf")

pheatmap(
  mat,
  scale = "row",
  annotation_col = coldata["condition"],
  show_rownames = FALSE,
  main = "Top DE genes (AraC 1 day)"
)

dev.off()

#################################
res_1d <- res_1d_shrunk

res_1d_df <- as.data.frame(res_1d) %>%
  rownames_to_column("gene") %>%
  filter(!is.na(padj))

res_1d_df$ensembl <- gsub("\\..*$", "", res_1d_df$gene)

library(org.Hs.eg.db)
library(AnnotationDbi)

symbol_map <- mapIds(
  org.Hs.eg.db,
  keys = res_1d_df$ensembl,
  keytype = "ENSEMBL",
  column = "SYMBOL",
  multiVals = "first"
)

res_1d_df$symbol <- symbol_map

res_1d_df_sym <- res_1d_df %>%
  filter(!is.na(symbol)) %>%
  distinct(symbol, .keep_all = TRUE)

gene_ranks_sym <- res_1d_df_sym$log2FoldChange
names(gene_ranks_sym) <- res_1d_df_sym$symbol

gene_ranks_sym <- sort(gene_ranks_sym, decreasing = TRUE)

library(fgsea)
library(msigdbr)

msig_h <- msigdbr(
  species = "Homo sapiens",
  category = "H"
)

hallmark_sets <- split(msig_h$gene_symbol, msig_h$gs_name)

fgsea_h <- fgsea(
  pathways = hallmark_sets,
  stats    = gene_ranks_sym
)

pdf("04.pdf")

top_h <- fgsea_h %>%
  filter(padj < 0.05) %>%
  arrange(desc(NES)) %>%
  head(15)

ggplot(top_h,
       aes(x = NES, y = reorder(pathway, NES))) +
  geom_point(aes(size = -log10(padj), color = NES)) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red") +
  theme_classic() +
  labs(
    title = "Hallmark pathways enriched after 1 day AraC treatment",
    x = "Normalized Enrichment Score",
    y = ""
  )

plotEnrichment(
  hallmark_sets[["HALLMARK_ALLOGRAFT_REJECTION"]],
  gene_ranks_sym
) +
  labs(title = "Inflammatory response pathway enrichment (AraC 1 day)")

dev.off()

























