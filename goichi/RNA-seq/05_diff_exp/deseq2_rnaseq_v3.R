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

gene_lengths <- counts %>%
  dplyr::select(Geneid, Length) %>%
  distinct() %>%
  column_to_rownames("Geneid")

counts_to_tpm <- function(counts, lengths) {
  lengths <- lengths[rownames(counts), , drop = FALSE]
  rpk <- counts / (lengths$Length / 1000)
  tpm <- sweep(rpk, 2, colSums(rpk), "/") * 1e6
  return(tpm)
}

tpm <- counts_to_tpm(count_matrix, gene_lengths)
log_tpm <- log2(tpm + 1)

pca_df <- prcomp(t(log_tpm), scale. = FALSE)$x %>%
  as.data.frame() %>%
  rownames_to_column("sample")

pca_df <- pca_df %>%
  left_join(
    coldata %>% rownames_to_column("sample"),
    by = "sample"
  )

pdf("v3/PCA.pdf")

  ggplot(pca_df, aes(PC1, PC2, color = condition, shape = protocol)) +
    geom_point(size = 3) +
    theme_classic()

dev.off()

# BiocManager::install("GSVA")
library(GSVA)

hallmark_sets <- msigdbr(
  species  = "Homo sapiens",
  category = "H"
) %>%
  split(x = .$ensembl_gene, f = .$gs_name)

# strip version numbers
rownames(log_tpm) <- sub("\\..*$", "", rownames(log_tpm))

ssgsea_param <- ssgseaParam(
  exprData = as.matrix(log_tpm),
  geneSets = hallmark_sets
)

ssgsea_scores <- gsva(param = ssgsea_param)

# Trajectory plots 
ssgsea_df <- t(ssgsea_scores) %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  left_join(
    coldata %>% rownames_to_column("sample"),
    by = "sample"
  )

pathways_to_plot <- c(
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_COAGULATION",
  "HALLMARK_APOPTOSIS",
  "HALLMARK_COMPLEMENT",
  "HALLMARK_P53_PATHWAY",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_ALLOGRAFT_REJECTION",
  "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_MYC_TARGETS_V2",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_FATTY_ACID_METABOLISM",
  "HALLMARK_GLYCOLYSIS",
  "HALLMARK_CHOLESTEROL_HOMEOSTASIS",
  "HALLMARK_BILE_ACID_METABOLISM"
)

ssgsea_long <- ssgsea_df %>%
  dplyr::select(sample, condition, protocol, all_of(pathways_to_plot)) %>%
  pivot_longer(
    cols = all_of(pathways_to_plot),
    names_to = "pathway",
    values_to = "score"
  )

# enforce ordering
ssgsea_long$condition <- factor(
  ssgsea_long$condition,
  levels = c("control", "AraC_1_day", "AraC_3_day")
)

plot_trajectory <- function(df, pathway_name) {
  ggplot(
    df %>% filter(pathway == pathway_name),
    aes(x = condition, y = score,
        color = protocol, group = interaction(protocol, sample))
  ) +
    geom_point(size = 3, alpha = 0.9) +
    geom_line(alpha = 0.4) +
    stat_summary(
      aes(group = protocol),
      fun = mean,
      geom = "line",
      linewidth = 1.2
    ) +
    stat_summary(
      aes(group = protocol),
      fun = mean,
      geom = "point",
      size = 4
    ) +
    theme_classic(base_size = 12) +
    labs(
      title = pathway_name,
      x = NULL,
      y = "ssGSEA score"
    )
}

pdf("v3/All_Trajectory.pdf", width = 15, height = 15)

  ggplot(
    ssgsea_long,
    aes(x = condition, y = score,
        color = protocol, group = interaction(protocol, sample))
  ) +
    geom_point(size = 2.5, alpha = 0.9) +
    geom_line(alpha = 0.3) +
    stat_summary(
      aes(group = protocol),
      fun = mean,
      geom = "line",
      linewidth = 1.2
    ) +
    stat_summary(
      aes(group = protocol),
      fun = mean,
      geom = "point",
      size = 3.5
    ) +
    facet_wrap(~ pathway, scales = "free_y", ncol = 3) +
    theme_classic(base_size = 12) +
    labs(
      x = NULL,
      y = "ssGSEA score",
      color = "Protocol",
      title = "Trajectory of AraC-induced transcriptional programs"
    )

dev.off()

pdf("v3/Selected_Trajectory.pdf")

  p1 <- plot_trajectory(ssgsea_long, "HALLMARK_INFLAMMATORY_RESPONSE")
  p2 <- plot_trajectory(ssgsea_long, "HALLMARK_MYC_TARGETS_V1")
  p3 <- plot_trajectory(ssgsea_long, "HALLMARK_OXIDATIVE_PHOSPHORYLATION")

  print(p1)
  print(p2)
  print(p3)

dev.off()














