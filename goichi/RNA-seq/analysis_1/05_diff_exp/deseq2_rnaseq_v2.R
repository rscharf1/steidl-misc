suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(data.table)
  library(apeglm)
  library(ggrepel)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(pheatmap)
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
  dplyr::select(-Chr, -Start, -End, -Strand, -Length) %>%
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

keep_samples <- rownames(coldata)[
  coldata$condition %in% c("control", "AraC_1_day")
]

count_matrix_sub <- count_matrix[, keep_samples]
coldata_sub      <- coldata[keep_samples, ]

coldata_sub$condition <- factor(
  coldata_sub$condition,
  levels = c("control", "AraC_1_day")
)

coldata_sub$protocol <- factor(coldata_sub$protocol)

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix_sub,
  colData   = coldata_sub,
  design    = ~ condition
)

dds <- dds[rowSums(counts(dds)) >= 10, ]

dds <- estimateSizeFactors(dds)

############
# PCA Plot 
  # Do AraC-treated samples globally shift away from control?
  # Ideally the two treated smaples cluster together 

vsd <- vst(dds, blind = TRUE)

pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

pdf("PCA.pdf")

  plotPCA(vsd, intgroup = "condition")

dev.off()

#########
# MA Plot
dds <- DESeq(dds)

res <- results(
  dds,
  contrast = c("condition", "AraC_1_day", "control")
)

summary(res)

resLFC <- lfcShrink(
  dds,
  coef = "condition_AraC_1_day_vs_control",
  type = "apeglm"
)

resLFC <- as.data.frame(resLFC) %>%
  rownames_to_column("gene") %>%
  filter(!is.na(padj))

resLFC$ensembl <- gsub("\\..*$", "", resLFC$gene)

symbol_map <- mapIds(
  org.Hs.eg.db,
  keys = resLFC$ensembl,
  keytype = "ENSEMBL",
  column = "SYMBOL",
  multiVals = "first"
)

resLFC$symbol <- symbol_map

resLFC <- resLFC %>%
  filter(!is.na(symbol)) %>%
  distinct(symbol, .keep_all = TRUE)

ma_df <- data.frame(
  symbol = resLFC$symbol,
  baseMean = resLFC$baseMean,
  log2FC   = resLFC$log2FoldChange,
  padj     = resLFC$padj
)

ma_df$significant <- ifelse(
  !is.na(ma_df$padj) & ma_df$padj < 0.1,
  "FDR < 0.1",
  "NS"
)

pdf("MA.pdf")

  ggplot(ma_df, aes(x = log10(baseMean + 1), y = log2FC)) +
    geom_point(aes(color = significant), size = 0.6, alpha = 0.6) +
    scale_color_manual(values = c("grey70", "red")) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    ylim(-5, 5) +
    labs(
      x = "log10 mean normalized counts",
      y = "Shrunken log2 fold change (AraC 1-day vs control)",
      title = "MA plot of early AraC response"
    ) +
    theme_classic()

dev.off()

# Volcano plot
  # Label the TFs that Goichi identified

volcano_df <- as.data.frame(resLFC) %>%
  mutate(
    negLog10Padj = -log10(padj),
    significant = case_when(
      !is.na(padj) & padj < 0.1 & abs(log2FoldChange) >= 1 ~ "FDR < 0.1 & |LFC| ≥ 1",
      !is.na(padj) & padj < 0.1 ~ "FDR < 0.1",
      TRUE ~ "NS"
    )
  )

top_genes <- volcano_df %>%
  filter(significant == "FDR < 0.1 & |LFC| ≥ 1") %>%
  arrange(padj) %>%
  head(10)

notable <- c("SPI1", "GATA1", "RUNX1", "CEBPA")

goichi_tfs <- volcano_df %>% 
  filter(symbol %in% notable)

volcano_df$symbol
grep("GATA", volcano_df$symbol) %>% table

volcano_df$symbol[grep("GATA", volcano_df$symbol)]

pdf("v2/Volcano.pdf")

ggplot(volcano_df, aes(x = log2FoldChange, y = negLog10Padj)) +
  geom_point(aes(color = significant), alpha = 0.6, size = 0.7) +
  geom_text_repel(
    data = top_genes,
    aes(label = symbol),
    size = 3,
    max.overlaps = 10
  ) +
  scale_color_manual(
    values = c("grey70", "steelblue", "firebrick")
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.1), linetype = "dashed") +
  theme_classic()

ggplot(volcano_df, aes(x = log2FoldChange, y = negLog10Padj)) +
  geom_point(aes(color = significant), alpha = 0.6, size = 0.7) +
  geom_text_repel(
    data = goichi_tfs,
    aes(label = symbol),
    size = 3,
    max.overlaps = 10
  ) +
  scale_color_manual(
    values = c("grey70", "steelblue", "firebrick")
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.1), linetype = "dashed") +
  theme_classic()

dev.off()

volcano_df <- as.data.frame(resLFC) %>%
  mutate(
    negLog10Padj = -log10(padj),
    significant = case_when(
      !is.na(padj) & padj < 0.1 & abs(log2FoldChange) >= 0.5 ~ "FDR < 0.1 & |LFC| ≥ 0.5",
      !is.na(padj) & padj < 0.1 ~ "FDR < 0.1",
      TRUE ~ "NS"
    )
  )

top_genes <- volcano_df %>%
  filter(significant == "FDR < 0.1 & |LFC| ≥ 0.5") %>%
  arrange(padj) %>%
  head(10)

notable <- c("SPI1", "GATA1", "RUNX1", "CEBPA")

goichi_tfs <- volcano_df %>% 
  filter(symbol %in% notable)

volcano_df$symbol
grep("GATA", volcano_df$symbol) %>% table

volcano_df$symbol[grep("GATA", volcano_df$symbol)]

#########
# Heatmap of top DE genes 
gene_map <- resLFC %>%
  dplyr::select(gene, symbol) %>%
  dplyr::distinct()

top_genes <- resLFC %>%
  arrange(padj) %>%
  filter(!is.na(padj)) %>%
  slice_head(n = 40) %>%
  pull(gene)

# Extract VST expression
mat <- assay(vsd)[top_genes, ]

new_names <- gene_map$symbol[
  match(rownames(mat), gene_map$gene)
]

# Replace rownames, keeping Ensembl if symbol is missing
rownames(mat) <- ifelse(
  is.na(new_names) | new_names == "",
  rownames(mat),
  new_names
)

mat_scaled <- t(scale(t(mat)))

# Annotation
annotation_col <- data.frame(
  condition = coldata_sub$condition
)

colnames(mat_scaled) <- paste0(
  coldata_sub$condition, "_", seq_len(ncol(mat_scaled))
)

annotation_col <- data.frame(
  condition = coldata_sub$condition
)
rownames(annotation_col) <- colnames(mat_scaled)

pdf("heatmap.pdf")

  pheatmap(
    mat_scaled,
    annotation_col = annotation_col,
    cluster_cols = FALSE,
    show_colnames = FALSE,
    show_rownames = TRUE,
    fontsize_row = 7,
    main = "Top differentially expressed genes (AraC 1-day vs control)"
  )

dev.off()

# Heatmap of early TFs 
gene_map <- resLFC %>%
  dplyr::select(gene, symbol) %>%
  dplyr::distinct()

tfs <- c(
  "SPI1",
  "RUNX1",
  "CEBPA"
)

"PU.1" %in% resLFC$symbol

resLFC$symbol[grep("GATA", resLFC$symbol)]

top_genes <- resLFC %>%
  filter(symbol %in% tfs) %>%
  pull(gene)

# Extract VST expression
mat <- assay(vsd)[top_genes, ]

new_names <- gene_map$symbol[
  match(rownames(mat), gene_map$gene)
]

# Replace rownames, keeping Ensembl if symbol is missing
rownames(mat) <- ifelse(
  is.na(new_names) | new_names == "",
  rownames(mat),
  new_names
)

mat_scaled <- t(scale(t(mat)))

# Annotation
annotation_col <- data.frame(
  condition = coldata_sub$condition
)

colnames(mat_scaled) <- paste0(
  coldata_sub$condition, "_", seq_len(ncol(mat_scaled))
)

annotation_col <- data.frame(
  condition = coldata_sub$condition
)
rownames(annotation_col) <- colnames(mat_scaled)

pdf("v2/heatmap_TFs.pdf")

  pheatmap(
    mat_scaled,
    annotation_col = annotation_col,
    cluster_cols = FALSE,
    show_colnames = FALSE,
    show_rownames = TRUE,
    fontsize_row = 7,
    main = "Notable TFs"
  )

dev.off()

########
# More general heatmap 
mat <- assay(vsd)

gene_var <- apply(mat, 1, var)

top_n <- 2000  # 500–2000 is reasonable
top_genes <- names(sort(gene_var, decreasing = TRUE))[1:top_n]

mat_top <- mat[top_genes, ]

mat_z <- t(scale(t(mat_top)))

annotation_col <- as.data.frame(coldata_sub["condition"])
rownames(annotation_col) <- rownames(coldata_sub)

pdf("v2/heatmap2.pdf")

  pheatmap(
    mat_z,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = FALSE,
    show_colnames = FALSE,
    annotation_col = annotation_col,
    color = colorRampPalette(c("blue", "white", "red"))(100),
    fontsize_col = 10,
    main = "Global expression heatmap (top 2,000 variable genes)"
  )

dev.off()

#########
# GSEA 
  # What biological programs are changing?

library(msigdbr) 
library(fgsea)

ranks <- resLFC %>%
  filter(!is.na(log2FoldChange)) %>%
  arrange(desc(log2FoldChange)) %>%
  dplyr::select(gene, log2FoldChange)

# Named numeric vector (required by fgsea)
ranked_genes <- ranks$log2FoldChange
names(ranked_genes) <- ranks$gene
names(ranked_genes) <- sub("\\..*$", "", names(ranked_genes))

hallmark_sets <- msigdbr(
  species  = "Homo sapiens",
  category = "H"
) %>%
  split(x = .$ensembl_gene, f = .$gs_name)

fgsea_res <- fgsea(
  pathways = hallmark_sets,
  stats    = ranked_genes
)

fgsea_res[order(fgsea_res$padj), ][, c("pathway", "NES", "padj")][1:10, ]

top_fgsea <- fgsea_res %>%
  filter(padj < 0.1)

pdf("GSEA.pdf")

  ggplot(top_fgsea,
         aes(x = reorder(pathway, NES), y = NES, fill = NES > 0)) +
    geom_col() +
    coord_flip() +
    scale_fill_manual(values = c("steelblue", "firebrick")) +
    labs(
      x = "Hallmark pathway",
      y = "Normalized Enrichment Score (NES)",
      title = "GSEA of early AraC response (1 day)"
    ) +
    theme_classic()

  plotEnrichment(
    hallmark_sets[["HALLMARK_INFLAMMATORY_RESPONSE"]],
    ranked_genes
  ) +
    labs(title = "GSEA: Inflammatory Response")

  plotEnrichment(
    hallmark_sets[["HALLMARK_MYC_TARGETS_V1"]],
    ranked_genes
  ) +
    labs(title = "GSEA: MYC Targets V1")

  plotEnrichment(
    hallmark_sets[["HALLMARK_MYC_TARGETS_V2"]],
    ranked_genes
  ) +
    labs(title = "GSEA: MYC Targets V2")

dev.off()

########
# Other signatures 

myeloid_sets <- msigdbr(
  species  = "Homo sapiens",
  category = "C7"
) %>%
  filter(
    grepl(
      "MYELOID",
      gs_name
    )
  ) %>%
  split(x = .$ensembl_gene, f = .$gs_name)

c("GSE10325_BCELL_VS_MYELOID_DN", "GSE10325_BCELL_VS_MYELOID_UP", 
"GSE10325_CD4_TCELL_VS_MYELOID_DN", "GSE10325_CD4_TCELL_VS_MYELOID_UP"
)

pdf("v2/myeloid.pdf")

  plotEnrichment(
    myeloid_sets[["GSE10325_BCELL_VS_MYELOID_DN"]],
    ranked_genes
  ) +
    labs(title = "GSEA: GSE10325_BCELL_VS_MYELOID_DN")

  plotEnrichment(
    myeloid_sets[["GSE10325_BCELL_VS_MYELOID_UP"]],
    ranked_genes
  ) +
    labs(title = "GSEA: GSE10325_BCELL_VS_MYELOID_UP")

  plotEnrichment(
    myeloid_sets[["GSE10325_CD4_TCELL_VS_MYELOID_DN"]],
    ranked_genes
  ) +
    labs(title = "GSEA: GSE10325_CD4_TCELL_VS_MYELOID_DN")

  plotEnrichment(
    myeloid_sets[["GSE10325_CD4_TCELL_VS_MYELOID_UP"]],
    ranked_genes
  ) +
    labs(title = "GSEA: GSE10325_CD4_TCELL_VS_MYELOID_UP")

dev.off()

lsc <- msigdbr(
  species  = "Homo sapiens",
  category = "C2"
) %>%
  filter(
    grepl(
      "LSC",
      gs_name
    )
  ) %>%
  split(x = .$ensembl_gene, f = .$gs_name)

c("EPPERT_CE_HSC_LSC", "EPPERT_LSC_R", "HUANG_AML_LSC47")

pdf("v2/LSC.pdf")

  plotEnrichment(
    lsc[["EPPERT_CE_HSC_LSC"]],
    ranked_genes
  ) +
    labs(title = "GSEA: EPPERT_CE_HSC_LSC")

  plotEnrichment(
    lsc[["EPPERT_LSC_R"]],
    ranked_genes
  ) +
    labs(title = "GSEA: EPPERT_LSC_R")

  plotEnrichment(
    lsc[["HUANG_AML_LSC47"]],
    ranked_genes
  ) +
    labs(title = "GSEA: HUANG_AML_LSC47")

dev.off()

hsc <- msigdbr(
  species  = "Homo sapiens",
  category = "C7"
) %>%
  filter(
    grepl(
      "HSC|HEMATOPOIETIC_STEM",
      gs_name
    )
  ) %>%
  split(x = .$ensembl_gene, f = .$gs_name)

fgsea_res <- fgsea(
  pathways = hsc,
  stats    = ranked_genes
)

fgsea_res[order(fgsea_res$padj), ][, c("pathway", "NES", "padj")][1:10, ]

top_fgsea <- fgsea_res %>%
  filter(padj < 0.1)

pdf("v2/HSC.pdf", width = 15)

  ggplot(top_fgsea,
         aes(x = reorder(pathway, NES), y = NES, fill = NES > 0)) +
    geom_col() +
    coord_flip() +
    scale_fill_manual(values = c("steelblue", "firebrick")) +
    labs(
      x = "Hallmark pathway",
      y = "Normalized Enrichment Score (NES)",
      title = "GSEA of early AraC response (1 day)"
    ) +
    theme_classic()

dev.off()

#### Goichi GSEA 
msig_go_bp <- msigdbr(
  species = "Homo sapiens",
  category = "C5"
) %>%
  split(x = .$ensembl_gene, f = .$gs_name)

names(msig_go_bp) %>% head()

names(msig_go_bp)[grep("DNA_BINDING_TRANSCRIPTION_FACTOR_ACTIVITY", names(msig_go_bp))]


pdf("v2/goichi.pdf")

  plotEnrichment(
    msig_go_bp[["GOBP_REGULATION_OF_DNA_BINDING_TRANSCRIPTION_FACTOR_ACTIVITY"]],
    ranked_genes
  ) +
    labs(title = "GSEA: Regulation of DNA Binding TF Activity")

  plotEnrichment(
    msig_go_bp[["GOBP_POSITIVE_REGULATION_OF_DNA_BINDING_TRANSCRIPTION_FACTOR_ACTIVITY"]],
    ranked_genes
  ) +
    labs(title = "GSEA: Positive Regulation of DNA Binding TF Activity")

dev.off()

msig_tf <- msigdbr(
  species = "Homo sapiens",
  category = "C3"
) %>%
  split(x = .$ensembl_gene, f = .$gs_name)

names(msig_tf) %>% head

names(msig_tf)[grep("GATA", names(msig_tf))]

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

names(msig_go_bp)[grep("FOXO4", names(msig_go_bp))]

sets <- c(
  "PU1" = "PU1_Q6",
  "GATA1" = "GATA1_01",
  # "RUNX1",
  "CEBPA" = "CEBPA_01",
  "MYCv1" = "HALLMARK_MYC_TARGETS_V1",
  "MYCv2" = "HALLMARK_MYC_TARGETS_V2",
  "FOXO4" = "FOXO4_TARGET_GENES"
)


runx_targets <- fread("v2/RUNX1_targets_JASPAR.csv")

symbol_df <- enframe(symbol_map, name = "ENSG", value = "Symbol") %>%
  filter(!is.na(Symbol))

# Join to your target list
runx_targets_ensg <- runx_targets %>%
  left_join(symbol_df, by = "Symbol")

runx_targets_ensg

as.vector(na.omit(runx_targets_ensg$ENSG))

pdf("v2/goichi_TFs.pdf")

  plotEnrichment(
    TF_sets[[sets[1]]],
    ranked_genes
  ) +
    labs(title = "GSEA: PU1 Targets")

  plotEnrichment(
    TF_sets[[sets[2]]],
    ranked_genes
  ) +
    labs(title = "GSEA: GATA1 Targets")

  plotEnrichment(
    TF_sets[[sets[3]]],
    ranked_genes
  ) +
    labs(title = "GSEA: CEBPA Targets")

  plotEnrichment(
    TF_sets[[sets[6]]],
    ranked_genes
  ) +
    labs(title = "GSEA: FOXO4 Targets")

  plotEnrichment(
    hallmark_sets[["HALLMARK_MYC_TARGETS_V1"]],
    ranked_genes
  ) +
    labs(title = "GSEA: MYC V1 Targets")

  plotEnrichment(
    hallmark_sets[["HALLMARK_MYC_TARGETS_V2"]],
    ranked_genes
  ) +
    labs(title = "GSEA: MYC V2 Targets")

  plotEnrichment(
    as.vector(na.omit(runx_targets_ensg$ENSG)),
    ranked_genes
  ) +
    labs(title = "JASPAR: RUNX1 Targets")

dev.off()







