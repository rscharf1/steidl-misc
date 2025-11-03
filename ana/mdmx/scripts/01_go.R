library(dplyr)
library(data.table)
library(ggplot2)
library(ggrepel)

# Columns 
	# baseMean: avgnormalized expression across all samples
	# lfcSE: standard error ofFC
	# stat: Wald test stat
	# padj: use <0.05

# WT v. MDMX-Tg

dt <- fread("inputs/wt_v_mdmx.csv", skip = 2, header = TRUE)

dt[, neglog10_padj := -log10(padj)]

pdf("out.pdf")

dt %>% 
	ggplot(., aes(x = log2FoldChange, y = neglog10_padj)) + 
	geom_point(alpha = 0.7, size = 2) +
	geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed") +
	theme_minimal() +
	labs(
	x = "log2FC(MDMX-Tg/WT)",
	y = expression(-log[10](adjusted~p~value)),
	color = "Gene Biotype",
	title = "Volcano Plot Colored by Gene Biotype"
	) +
	theme(
	plot.title = element_text(hjust = 0.5),
	legend.position = "right"
	)

dt %>%
	ggplot(aes(x = gene_biotype)) +
	geom_bar(fill = "steelblue") +
	theme_minimal() +
	labs(
	x = "Gene Biotype",
	y = "Number of Genes",
	title = "Distribution of Genes by Biotype"
	) +
	theme(
	axis.text.x = element_text(angle = 45, hjust = 1),
	plot.title = element_text(hjust = 0.5)
	)

dt %>%
  ggplot(aes(x = log2FoldChange, y = neglog10_padj, color = gene_biotype)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed") +
  theme_minimal() +
  labs(
    x = "log2FC(MDMX-Tg/WT)",
    y = expression(-log[10](adjusted~p~value)),
    color = "Gene Biotype",
    title = "Volcano Plot Colored by Gene Biotype"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right"
  )

dt %>%
  ggplot(aes(x = log2FoldChange, y = neglog10_padj, color = gene_biotype)) +
  geom_point(alpha = 0.7, size = 1.8) +
  geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed") +
  facet_wrap(~ gene_biotype, scales = "free") +
  theme_minimal() +
  labs(
    x = "log2FC(MDMX-Tg/WT)",
    y = expression(-log[10](adjusted~p~value)),
    color = "Gene Biotype",
    title = "Volcano Plots Split by Gene Biotype"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",   # hide legend since each facet is labeled
    strip.text = element_text(face = "bold", size = 9)
  )

count_cols <- c("M351WT", "M352Mdmx", "M371Mdmx", "M380Mdmx", "M383WT", "M396WT")

dt[, n_samples_above0 := rowSums(.SD > 0), .SDcols = count_cols]

dt[gene_biotype == "protein_coding"][!is.na(neglog10_padj)][n_samples_above0 > 2] %>% 
	ggplot(., aes(x = log2FoldChange, y = neglog10_padj)) + 
	geom_point(alpha = 0.7, size = 2) +
	geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed") +
	theme_minimal() +
	labs(
	x = "log2FC(MDMX-Tg/WT)",
	y = expression(-log[10](adjusted~p~value)),
	color = "Gene Biotype",
	title = "Volcano Plot Protein Coding Genes"
	) +
	theme(
	plot.title = element_text(hjust = 0.5),
	legend.position = "right"
	)


dt_sub <- dt[gene_biotype == "protein_coding" &
             !is.na(neglog10_padj) &
             n_samples_above0 > 2]

# Then plot using that explicitly
ggplot(dt_sub, aes(x = log2FoldChange, y = neglog10_padj)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed") +
  geom_text_repel(
    data = dt_sub[abs(log2FoldChange) > 9 | neglog10_padj > 8],
    aes(label = mgi_symbol),
    size = 3,
    max.overlaps = 50,
    box.padding = 0.5,
    point.padding = 0.2,
    segment.color = "grey50"
  ) +
  theme_minimal() +
  labs(
    x = "log2FC(MDMX-Tg/WT)",
    y = expression(-log[10](adjusted~p~value)),
    title = "Volcano Plot Protein Coding Genes"
  ) +
  theme(plot.title = element_text(hjust = 0.5))

dev.off()

# Intersect with surface data 
surface <- fread("inputs/surface.csv")

surface[, gene_symbol := trimws(`ENTREZ gene symbol`)]
dt[, mgi_symbol := trimws(mgi_symbol)]

surface[, is_surface := TRUE]

surface_unique <- unique(surface[, .(gene_symbol = `ENTREZ gene symbol`,
                                     is_surface = TRUE,
                                     CSPA_category = `CSPA category`,
                                     UP_Protein_name,
                                     CD)])

m <- merge(
  dt,
  surface_unique,
  by.x = "mgi_symbol",
  by.y = "gene_symbol",
  all.x = TRUE
)

pdf("out2.pdf")

m[gene_biotype == "protein_coding"][!is.na(neglog10_padj)][n_samples_above0 > 2][is_surface == TRUE] %>% 
	ggplot(., aes(x = log2FoldChange, y = neglog10_padj)) + 
	geom_point(alpha = 0.7, size = 2) +
	geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed") +
	theme_minimal() +
	labs(
	x = "log2FC(MDMX-Tg/WT)",
	y = expression(-log[10](adjusted~p~value)),
	color = "Gene Biotype",
	title = "Volcano Plot Protein Coding Genes"
	) +
	theme(
	plot.title = element_text(hjust = 0.5),
	legend.position = "right"
	)


m_sub <- m[gene_biotype == "protein_coding" &
             !is.na(neglog10_padj) &
             n_samples_above0 > 2 & 
             is_surface == TRUE]

# Then plot using that explicitly
ggplot(m_sub, aes(x = log2FoldChange, y = neglog10_padj)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed") +
  geom_text_repel(
    # data = m_sub[abs(log2FoldChange) > 9 | neglog10_padj > 8],
    data = m_sub[neglog10_padj > -log10(0.05)][order(-abs(log2FoldChange))][1:20],
    aes(label = mgi_symbol),
    size = 3,
    max.overlaps = 50,
    box.padding = 0.5,
    point.padding = 0.2,
    segment.color = "grey50"
  ) +
  theme_minimal() +
  labs(
    x = "log2FC(MDMX-Tg/WT)",
    y = expression(-log[10](adjusted~p~value)),
    title = "Volcano Plot Protein Coding Genes"
  ) +
  theme(plot.title = element_text(hjust = 0.5))

dev.off()