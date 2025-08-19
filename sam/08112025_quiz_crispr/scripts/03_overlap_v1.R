library(data.table)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(patchwork)
library(forcats)

# Read in paper data
paper <- readRDS("intermediates/paper_drugs_genes.rds")

# Read in Sam's data
sam <- readRDS("intermediates/sams_data_clean.rds") # neg logFC means drug caused downregulation

# 
sam[significant == TRUE][abs(log2FoldChange) > 3][order(pvalue)][1:20]

drugs <- c(
	"CDDP" = "Cisplatin",
	"pm" = "Phosphoramide Mustard",
	"Doxorubicin",
	"Etoposide",
	"Topotecan",
	"Vincristine",
	"ATRA",
	"JQAD1"
)

plots <- lapply(seq_along(drugs), function(x) {
	my_drug <- drugs[x]
	paper_sub <- paper[drug == my_drug]

	to_label <- sam[significant == TRUE][gene_name_upper %in% paper_sub$gene][abs(log2FoldChange) > 1][order(pvalue)][1:20]$gene_name_upper

	paper_sub[, is_label := gene %in% to_label]
	paper_sub[, side := ifelse(logFC < 0, "left", "right")]

	x_max <- paper_sub$logFC %>% range %>% abs %>% max %>% {. + 0.1} %>% { ceiling(. / 0.1) * 0.1 }

	xlim <- c(x_max * -1, x_max)
	ylim <- c(0, 6)

	p <- ggplot() +
	  # background points
	  geom_point(
	    data = paper_sub[is_label == FALSE & side == "left"],
	    aes(x = logFC, y = logP),
	    color = "grey70", size = 1.4, alpha = 0.9, shape = 16
	  ) +
	  geom_point(
	    data = paper_sub[is_label == FALSE & side == "right"],
	    aes(x = logFC, y = logP),
	    color = "grey70", size = 1.4, alpha = 0.9, shape = 17
	  ) +
	  # highlighted (left = red circles)
	  geom_point(
	    data = paper_sub[is_label == TRUE & side == "left"],
	    aes(x = logFC, y = logP),
	    color = "red2", size = 3, shape = 16
	  ) +
	  # highlighted (right = blue triangles)
	  geom_point(
	    data = paper_sub[is_label == TRUE & side == "right"],
	    aes(x = logFC, y = logP),
	    color = "steelblue3", size = 3, shape = 17
	  ) +
	  # labels
	  geom_text_repel(
	    data = paper_sub[is_label == TRUE],
	    aes(x = logFC, y = logP, label = gene),
	    size = 4,
	    nudge_x = ifelse(paper_sub[is_label == TRUE]$logFC < 0, -0.015, 0.015),
	    min.segment.length = 0,
	    box.padding = 0.25,
	    point.padding = 0.15,
	    max.overlaps = Inf,
	    segment.size = 0.4
	  ) +
	  geom_segment(aes(x = 0, xend = 0, y = 0, yend = ylim[2]), linewidth = 1) +
	  theme_classic(base_size = 14) +
	  labs(x = "LogFC", y = "-log(P)", title = my_drug) + 
	  xlim(xlim)

	  p

})

pdf("outputs/test.pdf", width = 16, height = 14)
	(plots[[1]] + plots[[2]]) / 
	(plots[[3]] + plots[[4]]) / 
	(plots[[5]] + plots[[6]]) / 
	(plots[[7]] + plots[[8]]) 
dev.off()

# Genes that are:
	# Significantly down-regulated in Sams data
	# Repeatedely show up in the "sensitivity" section of the papers volcano plots 

sig <- sam[significant == TRUE][log2FoldChange < -1][gene_name_upper %in% unique(paper$gene)][order(pvalue)]

res <- lapply(seq_len(nrow(sig)), function(x) {
	my_gene <- sig$gene_name_upper[x]

	data.table(
		gene = my_gene,
		n = nrow(paper[gene == my_gene][logFC < 0])
	)
}) %>% rbindlist()

pdf("outputs/gene_freq.pdf")
	res %>%
	  mutate(gene = fct_reorder(gene, n, .desc = TRUE)) %>%
	  ggplot(., aes(x = gene, y = n, fill = n)) +
	  geom_col(width = 0.7, show.legend = FALSE) +
	  scale_fill_gradient(low = "#A6CEE3", high = "#1F78B4") +
	  labs(
	    x = NULL, 
	    y = "Sensitization in n of 8 drugs", 
	    title = "Gene Sensitization Frequency"
	  ) +
	  theme_minimal(base_size = 13) +
	  theme(
	    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1.65),
	    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
	    panel.grid.major.x = element_blank(),
	    panel.grid.minor = element_blank()
	  )
dev.off()

to_label <- res[n >= 7]$gene

plots <- lapply(seq_along(drugs), function(x) {
	my_drug <- drugs[x]
	paper_sub <- paper[drug == my_drug]

	paper_sub[, is_label := gene %in% to_label]
	paper_sub[, side := ifelse(logFC < 0, "left", "right")]

	x_max <- paper_sub$logFC %>% range %>% abs %>% max %>% {. + 0.1} %>% { ceiling(. / 0.1) * 0.1 }

	xlim <- c(x_max * -1, x_max)
	ylim <- c(0, 6)

	p <- ggplot() +
	  # background points
	  geom_point(
	    data = paper_sub[is_label == FALSE & side == "left"],
	    aes(x = logFC, y = logP),
	    color = "grey70", size = 1.4, alpha = 0.9, shape = 16
	  ) +
	  geom_point(
	    data = paper_sub[is_label == FALSE & side == "right"],
	    aes(x = logFC, y = logP),
	    color = "grey70", size = 1.4, alpha = 0.9, shape = 17
	  ) +
	  # highlighted (left = red circles)
	  geom_point(
	    data = paper_sub[is_label == TRUE & side == "left"],
	    aes(x = logFC, y = logP),
	    color = "red2", size = 3, shape = 16
	  ) +
	  # highlighted (right = blue triangles)
	  geom_point(
	    data = paper_sub[is_label == TRUE & side == "right"],
	    aes(x = logFC, y = logP),
	    color = "steelblue3", size = 3, shape = 17
	  ) +
	  # labels
	  geom_text_repel(
	    data = paper_sub[is_label == TRUE],
	    aes(x = logFC, y = logP, label = gene),
	    size = 4,
	    nudge_x = ifelse(paper_sub[is_label == TRUE]$logFC < 0, -0.015, 0.015),
	    min.segment.length = 0,
	    box.padding = 0.25,
	    point.padding = 0.15,
	    max.overlaps = Inf,
	    segment.size = 0.4
	  ) +
	  geom_segment(aes(x = 0, xend = 0, y = 0, yend = ylim[2]), linewidth = 1) +
	  theme_classic(base_size = 14) +
	  labs(x = "LogFC", y = "-log(P)", title = my_drug) + 
	  xlim(xlim)

	  p

})

pdf("outputs/test.pdf", width = 16, height = 14)
	(plots[[1]] + plots[[2]]) / 
	(plots[[3]] + plots[[4]]) / 
	(plots[[5]] + plots[[6]]) / 
	(plots[[7]] + plots[[8]]) 
dev.off()

pdf("outputs/sam_volc.pdf")
	ggplot() + 
		geom_point(
			data = sam[significant == FALSE], 
			aes(x = log2FoldChange, y = logP),
			color = "grey70"
		) + 
		geom_point(
			data = sam[significant == TRUE & !(gene_name_upper %in% to_label)], 
			aes(x = log2FoldChange, y = logP),
			color = "lightskyblue2"
		) + 
		geom_point(
			data = sam[significant == TRUE & gene_name_upper %in% to_label], 
			aes(x = log2FoldChange, y = logP),
			color = "firebrick3"
		) +
		geom_text_repel(
	      data = sam[gene_name_upper %in% to_label],
	      aes(x = log2FoldChange, y = logP, label = gene_name_upper),
	      size = 4,
	      # nudge_x = ifelse(dt_m[!is.na(pathway)]$log2FoldChange < 0, -0.015, 0.015),
	      min.segment.length = 0,
	      box.padding = 0.25,
	      point.padding = 0.15,
	      max.overlaps = Inf,
	      segment.size = 0.4
	    ) + 
		labs(x = "Log2(Fold Change)", y = "-log(P)", title = "Quizartinib v. Vehicle")
dev.off()




