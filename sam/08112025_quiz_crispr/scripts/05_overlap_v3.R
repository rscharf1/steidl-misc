library(data.table)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(patchwork)

# Take Sam's list of DNA Repair genes 
	# Plot where they fall in Sam's volcano plot 
	# Take the winners, and look at where they fall in the paper's plots 

repair <- fread("inputs/DNA_repair_genes_mouse.csv")

# Read in Sam's data
sam <- readRDS("intermediates/sams_data_clean.rds") # neg logFC means drug caused downregulation

table(repair$gene_name %in% sam$gene_name_upper)

sam$dna_repair <- sam$gene_name_upper %in% repair$gene_name

# Volcano plot
pdf("outputs/sam_volc.pdf")
	ggplot() + 
		geom_point(
			data = sam[significant == TRUE], 
			aes(x = log2FoldChange, y = logP),
			color = "firebrick3"
		) + 
		geom_point(
			data = sam[significant == FALSE], 
			aes(x = log2FoldChange, y = logP),
			color = "grey70"
		) + 
		labs(x = "Log2(Fold Change)", y = "-log(P)", title = "Quizartinib v. Vehicle")
dev.off()

sam$group <- "Not significant"
sam$group[sam$significant] <- "Significant"
sam$group[sam$dna_repair] <- "DNA repair"

pdf("outputs/volc_dna_repair2.pdf")
	ggplot() +
	  geom_point(
	    data = sam[group == "Not significant"],
	    aes(x = log2FoldChange, y = logP, color = group)
	  ) +
	  geom_point(
	    data = sam[group == "Significant"],
	    aes(x = log2FoldChange, y = logP, color = group)
	  ) +
	  geom_point(  # plotted last = on top
	    data = sam[group == "DNA repair"],
	    aes(x = log2FoldChange, y = logP, color = group)
	  ) +
	  scale_color_manual(
	    values = c(
	      "Not significant" = "grey70",
	      "Significant" = "lightskyblue2",
	      "DNA repair" = "firebrick3"
	    )
	  ) +
	  labs(
	    x = "Log2(Fold Change)",
	    y = "-log(P)",
	    title = "Quizartinib v. Vehicle",
	    color = "Category"
	  )
dev.off()

paper <- readRDS("intermediates/paper_drugs_genes.rds")

sam$in_screen <- sam$gene_name_upper %in% unique(paper$gene)

pdf("outputs/volc_dna_repair3.pdf")
	ggplot() +
	  geom_point(
	    data = sam[in_screen == TRUE][group == "Not significant"],
	    aes(x = log2FoldChange, y = logP, color = group)
	  ) +
	  geom_point(
	    data = sam[in_screen == TRUE][group == "Significant"],
	    aes(x = log2FoldChange, y = logP, color = group)
	  ) +
	  geom_point(  # plotted last = on top
	    data = sam[in_screen == TRUE][group == "DNA repair"],
	    aes(x = log2FoldChange, y = logP, color = group)
	  ) +
	  scale_color_manual(
	    values = c(
	      "Not significant" = "grey70",
	      "Significant" = "lightskyblue2",
	      "DNA repair" = "firebrick3"
	    )
	  ) +
	  labs(
	    x = "Log2(Fold Change)",
	    y = "-log(P)",
	    title = "Quizartinib v. Vehicle",
	    color = "Category"
	  )
dev.off()





sam$to_label <- FALSE
# sam[group == "DNA repair"][log2FoldChange < 0][order(-logP)][1:10]$to_label <- TRUE

# sam[group == "DNA repair"][log2FoldChange < 0][logP > 5][gene_name_upper %in% unique(paper$gene)]$to_label <- TRUE
sam[group == "DNA repair"][log2FoldChange < 0][significant == TRUE][gene_name_upper %in% unique(paper$gene)]$to_label <- TRUE


pdf("outputs/volc_dna_repair3.pdf")
	ggplot() +
	  geom_point(
	    data = sam[in_screen == TRUE][group == "Not significant"],
	    aes(x = log2FoldChange, y = logP, color = group)
	  ) +
	  geom_point(
	    data = sam[in_screen == TRUE][group == "Significant"],
	    aes(x = log2FoldChange, y = logP, color = group)
	  ) +
	  geom_point(  # plotted last = on top
	    data = sam[in_screen == TRUE][group == "DNA repair"],
	    aes(x = log2FoldChange, y = logP, color = group)
	  ) +
	  scale_color_manual(
	    values = c(
	      "Not significant" = "grey70",
	      "Significant" = "lightskyblue2",
	      "DNA repair" = "firebrick3"
	    )
	  ) +
      geom_label_repel(
		  data = sam[in_screen == TRUE][to_label == TRUE],
		  aes(x = log2FoldChange, y = logP, label = gene_name_upper),
		  size = 4,
		  min.segment.length = 0,
		  box.padding = 0.25,
		  point.padding = 0.15,
		  max.overlaps = Inf,
		  segment.size = 0.4
		) +
	  labs(
	    x = "Log2(Fold Change)",
	    y = "-log(P)",
	    title = "Quizartinib v. Vehicle",
	    color = "Category"
	  )
dev.off()

# Read in drug data
	# Create the volcano plot, highlighting these 10 genes on each 

paper <- readRDS("intermediates/paper_drugs_genes.rds")

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

	to_label <- sam[to_label == TRUE]$gene_name_upper

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
	  # geom_text_repel(
	  #   data = paper_sub[is_label == TRUE],
	  #   aes(x = logFC, y = logP, label = gene),
	  #   size = 4,
	  #   nudge_x = ifelse(paper_sub[is_label == TRUE]$logFC < 0, -0.015, 0.015),
	  #   min.segment.length = 0,
	  #   box.padding = 0.25,
	  #   point.padding = 0.15,
	  #   max.overlaps = Inf,
	  #   segment.size = 0.4
	  # ) +
    	geom_label_repel(
		  data = paper_sub[is_label == TRUE],
		  aes(x = logFC, y = logP, label = gene),
		  size = 4,
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





