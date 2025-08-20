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

# For each drug:
	# Grab top 20 sensitizing knockouts
	# Show the volcano plot with them highlighted 

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

	# distance to -0.2, 4

	paper_sub$dist <- sqrt((paper_sub$logFC - -0.2)**2 + (paper_sub$logP - 4)**2)

	scaled <- paper_sub %>%
	  mutate(
	    logFC_scaled = (logFC - min(logFC)) / (max(logFC) - min(logFC)),
	    logP_scaled  = (logP  - min(logP))  / (max(logP)  - min(logP))
	  )

	corner <- c(min(scaled$logFC_scaled), max(scaled$logP_scaled))

	scaled <- scaled %>%
	  mutate(dist = sqrt((logFC_scaled - corner[1])^2 +
	                     (logP_scaled  - corner[2])^2))

	scaled[order(dist)]

	to_label <- scaled[order(dist)][1:20]$gene

	paper_sub[, is_label := gene %in% to_label]

	x_max <- paper_sub$logFC %>% range %>% abs %>% max %>% {. + 0.1} %>% { ceiling(. / 0.1) * 0.1 }

	xlim <- c(x_max * -1, x_max)
	ylim <- c(0, 6)

	p1 <- ggplot() +
	  geom_point(
	    data = paper_sub[is_label == FALSE],
	    aes(x = logFC, y = logP),
	    color = "grey70", size = 1.4, alpha = 0.9, shape = 16
	  ) +
	  geom_point(
	    data = paper_sub[is_label == TRUE ],
	    aes(x = logFC, y = logP),
	    color = "red2", size = 3, shape = 16
	  ) +
	  geom_segment(aes(x = 0, xend = 0, y = 0, yend = ylim[2]), linewidth = 1) +
	  theme_classic(base_size = 14) +
	  labs(x = "LogFC", y = "-log(P)", title = my_drug) + 
	  xlim(xlim)


	 p2 <- ggplot() + 
		geom_point(
			data = sam[significant == FALSE], 
			aes(x = log2FoldChange, y = logP),
			color = "grey70"
		) + 
		geom_point(
			data = sam[significant == TRUE], 
			aes(x = log2FoldChange, y = logP),
			color = "lightskyblue2"
		) + 
		geom_point(
			data = sam[gene_name_upper %in% to_label], 
			aes(x = log2FoldChange, y = logP),
			color = "firebrick3"
		) +
		labs(x = "Log2(Fold Change)", y = "-log(P)", title = "Sam: Quizartinib v. Vehicle")

		list(p1, p2)
})

pdf("outputs/test.pdf", width = 10)
	for(i in 1:8) {
		print(plots[[i]][[1]] + plots[[i]][[2]])
	}
dev.off()








