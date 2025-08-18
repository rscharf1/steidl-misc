library(data.table)
library(ggplot2)
library(dplyr)
library(readxl)
library(ggrepel)
library(patchwork)

# Recapitulate data from the paper
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

drug_abbr_table <- data.table(
	abbr = names(drugs),
	drug = drugs
)

drug_abbr_table$abbr <- ifelse(nchar(drug_abbr_table$abbr) > 0, drug_abbr_table$abbr, drug_abbr_table$drug)

# Assemble gene by drug table 
dt <- lapply(seq_len(nrow(drug_abbr_table)), function(x) {
	drug <- drug_abbr_table$drug[x]
	dt <- read_excel("inputs/screen_data.xlsx", sheet = drug_abbr_table$abbr[x]) %>% as.data.table()

	# dt <- read_excel("inputs/screen_data.xlsx", sheet = "Etoposide") %>% as.data.table()

	dt <- dt[, c(13,14,15,16), with = FALSE]
	colnames(dt) <- c("gene", "logFC", "logP_neg", "logP_pos")

	dt$logP <- ifelse(dt$logFC < 0, dt$logP_neg, dt$logP_pos)

	dt <- dt[, c(1,2,5), with = FALSE]

	dt$drug <- drug

	dt	
}) %>% rbindlist()

saveRDS(dt, "intermediates/paper_drugs_genes.rds")

# Look at one drug first
dt_sub <- dt[drug == "Etoposide"]

pdf("outputs/volcano.pdf")
	ggplot(dt_sub, aes(x = logFC, y = logP)) + 
		geom_point() + 
		xlim(-0.3, 0.3) + 
		labs(x = "LogFC", y = "-log(P)", title = "Etoposide")
dev.off()

	# Label some of those points
to_label <- c(
	"CFTR",
	"HDAC9",
	"MET",
	"MAP3K12",
	"RARA",
	"CDK4",
	"HDAC3",
	"TOP2A",
	"AMY2A"
)

dt_sub[, is_label := gene %in% to_label]
dt_sub[, side := ifelse(logFC < 0, "left", "right")]

xlim <- c(-0.3, 0.3)
ylim <- c(0, 6)

pdf("outputs/volcano2.pdf")

ggplot() +
  # background points
  geom_point(
    data = dt_sub[is_label == FALSE & side == "left"],
    aes(x = logFC, y = logP),
    color = "grey70", size = 1.4, alpha = 0.9, shape = 16
  ) +
  geom_point(
    data = dt_sub[is_label == FALSE & side == "right"],
    aes(x = logFC, y = logP),
    color = "grey70", size = 1.4, alpha = 0.9, shape = 17
  ) +
  # highlighted (left = red circles)
  geom_point(
    data = dt_sub[is_label == TRUE & side == "left"],
    aes(x = logFC, y = logP),
    color = "red2", size = 3, shape = 16
  ) +
  # highlighted (right = blue triangles)
  geom_point(
    data = dt_sub[is_label == TRUE & side == "right"],
    aes(x = logFC, y = logP),
    color = "steelblue3", size = 3, shape = 17
  ) +
  # labels
  geom_text_repel(
    data = dt_sub[is_label == TRUE],
    aes(x = logFC, y = logP, label = gene),
    size = 4,
    nudge_x = ifelse(dt_sub[is_label == TRUE]$logFC < 0, -0.015, 0.015),
    min.segment.length = 0,
    box.padding = 0.25,
    point.padding = 0.15,
    max.overlaps = Inf,
    segment.size = 0.4
  ) +
  # thick axes at 0
  # geom_segment(aes(x = xlim[1], xend = xlim[2], y = 0, yend = 0), linewidth = 1) +
  geom_segment(aes(x = 0, xend = 0, y = 0, yend = ylim[2]), linewidth = 1) +
  # coord_cartesian(xlim = xlim, ylim = ylim, clip = "off") +
  theme_classic(base_size = 14) +
  # theme(
  #   plot.margin = margin(10, 20, 10, 30),   # extra room so leftmost label isn’t clipped
  #   axis.line = element_blank()             # we drew custom axes above
  # ) +
  labs(x = "LogFC", y = "-log(P)", title = "Etoposide") + 
  xlim(xlim)

dev.off()

# 4 drugs with labeled points to show recapitulation
genes_to_label <- fread("inputs/genes_to_label.csv")

drugs <- genes_to_label$Drug %>% unique()

plots <- lapply(seq_along(drugs), function(x) {
	my_drug <- drugs[x]
	dt_sub <- dt[drug == my_drug]

	to_label <- genes_to_label[Drug == my_drug]$Gene

	dt_sub[, is_label := gene %in% to_label]
	dt_sub[, side := ifelse(logFC < 0, "left", "right")]

	xlim <- c(-0.3, 0.3)
	ylim <- c(0, 6)

	ggplot() +
	  # background points
	  geom_point(
	    data = dt_sub[is_label == FALSE & side == "left"],
	    aes(x = logFC, y = logP),
	    color = "grey70", size = 1.4, alpha = 0.9, shape = 16
	  ) +
	  geom_point(
	    data = dt_sub[is_label == FALSE & side == "right"],
	    aes(x = logFC, y = logP),
	    color = "grey70", size = 1.4, alpha = 0.9, shape = 17
	  ) +
	  # highlighted (left = red circles)
	  geom_point(
	    data = dt_sub[is_label == TRUE & side == "left"],
	    aes(x = logFC, y = logP),
	    color = "red2", size = 3, shape = 16
	  ) +
	  # highlighted (right = blue triangles)
	  geom_point(
	    data = dt_sub[is_label == TRUE & side == "right"],
	    aes(x = logFC, y = logP),
	    color = "steelblue3", size = 3, shape = 17
	  ) +
	  # labels
	  geom_text_repel(
	    data = dt_sub[is_label == TRUE],
	    aes(x = logFC, y = logP, label = gene),
	    size = 4,
	    nudge_x = ifelse(dt_sub[is_label == TRUE]$logFC < 0, -0.015, 0.015),
	    min.segment.length = 0,
	    box.padding = 0.25,
	    point.padding = 0.15,
	    max.overlaps = Inf,
	    segment.size = 0.4
	  ) +
	  # thick axes at 0
	  # geom_segment(aes(x = xlim[1], xend = xlim[2], y = 0, yend = 0), linewidth = 1) +
	  geom_segment(aes(x = 0, xend = 0, y = 0, yend = ylim[2]), linewidth = 1) +
	  # coord_cartesian(xlim = xlim, ylim = ylim, clip = "off") +
	  theme_classic(base_size = 14) +
	  # theme(
	  #   plot.margin = margin(10, 20, 10, 30),   # extra room so leftmost label isn’t clipped
	  #   axis.line = element_blank()             # we drew custom axes above
	  # ) +
	  labs(x = "LogFC", y = "-log(P)", title = my_drug) + 
	  xlim(xlim)

})

pdf("outputs/volcanos.pdf", width = 12, height = 12)
	(plots[[1]] + plots[[2]]) / 
	(plots[[3]] + plots[[4]])
dev.off()


