library(dplyr)
library(data.table)
library(ggplot2)

dt <- fread("diploidSV_table.tsv")

colnames(dt) <- c(
  "chr", "pos", "svtype", "filter",
  "HL60_ctrl", "PY_low", "PY_high"
)

dt_long <- melt(
  dt,
  id.vars = c("chr", "pos", "svtype", "filter"),
  variable.name = "sample",
  value.name = "gt"
)

# Plot A
panelA <- dt_long[
  filter == "PASS" & gt != "0/0",
  .N,
  by = .(sample, svtype)
]

panelA <- merge(
  CJ(
    sample = c("HL60_ctrl", "PY_low", "PY_high"),
    svtype = c("DEL", "DUP", "INV", "BND")
  ),
  panelA,
  by = c("sample", "svtype"),
  all.x = TRUE
)

panelA[is.na(N), N := 0]

panelA[, svtype_label := fifelse(
  svtype == "BND",
  "One-sided breakpoints (BND)",
  svtype
)]

panelA[, svtype_label := factor(
  svtype_label,
  levels = c(
    "DEL",
    "DUP",
    "INV",
    "One-sided breakpoints (BND)"
  )
)]

panelA[, svtype_full := fcase(
  svtype == "DEL", "Deletion",
  svtype == "DUP", "Duplication",
  svtype == "INV", "Inversion",
  svtype == "BND", "One-sided breakpoints"
)]

pdf("plotA.pdf")

	ggplot(panelA, aes(x = sample, y = N, fill = svtype)) +
	  geom_col(color = "black", width = 0.7) +
	  labs(
	    x = NULL,
	    y = "SV count",
	    fill = "SV type",
	    title = "High-confidence SV calls are dominated by one-sided breakpoints"
	  ) +
	  theme_classic(base_size = 14)

	ggplot(panelA, aes(x = sample, y = N, fill = svtype_full)) +
	  geom_col(color = "black", width = 0.7) +
	  labs(
	    x = NULL,
	    y = "SV count (PASS)",
	    fill = "SV class",
	    title = "SV calls are dominated by one-sided breakpoints"
	  ) +
	  theme_classic(base_size = 14) +
	  theme(
	    legend.position = "right",
	    plot.title = element_text(face = "bold")
	  )

	grid::grid.newpage()

dev.off()

# Plot B
panelB <- dt_long[
  filter != "PASS" & gt != "0/0"
][
  ,
  .(filter_reason = unlist(strsplit(filter, ";"))),
  by = .(chr, pos, svtype, sample)
][
  ,
  .N,
  by = .(sample, filter_reason)
]

panelB[, filter_reason := factor(
  filter_reason,
  levels = panelB[, .(total = sum(N)), by = filter_reason][
    order(-total), filter_reason
  ]
)]

pdf("plotB.pdf")

	ggplot(panelB, aes(x = sample, y = N, fill = filter_reason)) +
	  geom_col(color = "black", width = 0.7) +
	  labs(
	    x = NULL,
	    y = "Filtered SV count",
	    fill = "Filter reason",
	    title = "Low-confidence candidate SVs failing quality filters"
	  ) +
	  theme_classic(base_size = 14)

dev.off()

##########
dt_long[, confidence := fifelse(filter == "PASS", "High confidence (PASS)", "Low confidence (filtered)")]

# count SVs present per sample
panel <- dt_long[
  gt != "0/0",
  .N,
  by = .(sample, confidence)
]

# ensure both confidence classes appear
panel <- merge(
  CJ(
    sample = c("HL60_ctrl", "PY_low", "PY_high"),
    confidence = c("Low confidence (filtered)", "High confidence (PASS)")
  ),
  panel,
  by = c("sample", "confidence"),
  all.x = TRUE
)

panel[is.na(N), N := 0]

pdf("plotC.pdf")

	ggplot(panel, aes(x = sample, y = N, fill = confidence)) +
	  geom_col(color = "black", width = 0.7) +
	  scale_fill_manual(
	    values = c(
	      "Low confidence (filtered)" = "grey80",
	      "High confidence (PASS)"    = "#D55E00"
	    )
	  ) +
	  labs(
	    x = NULL,
	    y = "Structural variant count",
	    fill = NULL
	  ) +
	  theme_classic(base_size = 14)

dev.off()

