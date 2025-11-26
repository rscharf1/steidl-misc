library(dplyr)
library(data.table)
library(ggplot2)
library(ggrepel)

wt <- fread("outputs/wt_v_mdmx.csv")
ure <- fread("outputs/ure_v_ure_mdmx.csv")

pdf("outputs/out5.pdf")

genes <- wt[neglog10_padj > -log10(0.05)]$mgi_symbol

ggplot(ure, aes(x = log2FoldChange, y = neglog10_padj)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed") +
  geom_text_repel(
    data = ure[mgi_symbol %in% genes],
    aes(label = mgi_symbol),
    size = 3,
    max.overlaps = 50,
    box.padding = 0.5,
    point.padding = 0.2,
    segment.color = "grey50"
  ) +
  theme_minimal() +
  labs(
    x = "log2FC(URE;MDMX-Tg/URE)",
    y = expression(-log[10](adjusted~p~value)),
    title = "Differentially Expressed Genes in MDMX-Tg v. WT"
  ) +
  theme(plot.title = element_text(hjust = 0.5))

genes <- ure[neglog10_padj > -log10(0.05)]$mgi_symbol

ggplot(wt, aes(x = log2FoldChange, y = neglog10_padj)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed") +
  geom_text_repel(
    data = wt[mgi_symbol %in% genes],
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
    title = "Differentially Expressed Genes in URE;MDMX-Tg/URE"
  ) +
  theme(plot.title = element_text(hjust = 0.5))

dev.off()

################
wt <- fread("inputs/wt_v_mdmx.csv", skip = 2, header = TRUE)

# mean of top 3
count_cols <- c("M351WT", "M352Mdmx", "M371Mdmx", "M380Mdmx", "M383WT", "M396WT")
wt[, top3_mean := apply(.SD, 1, function(x) mean(sort(x, decreasing = TRUE)[1:3])),
    .SDcols = count_cols]

wt_de <- fread("outputs/wt_v_mdmx.csv") %>% 
	.[neglog10_padj > -log10(0.05)] %>% 
	.$mgi_symbol

# Rank genes by top3_mean
wt_plot <- wt %>%
  arrange(top3_mean) %>%
  mutate(rank = row_number())

# Identify Esam row (if exists)
esam_row <- wt_plot %>% filter(mgi_symbol == "Esam")
tet2_row <- wt_plot %>% filter(mgi_symbol == "Tet2")

pdf("outputs/out6.pdf")

ggplot(wt_plot, aes(x = rank, y = top3_mean)) +
  geom_line(color = "gray40") +
  scale_y_log10() +
  labs(
    x = "Increasing Average Expression",
    y = "Mean Expression of Top 3 Samples",
    title = "Ranked gene expression (MDMX-Tg v. WT)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_blank()
  )

ggplot(wt_plot, aes(x = rank, y = top3_mean)) +
  geom_line(color = "gray40") +
  geom_point(data = esam_row, aes(x = rank, y = top3_mean),
             color = "red", size = 3) +
  geom_text_repel(
	  data = esam_row,
	  aes(x = rank, y = top3_mean, label = "Esam"),
	  color = "red",
	  size = 4,
	  nudge_y = 0.5,         # push upward
	  nudge_x = 5,           # push right slightly
	  box.padding = 0.5,
	  min.segment.length = 0,
	  max.overlaps = Inf
	) + 
  geom_point(data = tet2_row, aes(x = rank, y = top3_mean),
             color = "red", size = 3) +
  geom_text_repel(
	  data = tet2_row,
	  aes(x = rank, y = top3_mean, label = "TET2"),
	  color = "red",
	  size = 4,
	  nudge_y = 0.5,         # push upward
	  nudge_x = 5,           # push right slightly
	  box.padding = 0.5,
	  min.segment.length = 0,
	  max.overlaps = Inf
	) + 
  scale_y_log10() +
  labs(
    x = "Increasing Average Expression",
    y = "Mean Expression of Top 3 Samples",
    title = "Ranked gene expression (MDMX-Tg v. WT)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_blank()
  )

ggplot(wt_plot, aes(x = rank, y = top3_mean)) +
  geom_line(color = "gray40") +
  geom_point(
    data = wt_plot[mgi_symbol %in% wt_de],
    aes(x = rank, y = top3_mean),
    color = "blue",
    size = 2
  ) +
  geom_point(data = esam_row, aes(x = rank, y = top3_mean),
             color = "red", size = 3) +
  geom_text_repel(
	  data = esam_row,
	  aes(x = rank, y = top3_mean, label = "Esam"),
	  color = "red",
	  size = 4,
	  nudge_y = 0.5,         # push upward
	  nudge_x = 5,           # push right slightly
	  box.padding = 0.5,
	  min.segment.length = 0,
	  max.overlaps = Inf
	) + 
  geom_point(data = tet2_row, aes(x = rank, y = top3_mean),
             color = "red", size = 3) +
  geom_text_repel(
	  data = tet2_row,
	  aes(x = rank, y = top3_mean, label = "TET2"),
	  color = "red",
	  size = 4,
	  nudge_y = 0.5,         # push upward
	  nudge_x = 5,           # push right slightly
	  box.padding = 0.5,
	  min.segment.length = 0,
	  max.overlaps = Inf
	) + 
  scale_y_log10() +
  labs(
    x = "Increasing Average Expression",
    y = "Mean Expression of Top 3 Samples",
    title = "Ranked gene expression (MDMX-Tg v. WT)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_blank()
  )

ggplot(wt_plot, aes(x = rank, y = top3_mean)) +
  geom_line(color = "gray40") +
  geom_point(
    data = wt_plot[mgi_symbol %in% wt_de],
    aes(x = rank, y = top3_mean),
    color = "blue",
    size = 2
  ) +
  geom_text_repel(
    data = wt_plot[mgi_symbol %in% wt_de],
    aes(x = rank, y = top3_mean, label = mgi_symbol),
    color = "blue",
    size = 3,
    max.overlaps = Inf,
    min.segment.length = 0
  ) +
  geom_point(data = esam_row, aes(x = rank, y = top3_mean),
             color = "red", size = 3) +
  geom_text_repel(
	  data = esam_row,
	  aes(x = rank, y = top3_mean, label = "Esam"),
	  color = "red",
	  size = 4,
	  nudge_y = 0.5,         # push upward
	  nudge_x = 5,           # push right slightly
	  box.padding = 0.5,
	  min.segment.length = 0,
	  max.overlaps = Inf
	) + 
  geom_point(data = tet2_row, aes(x = rank, y = top3_mean),
             color = "red", size = 3) +
  geom_text_repel(
	  data = tet2_row,
	  aes(x = rank, y = top3_mean, label = "TET2"),
	  color = "red",
	  size = 4,
	  nudge_y = 0.5,         # push upward
	  nudge_x = 5,           # push right slightly
	  box.padding = 0.5,
	  min.segment.length = 0,
	  max.overlaps = Inf
	) + 
  scale_y_log10() +
  labs(
    x = "Increasing Average Expression",
    y = "Mean Expression of Top 3 Samples",
    title = "Ranked gene expression (MDMX-Tg v. WT)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_blank()
  )

dev.off()

############

ure <- fread("inputs/ure_v_ure_mdmx.csv", skip = 2, header = TRUE)
ure[mgi_symbol == "Esam"]
count_cols <- c(
  "M356URE.Mdmx",
  "M358URE",
  "M372URE",
  "M373URE",
  "M381URE.Mdmx",
  "M384URE.Mdmx"
)

ure[, top3_mean := apply(.SD, 1, function(x) mean(sort(x, decreasing = TRUE)[1:3])),
    .SDcols = count_cols]

ure_de <- fread("outputs/ure_v_ure_mdmx.csv") %>% 
	.[neglog10_padj > -log10(0.05)] %>% 
	.$mgi_symbol

# Rank genes by top3_mean
ure_plot <- ure %>%
  arrange(top3_mean) %>%
  mutate(rank = row_number())

# Identify Esam row (if exists)
esam_row <- ure_plot %>% filter(mgi_symbol == "Esam")
tet2_row <- ure_plot %>% filter(mgi_symbol == "Tet2")

pdf("outputs/out7.pdf")

ggplot(ure_plot, aes(x = rank, y = top3_mean)) +
  geom_line(color = "gray40") +
  scale_y_log10() +
  labs(
    x = "Increasing Average Expression",
    y = "Mean Expression of Top 3 Samples",
    title = "Ranked gene expression (URE v. MDMX-Tg/URE)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_blank()
  )

ggplot(ure_plot, aes(x = rank, y = top3_mean)) +
  geom_line(color = "gray40") +
  geom_point(data = esam_row, aes(x = rank, y = top3_mean),
             color = "red", size = 3) +
  geom_text_repel(
	  data = esam_row,
	  aes(x = rank, y = top3_mean, label = "Esam"),
	  color = "red",
	  size = 4,
	  nudge_y = 0.5,         # push upward
	  nudge_x = 5,           # push right slightly
	  box.padding = 0.5,
	  min.segment.length = 0,
	  max.overlaps = Inf
	) + 
  geom_point(data = tet2_row, aes(x = rank, y = top3_mean),
             color = "red", size = 3) +
  geom_text_repel(
	  data = tet2_row,
	  aes(x = rank, y = top3_mean, label = "TET2"),
	  color = "red",
	  size = 4,
	  nudge_y = 0.5,         # push upward
	  nudge_x = 5,           # push right slightly
	  box.padding = 0.5,
	  min.segment.length = 0,
	  max.overlaps = Inf
	) + 
  scale_y_log10() +
  labs(
    x = "Increasing Average Expression",
    y = "Mean Expression of Top 3 Samples",
    title = "Ranked gene expression (URE v. MDMX-Tg/URE)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_blank()
  )

ggplot(ure_plot, aes(x = rank, y = top3_mean)) +
  geom_line(color = "gray40") +
  geom_point(
    data = ure_plot[mgi_symbol %in% ure_de],
    aes(x = rank, y = top3_mean),
    color = "blue",
    size = 2
  ) +
  geom_point(data = esam_row, aes(x = rank, y = top3_mean),
             color = "red", size = 3) +
  geom_text_repel(
	  data = esam_row,
	  aes(x = rank, y = top3_mean, label = "Esam"),
	  color = "red",
	  size = 4,
	  nudge_y = 0.5,         # push upward
	  nudge_x = 5,           # push right slightly
	  box.padding = 0.5,
	  min.segment.length = 0,
	  max.overlaps = Inf
	) + 
  geom_point(data = tet2_row, aes(x = rank, y = top3_mean),
             color = "red", size = 3) +
  geom_text_repel(
	  data = tet2_row,
	  aes(x = rank, y = top3_mean, label = "TET2"),
	  color = "red",
	  size = 4,
	  nudge_y = 0.5,         # push upward
	  nudge_x = 5,           # push right slightly
	  box.padding = 0.5,
	  min.segment.length = 0,
	  max.overlaps = Inf
	) + 
  scale_y_log10() +
  labs(
    x = "Increasing Average Expression",
    y = "Mean Expression of Top 3 Samples",
    title = "Ranked gene expression (URE v. MDMX-Tg/URE)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_blank()
  )

ggplot(ure_plot, aes(x = rank, y = top3_mean)) +
  geom_line(color = "gray40") +
  geom_point(
    data = ure_plot[mgi_symbol %in% ure_de],
    aes(x = rank, y = top3_mean),
    color = "blue",
    size = 2
  ) +
  geom_text_repel(
    data = ure_plot[mgi_symbol %in% ure_de],
    aes(x = rank, y = top3_mean, label = mgi_symbol),
    color = "blue",
    size = 3,
    max.overlaps = Inf,
    min.segment.length = 0
  ) +
  geom_point(data = esam_row, aes(x = rank, y = top3_mean),
             color = "red", size = 3) +
  geom_text_repel(
	  data = esam_row,
	  aes(x = rank, y = top3_mean, label = "Esam"),
	  color = "red",
	  size = 4,
	  nudge_y = 0.8,         # push upward
	  nudge_x = 5,           # push right slightly
	  box.padding = 0.5,
	  min.segment.length = 0,
	  max.overlaps = Inf
	) + 
  geom_point(data = tet2_row, aes(x = rank, y = top3_mean),
             color = "red", size = 3) +
  geom_text_repel(
	  data = tet2_row,
	  aes(x = rank, y = top3_mean, label = "TET2"),
	  color = "red",
	  size = 4,
	  nudge_y = 0.8,         # push upward
	  nudge_x = 5,           # push right slightly
	  box.padding = 0.5,
	  min.segment.length = 0,
	  max.overlaps = Inf
	) + 
  scale_y_log10() +
  labs(
    x = "Increasing Average Expression",
    y = "Mean Expression of Top 3 Samples",
    title = "Ranked gene expression (URE v. MDMX-Tg/URE)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_blank()
  )

dev.off()
