library(Seurat)
library(data.table)
library(dplyr)
library(Matrix)
library(ggplot2)

# Pre-process
obj <- readRDS("outputs/combined_mapped.rds")

obj@meta.data$predicted_CellType_No_Abbreviations %>% unique()

hsc <- subset(
  obj, 
  subset = predicted_CellType_No_Abbreviations == "Hematopoietic Stem Cell"
)

# Normalize + UMAP
hsc <- NormalizeData(hsc)                                   # creates "data" layer in v5
hsc <- FindVariableFeatures(hsc, nfeatures = 3000)

# PCA → neighbors → UMAP → clusters
hsc <- RunPCA(hsc, features = VariableFeatures(hsc))
hsc <- FindNeighbors(hsc, dims = 1:20)
hsc <- FindClusters(hsc, resolution = 0.2)
hsc <- RunUMAP(hsc, dims = 1:20)

pdf("01.hsc-umap.pdf")
	DimPlot(hsc, reduction = "umap", 
		label = FALSE,
		pt.size = 0.4)
dev.off()

obj <- hsc
g1 <- "RUNX1"
g2 <- "IRF8"
title <- "Test"
thresh = 0.25

coexp_umap <- function(obj, g1, g2, thresh = 0.25, title = NULL) {
  DefaultAssay(obj) <- "RNA"
  if (!"data" %in% Layers(obj[["RNA"]])) obj <- NormalizeData(obj, verbose = FALSE)

  genes <- c(g1, g2)
  miss  <- setdiff(genes, rownames(obj))
  if (length(miss)) stop("Missing genes: ", paste(miss, collapse = ", "))

  E <- LayerData(obj[["RNA"]], layer = "data")[genes, , drop = FALSE]
  g1hi <- as.numeric(E[g1, ]) > thresh
  g2hi <- as.numeric(E[g2, ]) > thresh

  state <- c("neither", "g1-only", "g2-only", "both-high")[1 + g1hi + 2*g2hi]
  state <- factor(state, levels = c("neither","both-high", "g2-only", "g1-only"))
  obj$coexp_state <- state

  pal <- c(
  	"neither"="grey80", 
  	"both-high"="#b2182b",
  	"g2-only" = "#fdb863",
  	"g1-only" = "#1f78b4"
	)

  if (is.null(title)) title <- paste0("HSC — ", g1, "/", g2, " state (>", thresh, ")")

  p <- DimPlot(
    obj, reduction = "umap",
    group.by = "coexp_state",
    cols = pal, pt.size = 0.5,
    order = levels(state)  # draw both-high last
  ) + ggtitle(title)

  list(plot = p,
  	frac_both = mean(state == "both-high"),
  	n_both = sum(state == "both-high"),
  	frac_g1_high = mean(state == "g1-only"),
  	frac_g2_high = mean(state == "g2-only")
	)
}

# Positive Control Pairs
pos_candidates <- list(
  c("Gapdh","Actb"),
  c("Rpl13a","Rplp0"),
  c("Ppia","B2m")
) 

	# UMAPs
pdf("02.pos_controls.pdf")

for(p in pos_candidates) {
	res <- coexp_umap(obj, toupper(p[1]), toupper(p[2]))
	print(res$plot)
}

dev.off()
	
	# Bar plots
pos_dt <- lapply(seq_along(pos_candidates), function(p) {
	# print(p)
	g1 <- pos_candidates[[p]][1] %>% toupper()
	g2 <- pos_candidates[[p]][2] %>% toupper()

	res <- coexp_umap(obj, g1, g2)

	data.table(
		g1 = g1,
		g2 = g2,
		frac_both = res$frac_both,
		frac_g1_high = res$frac_g1_high,
		frac_g2_high = res$frac_g2_high,
		group = "pos ctrl"
	)
}) %>% rbindlist()

pos_dt$name <- paste0(pos_dt$g1, " & ", pos_dt$g2)

# Negative Control Pairs 
neg_candidates <- list(
	c("LYZ2", "KLF1"),
	c("MKI67", "TOP2A"),
	c("Procr","Flt3")
)

	# UMAPs
pdf("03.neg_controls.pdf")

for(p in neg_candidates) {
	res <- coexp_umap(obj, toupper(p[1]), toupper(p[2]))
	print(res$plot)
}

dev.off()

for(p in neg_candidates) {
	res <- coexp_umap(obj, toupper(p[1]), toupper(p[2]))
	print(res$frac_both)
}

neg_dt <- lapply(seq_along(pos_candidates), function(p) {
	# print(p)
	g1 <- neg_candidates[[p]][1] %>% toupper()
	g2 <- neg_candidates[[p]][2] %>% toupper()

	res <- coexp_umap(obj, g1, g2)

	data.table(
		g1 = g1,
		g2 = g2,
		frac_both = res$frac_both,
		frac_g1_high = res$frac_g1_high,
		frac_g2_high = res$frac_g2_high,
		group = "neg ctrl"
	)
}) %>% rbindlist()

neg_dt$name <- paste0(neg_dt$g1, " & ", neg_dt$g2)

dt <- rbind(pos_dt, neg_dt)

pdf("04.all_ctrls.pdf")
	ggplot(dt, aes(x = name, y = frac_both)) + 
		geom_col()
dev.off()

dt[, group := factor(group, levels = c("pos ctrl","neg ctrl"))]
setorder(dt, group, -frac_both)                 # or: setorder(dt, group, name)
dt[, name := factor(name, levels = dt$name)]     # lock x-axis order

n_pos <- sum(dt$group == "pos ctrl")

pdf("05.all_ctrls_grouped.pdf", width = 7, height = 3.5)
ggplot(dt, aes(x = name, y = frac_both, fill = group)) +
  geom_col() +
  geom_vline(xintercept = n_pos + 0.5, linetype = 2) +  # visual split
  scale_fill_manual(values = c("pos ctrl" = "#1f78b4", "neg ctrl" = "#b2df8a")) +
  labs(x = NULL, y = "fraction both-high") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top")
dev.off()

# IRF8 and RUNX1

res_irf8_runx1 <- coexp_umap(hsc, "IRF8", "RUNX1", thresh = 0.25, title = "IRF8/RUNX1 state")

pdf("06.RI_state.pdf")
	print(res_irf8_runx1$plot)
dev.off()

my_pair <- data.table(
	g1 = "IRF8",
	g2 = "RUNX1",
	frac_both = res_irf8_runx1$frac_both,
	frac_g1_high = res_irf8_runx1$frac_g1_high,
	frac_g2_high = res_irf8_runx1$frac_g2_high,
	group = "exp",
	name = "IRF8 & RUNX1"
)

dt2 <- rbind(dt, my_pair)

dt2[, group := factor(group, levels = c("pos ctrl","neg ctrl", "exp"))]
setorder(dt2, group, -frac_both)                 # or: setorder(dt, group, name)
dt2[, name := factor(name, levels = dt2$name)]     # lock x-axis order

n_pos <- sum(dt2$group == "pos ctrl")

pdf("07.all.pdf", width = 7, height = 3.5)
ggplot(dt2, aes(x = name, y = frac_both, fill = group)) +
  geom_col() +
  geom_vline(xintercept = n_pos + 0.5, linetype = 2) +  # visual split
  geom_vline(xintercept = n_pos + n_pos + 0.5, linetype = 2) +  # visual split
  scale_fill_manual(values = c("pos ctrl" = "#1f78b4", "neg ctrl" = "#b2df8a", "exp" = "firebrick3")) +
  labs(x = NULL, y = "fraction both-high") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top")
dev.off()

# WT v. DNTM3A KO
wt <- c("GSM5171988", "GSM5171989", "GSM5257698")
ko <- c("GSM5171990", "GSM5171991", "GSM5257699", "GSM5257700")

	# UMAP
res1 <- coexp_umap(subset(hsc, subset = batch %in% wt), "IRF8", "RUNX1", thresh = 0.25, title = "IRF8/RUNX1 state")
res2 <- coexp_umap(subset(hsc, subset = batch %in% ko), "IRF8", "RUNX1", thresh = 0.25, title = "IRF8/RUNX1 state")

pdf("08.wt_ko_umap.pdf")
	print(res1$plot)
	print(res2$plot)
dev.off()

# WILD TYPE 
pos_candidates <- list(
  c("Gapdh","Actb"),
  c("Rpl13a","Rplp0"),
  c("Ppia","B2m")
) 

pos_dt <- lapply(seq_along(pos_candidates), function(p) {
	g1 <- pos_candidates[[p]][1] %>% toupper()
	g2 <- pos_candidates[[p]][2] %>% toupper()

	res <- coexp_umap(subset(hsc, subset = batch %in% wt), g1, g2)

	data.table(
		batch = "WT",
		g1 = g1,
		g2 = g2,
		frac_both = res$frac_both,
		frac_g1_high = res$frac_g1_high,
		frac_g2_high = res$frac_g2_high,
		group = "pos ctrl"
	)
}) %>% rbindlist()

pos_dt$name <- paste0(pos_dt$g1, " & ", pos_dt$g2)

neg_candidates <- list(
	c("LYZ2", "KLF1"),
	c("MKI67", "TOP2A"),
	c("Procr","Flt3")
)

neg_dt <- lapply(seq_along(pos_candidates), function(p) {
	g1 <- neg_candidates[[p]][1] %>% toupper()
	g2 <- neg_candidates[[p]][2] %>% toupper()

	res <- coexp_umap(subset(hsc, subset = batch %in% wt), g1, g2)

	data.table(
		batch = "WT",
		g1 = g1,
		g2 = g2,
		frac_both = res$frac_both,
		frac_g1_high = res$frac_g1_high,
		frac_g2_high = res$frac_g2_high,
		group = "neg ctrl"
	)
}) %>% rbindlist()

neg_dt$name <- paste0(neg_dt$g1, " & ", neg_dt$g2)

dt <- rbind(pos_dt, neg_dt)

res_irf8_runx1 <- coexp_umap(subset(hsc, subset = batch %in% wt), "IRF8", "RUNX1", thresh = 0.25, title = "IRF8/RUNX1 state")

my_pair <- data.table(
	batch = "WT",
	g1 = "IRF8",
	g2 = "RUNX1",
	frac_both = res_irf8_runx1$frac_both,
	frac_g1_high = res_irf8_runx1$frac_g1_high,
	frac_g2_high = res_irf8_runx1$frac_g2_high,
	group = "exp",
	name = "IRF8 & RUNX1"
)

dt2 <- rbind(dt, my_pair)

dt2[, group := factor(group, levels = c("pos ctrl","neg ctrl", "exp"))]
setorder(dt2, group, -frac_both)                 # or: setorder(dt, group, name)
dt2[, name := factor(name, levels = dt2$name)]     # lock x-axis order

wt_dt <- dt2

n_pos <- sum(dt2$group == "pos ctrl")

pdf("09.all_WT.pdf", width = 7, height = 3.5)
ggplot(dt2, aes(x = name, y = frac_both, fill = group)) +
  geom_col() +
  geom_text(aes(label = round(frac_both, digits = 3)), vjust = -0.3, size = 3) + 
  geom_vline(xintercept = n_pos + 0.5, linetype = 2) +  # visual split
  geom_vline(xintercept = n_pos + n_pos + 0.5, linetype = 2) +  # visual split
  scale_fill_manual(values = c("pos ctrl" = "#1f78b4", "neg ctrl" = "#b2df8a", "exp" = "firebrick3")) +
  labs(x = NULL, y = "fraction both-high") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top") + 
  labs(x = "Gene Pair", y = "Fraction both-high", title = "Wild Type")
dev.off()

# KNOCKOUT 
pos_candidates <- list(
  c("Gapdh","Actb"),
  c("Rpl13a","Rplp0"),
  c("Ppia","B2m")
) 

pos_dt <- lapply(seq_along(pos_candidates), function(p) {
	g1 <- pos_candidates[[p]][1] %>% toupper()
	g2 <- pos_candidates[[p]][2] %>% toupper()

	res <- coexp_umap(subset(hsc, subset = batch %in% ko), g1, g2)

	data.table(
		batch = "KO",
		g1 = g1,
		g2 = g2,
		frac_both = res$frac_both,
		frac_g1_high = res$frac_g1_high,
		frac_g2_high = res$frac_g2_high,
		group = "pos ctrl"
	)
}) %>% rbindlist()

pos_dt$name <- paste0(pos_dt$g1, " & ", pos_dt$g2)

neg_candidates <- list(
	c("LYZ2", "KLF1"),
	c("MKI67", "TOP2A"),
	c("Procr","Flt3")
)

neg_dt <- lapply(seq_along(pos_candidates), function(p) {
	g1 <- neg_candidates[[p]][1] %>% toupper()
	g2 <- neg_candidates[[p]][2] %>% toupper()

	res <- coexp_umap(subset(hsc, subset = batch %in% ko), g1, g2)

	data.table(
		batch = "KO",
		g1 = g1,
		g2 = g2,
		frac_both = res$frac_both,
		frac_g1_high = res$frac_g1_high,
		frac_g2_high = res$frac_g2_high,
		group = "neg ctrl"
	)
}) %>% rbindlist()

neg_dt$name <- paste0(neg_dt$g1, " & ", neg_dt$g2)

dt <- rbind(pos_dt, neg_dt)

res_irf8_runx1 <- coexp_umap(subset(hsc, subset = batch %in% ko), "IRF8", "RUNX1", thresh = 0.25, title = "IRF8/RUNX1 state")

my_pair <- data.table(
	batch = "KO",
	g1 = "IRF8",
	g2 = "RUNX1",
	frac_both = res_irf8_runx1$frac_both,
	frac_g1_high = res_irf8_runx1$frac_g1_high,
	frac_g2_high = res_irf8_runx1$frac_g2_high,
	group = "exp",
	name = "IRF8 & RUNX1"
)

dt2 <- rbind(dt, my_pair)

dt2[, group := factor(group, levels = c("pos ctrl","neg ctrl", "exp"))]
setorder(dt2, group, -frac_both)                 # or: setorder(dt, group, name)
dt2[, name := factor(name, levels = dt2$name)]     # lock x-axis order

ko_dt <- dt2

n_pos <- sum(dt2$group == "pos ctrl")

pdf("10.all_KO.pdf", width = 7, height = 3.5)
	ggplot(dt2, aes(x = name, y = frac_both, fill = group)) +
	  geom_col() +
	  geom_text(aes(label = round(frac_both, digits = 3)), vjust = -0.3, size = 3) + 
	  geom_vline(xintercept = n_pos + 0.5, linetype = 2) +  # visual split
	  geom_vline(xintercept = n_pos + n_pos + 0.5, linetype = 2) +  # visual split
	  scale_fill_manual(values = c("pos ctrl" = "#1f78b4", "neg ctrl" = "#b2df8a", "exp" = "firebrick3")) +
	  labs(x = NULL, y = "fraction both-high") +
	  theme_classic() +
	  theme(axis.text.x = element_text(angle = 45, hjust = 1),
	        legend.position = "top") + 
	  labs(x = "Gene Pair", y = "Fraction both-high", title = "DNMT3A Knockout")
dev.off()

# Fold Change 

dt_both <- merge(
	wt_dt[, c("batch", "group", "name", "frac_both"), with = FALSE],
	ko_dt[, c("batch", "group", "name", "frac_both"), with = FALSE],
	by = "name"
)

dt_both$logFC <- log2(dt_both$frac_both.x / dt_both$frac_both.y)

pdf("11.FC.pdf")
ggplot(dt_both, aes(x = name, y = logFC, fill = group.x)) + 
	geom_col() + 
	theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "top") + 
  geom_vline(xintercept = n_pos + 0.5, linetype = 2) +  # visual split
  geom_vline(xintercept = n_pos + n_pos + 0.5, linetype = 2) +  # visual split
  scale_fill_manual(values = c("pos ctrl" = "#1f78b4", "neg ctrl" = "#b2df8a", "exp" = "firebrick3")) + 
  labs(x = "Gene Pair", y = "Log Fold Change (Wildtype/Knockout)") + 
  ylim(-2, 2)
dev.off()


