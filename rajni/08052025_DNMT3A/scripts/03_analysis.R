library(Seurat)
library(data.table)
library(dplyr)
library(Matrix)
library(ggplot2)

# Now I have a RDS object that is mapped to cell types 
# Only really care about the HSCs and how the expression of these 5 genes looks in each of the samples 

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

# Plot
pdf("outputs/HSC_umap.pdf")
	DimPlot(hsc, group.by = "batch", reduction = "umap")
	FeaturePlot(hsc, features = c("IRF8","RUNX1"))
dev.off()

# One sample
sample_id <- "GSM5171988"

# 1) subset to just that sample (assumes metadata column 'batch')
s <- subset(hsc, subset = batch == sample_id)

# 2) make sure we have a normalized layer for plotting
if (!"data" %in% Layers(s[["RNA"]])) s <- NormalizeData(s)

# 3) side-by-side UMAP feature plots for Irf8 and Runx1
p <- FeaturePlot(
  s,
  features = c("IRF8","RUNX1"),
  reduction = "umap",
  # layer = "data",        # Seurat v5: plot normalized data layer
  order = TRUE,
  ncol = 2
)

pdf("outputs/runx_irf8.pdf")
  print(p)
dev.off()

####################

z_cutoff <- 1

genes <- c("IRF8","RUNX1")
DefaultAssay(hsc) <- "RNA"

# Ensure normalized layer exists for plotting
if (!"data" %in% Layers(hsc[["RNA"]])) hsc <- NormalizeData(hsc)

# Create z-scores for just these genes (writes to 'scale.data' layer)
hsc <- ScaleData(hsc, features = genes, verbose = FALSE)

# Pull the z-scores for the two genes
Z <- t(LayerData(hsc[["RNA"]], layer = "scale.data")[genes, , drop = FALSE])
colnames(Z) <- genes

# Flag cells with BOTH genes high (z >= 1)
hsc$IRF8_RUNX1_both_high <- (Z[, "IRF8"] >= z_cutoff) & (Z[, "RUNX1"] >= z_cutoff)

# (Optional) 4-way state: both-high / Irf8-only / Runx1-only / neither
hsc$IRF8_RUNX1_state <- c("neither","IRF8-only","RUNX1-only","both-high")[
  1 + (Z[, "IRF8"] >= z_cutoff) + 2*(Z[, "RUNX1"] >= z_cutoff)
]



sample_id <- "GSM5171988"
s <- subset(hsc, subset = batch == sample_id)

pdf("out.pdf")

  # gene-wise expression (normalized)
  FeaturePlot(s, features = genes, reduction = "umap", ncol = 2, order = TRUE)

  # highlight double-high cells on the UMAP
  cells_hi <- colnames(s)[s$IRF8_RUNX1_both_high]
  DimPlot(s, reduction = "umap",
          cells.highlight = cells_hi,
          cols = c("grey80","red"), sizes = c(0.2, 0.8)) +
    ggtitle(paste0(sample_id, ": IRF8 & RUNX1 both-high"))

dev.off()


############
# CONTROL
ctrl <- c("GSM5171988", "GSM5171989", "GSM5257698")

s <- subset(hsc, subset = batch %in% ctrl)

s$IRF8_RUNX1_state <- droplevels(factor(
  s$IRF8_RUNX1_state,
  levels = c("neither","IRF8-only","RUNX1-only","both-high")
))

pal <- c("neither"    = "grey80",
         "IRF8-only"  = "#1f77b4",
         "RUNX1-only" = "#ff7f0e",
         "both-high"  = "#d62728")

Z <- t(LayerData(s[["RNA"]], layer = "scale.data")[genes, , drop = FALSE])
df <- data.frame(IRF8 = Z[, "IRF8"], RUNX1 = Z[, "RUNX1"],
                 state = s$IRF8_RUNX1_state)

pdf("outputs/ctrl.pdf", width = 10)

  # gene-wise expression (normalized)
  FeaturePlot(s, features = genes, reduction = "umap", ncol = 2, order = TRUE)

  DimPlot(
    s, reduction = "umap",
    group.by = "IRF8_RUNX1_state",
    cols = pal,
    order = c("both-low","IRF8-only","RUNX1-only","both-high"),  # plot both-high on top
    pt.size = 0.5
  ) + ggtitle("CTRL — IRF8/RUNX1 state")

  ggplot(df, aes(x = IRF8, y = RUNX1)) + 
    geom_point() + 
    geom_vline(xintercept = 1, linetype = 2) +
    geom_hline(yintercept = 1, linetype = 2) + 
    labs(x = "IRF8 z-score", y = "RUNX1 z-score")

dev.off()

############
# EXPERIMENTAL
exp <- c("GSM5171990", "GSM5171991", "GSM5257699", "GSM5257700")

s <- subset(hsc, subset = batch %in% exp)

s$IRF8_RUNX1_state <- droplevels(factor(
  s$IRF8_RUNX1_state,
  levels = c("neither","IRF8-only","RUNX1-only","both-high")
))

pal <- c("neither"    = "grey80",
         "IRF8-only"  = "#1f77b4",
         "RUNX1-only" = "#ff7f0e",
         "both-high"  = "#d62728")

pdf("outputs/exp.pdf")

  # gene-wise expression (normalized)
  FeaturePlot(s, features = genes, reduction = "umap", ncol = 2, order = TRUE)

  DimPlot(
    s, reduction = "umap",
    group.by = "IRF8_RUNX1_state",
    cols = pal,
    order = c("both-low","IRF8-only","RUNX1-only","both-high"),  # plot both-high on top
    pt.size = 0.5
  ) + ggtitle("DNMT3a KO — IRF8/RUNX1 state")

dev.off()

#################

obj <- hsc
samples <- ctrl
g1 <- "KIT"
g2 <- "LY6A"
z_cut <- 1

plot_pair <- function(obj, samples, g1, g2, z_cut=1) {
  s <- subset(obj, subset = batch %in% samples)

  feats <- unique(c(g1,g2))

  s <- ScaleData(s, features = feats, verbose=FALSE)

  # if (!"scale.data" %in% Layers(s[["RNA"]])) s <- ScaleData(s, features = feats, verbose=FALSE)

  Z <- t(LayerData(s[["RNA"]], layer="scale.data")[feats, , drop=FALSE])
  state <- c("neither","g1-only","g2-only","both-high")[
    1 + (Z[, g1] >= z_cut) + 2*(Z[, g2] >= z_cut)
  ]

  state <- factor(state, levels=c("neither","g1-only","g2-only","both-high"))

  p_umap <- DimPlot(s, reduction="umap", group.by=NULL, cols="grey80") +
    ggtitle(paste0(samples, " — ", g1, " & ", g2))

  p_state <- DimPlot(s, reduction="umap", group.by=NULL,
                     cells.highlight = colnames(s)[state=="both-high"],
                     cols=c("grey80","red"))

  list(
    umap=p_umap, 
    both_high=p_state, 
    # scatter = scatter,
    frac_both_high=mean(state=="both-high")
  )
}

plot_pair2 <- function(obj, samples, g1, g2, z_cut = 1) {
  s <- subset(obj, subset = batch %in% samples)
  DefaultAssay(s) <- "RNA"

  feats <- unique(c(g1, g2))
  miss  <- setdiff(feats, rownames(s))
  if (length(miss)) stop("Missing genes: ", paste(miss, collapse = ", "))

  # use normalized data -> gene-wise z across cells (robust in Seurat v5)
  if (!"data" %in% Layers(s[["RNA"]])) s <- NormalizeData(s, verbose = FALSE)
  D <- LayerData(s[["RNA"]], layer = "data")[feats, , drop = FALSE]
  Z <- t(scale(t(as.matrix(D))))                 # cells x genes (z-scores)
  rownames(Z) <- feats
  Z <- t(Z)

  # 4-state call, labeled by gene names
  state <- c("neither", "g1-only", "g2-only", "both-high")[
    1 + (Z[, g1] >= z_cut) + 2 * (Z[, g2] >= z_cut)
  ]
  s$PAIR_state <- factor(
    state,
    # order legend like your example: neither, both-high, g2-only, g1-only
    levels = c("neither", "both-high", "g2-only", "g1-only")
  )

  # palette like the attached figure
  pal <- c(
    "neither" = "grey80",
    "both-high" = "firebrick3",          # red
    "g1-only" = "lightsteelblue2",  # orange
    "g2-only" = "darkseagreen2"   # blue
  )

  title_txt <- paste0(paste(samples, collapse = ", "),
                      " — ", toupper(g1), "/", toupper(g2), " state")

  # >>> THIS is the revised p_umap <<<
  p_umap <- DimPlot(
    s, reduction = "umap",
    group.by = "PAIR_state",
    cols = pal, 
    pt.size = 0.5,
    # draw order so both-high is plotted last (on top)
    # order = c("neither", paste0(g2, "-only"), paste0(g1, "-only"), "both-high")
  ) + ggtitle(title_txt)

  # keep your both-high highlight as a second panel (optional)
  p_state <- DimPlot(
    s, reduction = "umap", group.by = NULL,
    cells.highlight = colnames(s)[s$PAIR_state == "both-high"],
    cols = c("grey80", pal["both-high"]), pt.size = 0.5
  ) + ggtitle(paste0(title_txt, " — both-high"))

  list(
    umap = p_umap,
    both_high = p_state,
    frac_both_high = mean(s$PAIR_state == "both-high")
  )
}


ctrl <- c("GSM5171988", "GSM5171989", "GSM5257698")
exp <- c("GSM5171990", "GSM5171991", "GSM5257699", "GSM5257700")

res_pos <- plot_pair(hsc, ctrl, "KIT", "LY6A", z_cut=0.5)
res_pos$frac_both_high  # should be noticeably > your IRF8/RUNX1 fraction
res_pos <- plot_pair(hsc, ctrl, "GAPDH", "ACTB", z_cut=0.5)
res_pos$frac_both_high 
res_pos <- plot_pair(hsc, ctrl, "MECOM", "HLF", z_cut=0.5)
res_pos$frac_both_high 
res_pos <- plot_pair(hsc, ctrl, "PROCR", "HLF", z_cut=0.5)
res_pos$frac_both_high 
res_pos <- plot_pair(hsc, ctrl, "SLAMF1", "MECOM", z_cut=0.5)
res_pos$frac_both_high 
res_pos <- plot_pair(hsc, ctrl, "IRF8", "RUNX1", z_cut=0.5)
res_pos$frac_both_high



res_pos <- plot_pair2(hsc, ctrl, "IRF8", "RUNX1", z_cut=0.5)

pdf("pairs2.pdf")
  print(res_pos$umap)
  print(res_pos$both_high)
dev.off()

res_pos$both_high       # shows clear co-high cluster(s)















