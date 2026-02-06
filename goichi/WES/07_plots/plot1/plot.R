library(tidyverse)
library(data.table)
library(ggplot2)

dt <- fread("plot1/vaf.tsv")

dt <- dt[, !c(6,10)]

colnames(dt) <- c(
    "chr", "pos",
    "AD_ctrl", "AD_pyh", "AD_pyl",
    "DP_ctrl", "DP_pyh", "DP_pyl"
  )

calc_vaf <- function(ad, dp) {
  counts <- as.numeric(strsplit(ad, ",")[[1]])
  counts[2] / as.numeric(dp)
}

dt[, `:=`(
  vaf_ctrl = mapply(calc_vaf, AD_ctrl, DP_ctrl),
  vaf_pyh  = mapply(calc_vaf, AD_pyh,  DP_pyh),
  vaf_pyl  = mapply(calc_vaf, AD_pyl,  DP_pyl)
)]

dt <- dt[!(is.na(vaf_ctrl) | is.na(vaf_pyh) | is.na(vaf_pyl))]

r_pyh <- cor(dt$vaf_ctrl, dt$vaf_pyh)
r_pyl <- cor(dt$vaf_ctrl, dt$vaf_pyl)

pdf("plot1/out.pdf")

  ggplot(dt, aes(vaf_ctrl, vaf_pyh)) +
    geom_point(alpha = 0.25, size = 0.6) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    annotate(
      "text",
      x = 0.05, y = 0.95,
      label = paste0("Pearson r = ", round(r_pyh, 3)),
      hjust = 0
    ) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(
      title = "VAF concordance: Control vs PY_high",
      x = "Control VAF",
      y = "PY_high VAF"
    ) +
    theme_classic()

  ggplot(dt, aes(vaf_ctrl, vaf_pyl)) +
    geom_point(alpha = 0.25, size = 0.6) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    annotate(
      "text",
      x = 0.05, y = 0.95,
      label = paste0("Pearson r = ", round(r_pyl, 3)),
      hjust = 0
    ) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(
      title = "VAF concordance: Control vs PY_low",
      x = "Control VAF",
      y = "PY_low VAF"
    ) +
    theme_classic()

  ggplot(dt, aes(x = vaf_pyh - vaf_ctrl)) +
    geom_histogram(
      bins = 50,
      fill = "grey70",
      color = "black"
    ) +
    labs(
      title = "VAF difference (PY_high - Control)",
      x = expression(Delta*VAF),
      y = "Variant count"
    ) +
    theme_classic()

  ggplot(dt, aes(x = vaf_pyl - vaf_ctrl)) +
    geom_histogram(
      bins = 50,
      fill = "grey70",
      color = "black"
    ) +
    labs(
      title = "VAF difference (PY_low - Control)",
      x = expression(Delta*VAF),
      y = "Variant count"
    ) +
    theme_classic()

dev.off()






