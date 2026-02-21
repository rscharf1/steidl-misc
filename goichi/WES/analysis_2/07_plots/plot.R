library(tidyverse)
library(data.table)
library(ggplot2)

dt <- fread("vaf.tsv")

dt <- dt[, !c(5,8)]

colnames(dt) <- c(
    "chr", "pos",
    "AD_ctrl", "AD_treat",
    "DP_ctrl", "DP_treat"
  )

calc_vaf <- function(ad, dp) {
  counts <- as.numeric(strsplit(ad, ",")[[1]])
  counts[2] / as.numeric(dp)
}

dt[, `:=`(
  vaf_ctrl = mapply(calc_vaf, AD_ctrl, DP_ctrl),
  vaf_treat  = mapply(calc_vaf, AD_treat,  DP_treat)
)]

dt <- dt[!(is.na(vaf_ctrl) | is.na(vaf_treat))]

r_cor <- cor(dt$vaf_ctrl, dt$vaf_treat)

pdf("out.pdf")

  ggplot(dt, aes(vaf_ctrl, vaf_treat)) +
    geom_point(alpha = 0.25, size = 0.6) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    annotate(
      "text",
      x = 0.05, y = 0.95,
      label = paste0("Pearson r = ", round(r_cor, 3)),
      hjust = 0
    ) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(
      title = "VAF concordance: Control vs Treat",
      x = "Control VAF",
      y = "Treat VAF"
    ) +
    theme_classic()

dev.off()






