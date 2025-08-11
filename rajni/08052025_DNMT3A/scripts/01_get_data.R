suppressPackageStartupMessages({
  library(forcats)
  library(ggplot2)
  library(GEOquery)
  library(preprocessCore) # quantile normalization
  library(tidyverse)
  library(viridis)
  library(data.table)
})

gse <- getGEO("GSE168807",GSEMatrix=TRUE)

scrna_gsms <- c("GSM5171988","GSM5171989","GSM5171990","GSM5171991", "GSM5257698", "GSM5257699", "GSM5257700")
outdir <- "GSE168807_scRNA_mtx"
dir.create(outdir, showWarnings = FALSE)

# download per-GSM 10x files
dl <- rbindlist(lapply(scrna_gsms, function(gsm) {
  cat("Fetching", gsm, "...\n")
  df <- getGEOSuppFiles(gsm, baseDir = outdir, makeDirectory = TRUE)
  setDT(df, keep.rownames = "filename")[, gsm := gsm]
}), fill = TRUE)































