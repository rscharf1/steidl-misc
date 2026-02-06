suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(data.table)
  library(apeglm)
})

counts_file <- "../03_align/counts/counts.txt"
meta_file   <- "sample_metadata.csv"

counts <- read.delim(
  counts_file,
  comment.char = "#",
  check.names = FALSE
)

# featureCounts columns:
# Geneid Chr Start End Strand Length sample1 sample2 ...
count_matrix <- counts %>%
  select(-Chr, -Start, -End, -Strand, -Length) %>%
  column_to_rownames("Geneid")

colnames(count_matrix) <- colnames(count_matrix) %>% 
  basename() %>% 
  gsub(".Aligned.sortedByCoord.out.bam", "", .)

sample_names <- fread(meta_file)

coldata <- as.data.frame(sample_names)
rownames(coldata) <- coldata$sample_id
coldata$sample_id <- NULL

coldata$condition <- factor(coldata$condition)

count_matrix <- count_matrix[, rownames(coldata)]

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData   = coldata,
  design    = ~ condition
)

dds <- DESeq(dds)

res <- results(
  dds,
  contrast = c("condition", "AraC_1_day", "control")
)

summary(res)

resultsNames(dds)

resLFC <- lfcShrink(
  dds,
  coef = "condition_AraC_1_day_vs_control",
  type = "apeglm"
)

resLFC <- resLFC[order(resLFC$padj), ]

write.csv(
  as.data.frame(resLFC),
  file = "deseq2_AraC_vs_control.csv"
)






