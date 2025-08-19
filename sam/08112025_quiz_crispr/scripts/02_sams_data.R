library(data.table)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(patchwork)

dt <- fread("inputs/Sam_quiz_v_veh.tsv")

dt <- dt[, c("gene_name", "log2FoldChange", "pvalue", "padj"), with = FALSE]

dt[, logP := -log10(pvalue)]
dt[, significant := padj < 0.05]

dt$gene_name_upper <- dt$gene_name %>% toupper()

saveRDS(dt, "intermediates/sams_data_clean.rds")

# Volcano plot
pdf("outputs/sam_volc.pdf")
	ggplot() + 
		geom_point(
			data = dt[significant == TRUE], 
			aes(x = log2FoldChange, y = logP),
			color = "firebrick3"
		) + 
		geom_point(
			data = dt[significant == FALSE], 
			aes(x = log2FoldChange, y = logP),
			color = "grey70"
		) + 
		labs(x = "Log2(Fold Change)", y = "-log(P)", title = "Quizartinib v. Vehicle")
dev.off()

dna_repair_genes <- c(
  "AAAS","ADA","ADCY6","ADRM1","AK1","AK3","APRT","ARL6IP1","BCAM","BCAP31",
  "BOLA2","BRF2","MPC2","CANT1","CCNO","CDA","CETN2","CLP1","CMPK2","NELFB",
  "COX17","CSTF3","DAD1","DCTN4","DDB1","DDB2","GSDME","DGCR8","DGUOK","DUT",
  "EDF1","EIF1B","AGO4","ELL","ERCC1","ERCC2","ERCC3","ERCC4","ERCC5","ERCC8",
  "FEN1","GMPR2","GPX4","GTF2A2","GTF2B","GTF2F1","GTF2H1","GTF2H3","GTF2H5",
  "GTF3C5","GUK1","HCLS1","HPRT1","IMPDH2","ITPA","LIG1","MPG","MRPL40","NCBP2",
  "NFX1","NME1","NME3","NME4","NPR2","NT5C","NT5C3A","NUDT21","NUDT9","PCNA",
  "PDE4B","PDE6G","PNP","POLA1","POLA2","POLB","POLD1","POLD3","POLD4","POLE4",
  "POLH","POLL","POLR1C","POLR1D","POLR2A","POLR2C","POLR2D","POLR2E","POLR2F",
  "POLR2G","POLR2H","POLR2I","POLR2J","POLR2K","POLR3C","POLR3GL","POM121",
  "PRIM1","RAD51","RAD52","RAE1","RALA","RBX1","NELFE","REV3L","RFC2","RFC3",
  "RFC4","RFC5","RNMT","RPA2","RPA3","RRM2B","SAC3D1","SDCBP","SEC61A1","SF3A3",
  "SMAD5","SNAPC4","SNAPC5","SRSF6","SSRP1","STX3","SUPT4H1","SUPT5H","SURF1",
  "TAF10","TAF12","TAF13","TAF1C","TAF6","TAF9","TARBP2","ELOA","NELFCD",
  "ALYREF","TK2","TMED2","TP53","TSG101","TYMS","UMPS","UPF3B","USP11","VPS28",
  "VPS37B","VPS37D","XPC","ZNF707","POLR1H","ZWINT"
)

dna_repair_genes <- data.frame(
  gene = c(
    # BER
    "OGG1", "MUTYH", "APEX1",
    # NER
    "XPA", "XPC", "ERCC1", "ERCC2",
    # MMR
    "MSH2", "MSH6", "MLH1", "PMS2",
    # DSB Repair
    "BRCA1", "BRCA2", "RAD51", "XRCC6", "XRCC5",
    # Signaling
    "TP53", "ATM", "ATR"
  ),
  pathway = c(
    # BER
    rep("Base Excision Repair", 3),
    # NER
    rep("Nucleotide Excision Repair", 4),
    # MMR
    rep("Mismatch Repair", 4),
    # DSB Repair
    rep("Double-Strand Break Repair", 5),
    # Signaling
    rep("DNA Damage Signaling", 3)
  ),
  stringsAsFactors = FALSE
)

dt[gene_name_upper %in% dna_repair_genes$gene][significant == TRUE]

dt_m <- merge(
  dt, 
  dna_repair_genes,
  by.x = "gene_name_upper",
  by.y = "gene",
  all.x = TRUE
)

pdf("outputs/sam_volc2.pdf")
  ggplot() + 
    geom_point(
      data = dt_m[significant == TRUE & is.na(pathway)], 
      aes(x = log2FoldChange, y = logP),
      color = "firebrick3"
    ) + 
    geom_point(
      data = dt[significant == FALSE], 
      aes(x = log2FoldChange, y = logP),
      color = "grey70"
    ) + 
    geom_point(
      data = dt_m[!is.na(pathway)], 
      aes(x = log2FoldChange, y = logP, color = pathway)
      # color = pathway
    ) + 
        # labels
    geom_text_repel(
      data = dt_m[!is.na(pathway)],
      aes(x = log2FoldChange, y = logP, label = gene_name_upper),
      size = 4,
      nudge_x = ifelse(dt_m[!is.na(pathway)]$log2FoldChange < 0, -0.015, 0.015),
      min.segment.length = 0,
      box.padding = 0.25,
      point.padding = 0.15,
      max.overlaps = Inf,
      segment.size = 0.4
    ) + 
    labs(
      x = "Log2(Fold Change)", 
      y = "-log(P)", 
      title = "Quizartinib v. Vehicle -- DNA Repair Genes",
      color = "Repair Pathway"
    )
dev.off()




