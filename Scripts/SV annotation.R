###############################
# SVIM-asm Structural Variant
# Biological Impact Analysis
###############################

## ---- 1. Load required packages ----
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "VariantAnnotation",
  "GenomicRanges",
  "rtracklayer"
), ask = FALSE, update = FALSE)

install.packages(c("ggplot2"))

library(VariantAnnotation)
library(GenomicRanges)
library(rtracklayer)
library(ggplot2)

## ---- 2. Input files (EDIT PATHS IF NEEDED) ----
vcf_file  <- "variants.vcf"
gff_file  <- "GCF_000006945.2_ASM694v2_genomic.gff"

## ---- 3. Read SVIM-asm VCF ----
vcf <- readVcf(vcf_file)
sv_gr <- rowRanges(vcf)

## ---- 4. Load reference genome annotation ----
annotation <- import(gff_file)

# Keep only gene features
genes <- annotation[annotation$type == "gene"]

## ---- 5. Overlap SVs with genes ----
hits <- findOverlaps(sv_gr, genes)

# Default classification
sv_gr$impact <- "intergenic"
sv_gr$impact[queryHits(hits)] <- "genic"

# Attach gene names if available
sv_gr$gene <- NA
if ("Name" %in% names(mcols(genes))) {
  sv_gr$gene[queryHits(hits)] <- genes$Name[subjectHits(hits)]
}

## ---- 6. Build results table ----
sv_table <- data.frame(
  chromosome = as.character(seqnames(sv_gr)),
  start      = start(sv_gr),
  end        = end(sv_gr),
  sv_type    = info(vcf)$SVTYPE,
  sv_length  = info(vcf)$SVLEN,
  impact     = sv_gr$impact,
  gene       = sv_gr$gene,
  stringsAsFactors = FALSE
)

## ---- 7. Save results ----
write.csv(
  sv_table,
  file = "SV_biological_impact_summary.csv",
  row.names = FALSE
)

## ---- 8. Summary statistics ----
cat("\nSV counts by type:\n")
print(table(sv_table$sv_type))

cat("\nGenic vs intergenic SVs:\n")
print(table(sv_table$impact))

cat("\nGenes affected by SVs:\n")
print(unique(na.omit(sv_table$gene)))

## ---- 9. Visualization ----
p <- ggplot(sv_table, aes(x = sv_type, fill = impact)) +
  geom_bar(color = "black") +
  theme_minimal() +
  labs(
    title = "Structural Variants by Type and Genomic Impact",
    x = "Structural Variant Type",
    y = "Count"
  )
p
ggsave(
  filename = "SV_impact_barplot.png",
  plot = p,
  width = 7,
  height = 5
)

cat("\nAnalysis complete.\n")
