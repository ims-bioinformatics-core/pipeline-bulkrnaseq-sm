#!/usr/bin/env Rscript

# ---------------------------
# Load packages
# ---------------------------
suppressMessages(library(optparse))
suppressMessages(library(tidyverse))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(circlize))

# ---------------------------
# Command-line options
# ---------------------------
option_list <- list(
  make_option("--counts", type="character"),
  make_option("--metadata", type="character"),
  make_option("--annotation", type="character"),
  make_option("--out", type="character"),
  make_option("--group", type="character")
)

opt <- parse_args(OptionParser(option_list = option_list))

counts_file     <- opt$counts
metadata_file   <- opt$metadata
annot_file      <- opt$annotation
output_file     <- opt$out
group           <- opt$group

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

# ---------------------------
# 1. Metadata
# ---------------------------
metadata <- read.csv(metadata_file, stringsAsFactors = FALSE, check.names = FALSE)

metadata <- metadata[!is.na(metadata$Barcode) & metadata$Barcode != "", ]
metadata$Barcode <- make.unique(trimws(metadata$Barcode))
rownames(metadata) <- metadata$Barcode
metadata$Barcode <- NULL

# ---------------------------
# 2. Counts (Ensembl IDs)
# ---------------------------
counts <- read.csv(counts_file, stringsAsFactors = FALSE, check.names = FALSE)

rownames(counts) <- counts[[1]]
counts <- counts[, -1, drop = FALSE]
colnames(counts) <- trimws(colnames(counts))

# remove Ensembl version
rownames(counts) <- sub("\\..*$", "", rownames(counts))

# sample match
common_samples <- intersect(colnames(counts), rownames(metadata))

if (length(common_samples) < 2) {
  stop("❌ Not enough overlapping samples between counts and metadata")
}

counts <- counts[, common_samples, drop = FALSE]
metadata <- metadata[common_samples, , drop = FALSE]

# ---------------------------
# 3. Read annotation (NEW - replaces GTF)
# ---------------------------
cat("ℹ️ Reading gene annotation...\n")

annot <- read.delim(
  annot_file,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

colnames(annot) <- trimws(colnames(annot))

required <- c("ensembl_gene_id", "gene_symbol")

if (!all(required %in% colnames(annot))) {
  stop(
    "❌ Annotation must contain: ",
    paste(required, collapse = ", "),
    "\nFound: ",
    paste(colnames(annot), collapse = ", ")
  )
}

annot <- annot %>%
  mutate(
    ensembl_gene_id = trimws(ensembl_gene_id),
    gene_symbol = trimws(gene_symbol)
  ) %>%
  filter(ensembl_gene_id != "")

# keep only protein coding if column exists
if ("protein_coding" %in% colnames(annot)) {
  annot <- annot[annot$protein_coding == TRUE, ]
}

# ---------------------------
# 4. Filter counts to protein coding genes
# ---------------------------
protein_gene_ids <- unique(annot$ensembl_gene_id)

counts <- counts[rownames(counts) %in% protein_gene_ids, , drop = FALSE]

cat("✅ Counts filtered to", nrow(counts), "protein-coding genes.\n")

if (nrow(counts) == 0) {
  stop("❌ No protein-coding genes remaining after filtering.")
}

# ---------------------------
# 5. Normalize expression
# ---------------------------
counts_log <- log2(as.matrix(counts) + 1)

# ---------------------------
# 6. Correlation matrix
# ---------------------------
cor_mat <- cor(
  counts_log,
  method = "pearson",
  use = "pairwise.complete.obs"
)

# ---------------------------
# 7. Annotation (group)
# ---------------------------
if (is.null(group) || !(group %in% colnames(metadata))) {
  stop("❌ group column not found in metadata: ", group)
}

meta_group <- trimws(metadata[[group]])
names(meta_group) <- rownames(metadata)
meta_group <- factor(meta_group)

group_colors <- setNames(
  colorRampPalette(c(
    "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e"
  ))(length(levels(meta_group))),
  levels(meta_group)
)

ha <- HeatmapAnnotation(
  Group = meta_group,
  col = list(Group = group_colors),
  show_annotation_name = TRUE
)

# ---------------------------
# 8. Plot heatmap
# ---------------------------
png(output_file, width = 900, height = 800, res = 150)

min_cor <- min(cor_mat, na.rm = TRUE)
max_cor <- max(cor_mat, na.rm = TRUE)

Heatmap(
  cor_mat,
  name = "Pearson\ncorrelation",
  top_annotation = ha,
  col = colorRamp2(
    c(min_cor, (min_cor + max_cor)/2, max_cor),
    c("blue", "white", "red")
  ),
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 9),
  column_names_gp = gpar(fontsize = 9, rot = 45)
)
dev.off()

cat("✅ Correlation heatmap saved to:", output_file, "\n")