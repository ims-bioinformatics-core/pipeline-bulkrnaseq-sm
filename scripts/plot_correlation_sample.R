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

n_samples <- ncol(cor_mat)

# ---------------------------
# Label font size
# ---------------------------
#
# Use readable font sizes.
# Do not reduce to 4-5 pt.
#

if (n_samples <= 30) {

  label_size <- 10

} else if (n_samples <= 60) {

  label_size <- 8

} else if (n_samples <= 100) {

  label_size <- 7

} else {

  label_size <- 6

}

# ---------------------------
# Figure dimensions
# ---------------------------
#
# Moderate output size.
# Avoid excessive scaling with n_samples.
#

fig_width <- 3000
fig_height <- 2800

# ---------------------------
# Correlation colour scale
# ---------------------------

min_cor <- min(
  cor_mat,
  na.rm = TRUE
)

max_cor <- max(
  cor_mat,
  na.rm = TRUE
)

mid_cor <- (
  min_cor +
    max_cor
) / 2

# ---------------------------
# Output
# ---------------------------

png(
  filename = output_file,
  width = fig_width,
  height = fig_height,
  res = 150
)

# ---------------------------
# Heatmap
# ---------------------------

ht <- Heatmap(

  cor_mat,

  name = "Pearson\ncorrelation",

  # -------------------------
  # Group annotation
  # -------------------------

  top_annotation = ha,

  # -------------------------
  # Colour scale
  # -------------------------

  col = colorRamp2(
    c(
      min_cor,
      mid_cor,
      max_cor
    ),
    c(
      "blue",
      "white",
      "red"
    )
  ),

  # -------------------------
  # Clustering
  # -------------------------

  cluster_rows = TRUE,
  cluster_columns = TRUE,

  # -------------------------
  # BOTH ROW AND COLUMN LABELS
  # -------------------------

  show_row_names = TRUE,
  show_column_names = TRUE,

  # -------------------------
  # FULL LABELS
  # NO WRAPPING
  # -------------------------

  row_labels = rownames(cor_mat),
  column_labels = colnames(cor_mat),

  # -------------------------
  # ROW LABELS
  # -------------------------

  row_names_gp = gpar(
    fontsize = label_size
  ),

  row_names_rot = 0,

  # -------------------------
  # COLUMN LABELS
  # -------------------------

  column_names_gp = gpar(
    fontsize = label_size
  ),

  column_names_rot = 45,

  # -------------------------
  # Heatmap borders
  # -------------------------

  border = TRUE,

  # -------------------------
  # Legends
  # -------------------------

  heatmap_legend_param = list(

    title_gp = gpar(
      fontsize = 12,
      fontface = "bold"
    ),

    labels_gp = gpar(
      fontsize = 10
    ),

    grid_width = unit(
      6,
      "mm"
    ),

    grid_height = unit(
      6,
      "mm"
    )

  )

)

# ---------------------------
# Draw heatmap
# ---------------------------

draw(

  ht,

  heatmap_legend_side = "right",

  annotation_legend_side = "right",

  padding = unit(
    c(10, 35, 15, 15),
    "mm"
  )

)

dev.off()

cat(
  "✅ Correlation heatmap saved to:",
  output_file,
  "\n"
)