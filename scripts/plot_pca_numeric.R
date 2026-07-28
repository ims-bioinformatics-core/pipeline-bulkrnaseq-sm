#!/usr/bin/env Rscript

# ---------------------------
# Load packages
# ---------------------------
suppressMessages(library(optparse))
suppressMessages(library(ggplot2))
suppressMessages(library(viridis))

# ---------------------------
# Command-line options
# ---------------------------
option_list <- list(
  make_option("--pca", type="character", help="Path to PCA table CSV (must include PC1, PC2, Barcode)"),
  make_option("--metadata", type="character", help="Path to metadata CSV"),
  make_option("--pca_col", type="character", help="Metadata column(s), comma-separated"),
  make_option("--out", type="character", help="Output file pattern (use {col})")
)

opt <- parse_args(OptionParser(option_list = option_list))

pca_file <- opt$pca
metadata_file <- opt$metadata
pca_cols <- strsplit(opt$pca_col, ",")[[1]] |> trimws()
output_file <- opt$out

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

# ---------------------------
# Read PCA table
# ---------------------------
pca_df <- read.csv(pca_file, stringsAsFactors = FALSE, check.names = FALSE)

if (!all(c("PC1", "PC2", "Barcode") %in% colnames(pca_df))) {
  stop("PCA file must contain columns: PC1, PC2, Barcode")
}

# ---------------------------
# Read metadata
# ---------------------------

#lines <- readLines(metadata_file)
#start_idx <- which(grepl("Pool", lines))[1]

#if (is.na(start_idx)) {
#  stop("No line containing 'Pool' found in metadata file!")
#}

metadata <- read.csv(
  metadata_file,
#  skip = start_idx,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Remove empty Barcode rows
metadata <- metadata[!is.na(metadata$Barcode) & metadata$Barcode != "", ]

# Clean and set rownames
metadata$Barcode <- make.unique(trimws(metadata$Barcode))
rownames(metadata) <- metadata$Barcode
metadata$Barcode <- NULL

# Replace spaces in column names with underscores
colnames(metadata) <- gsub(" ", "_", colnames(metadata))

# ---------------------------PCA file must contain
# Plot PCA per metadata column
# ---------------------------
for (col in pca_cols) {

  col_safe <- gsub(" ", "_", col)

  if (!(col_safe %in% colnames(metadata))) {
    stop(paste("Missing metadata column:", col))
  }

  # join metadata to PCA
  pca_df[[col_safe]] <- metadata[pca_df$Barcode, col_safe]

  # convert to numeric if possible
  suppressWarnings(is_num <- !any(is.na(as.numeric(pca_df[[col_safe]])) & !is.na(pca_df[[col_safe]])))

  if (is_num) {
    pca_df[[col_safe]] <- as.numeric(pca_df[[col_safe]])

    fill_scale <- scale_color_viridis_c(option = "viridis", na.value = "grey50")

  } else {
    pca_df[[col_safe]] <- as.factor(pca_df[[col_safe]])

    fill_scale <- scale_color_viridis_d(option = "viridis", na.value = "grey50")
  }

  out_file <- gsub("\\{col\\}", col_safe, output_file)

  p <- ggplot(pca_df, aes(x = PC1, y = PC2)) +
    geom_point(aes(color = .data[[col_safe]]), size = 3) +
    geom_text(aes(label = Barcode), vjust = -0.5, size = 2.5) +
    fill_scale +
    theme_minimal() +
    labs(
      title = paste("PCA colored by", col),
      x = "PC1",
      y = "PC2",
      color = col
    ) +
    theme(legend.position = "bottom")

  ggsave(out_file, p, width = 7, height = 8, dpi = 150, bg = "white")

  cat("✅ PCA saved to:", out_file, "\n")
}