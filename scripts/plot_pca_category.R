#!/usr/bin/env Rscript

# ---------------------------
# Load packages
# ---------------------------
suppressMessages(library(optparse))
suppressMessages(library(ggplot2))
suppressMessages(library(ggrepel))

# ---------------------------
# Options
# ---------------------------
option_list <- list(
  make_option("--pca", type="character"),
  make_option("--metadata", type="character"),
  make_option("--pca_col", type="character"),
  make_option("--out", type="character")
)

opt <- parse_args(OptionParser(option_list = option_list))

# ---------------------------
# Inputs
# ---------------------------
pca_file <- opt$pca
metadata_file <- opt$metadata

# ---------------------------
# Parse inputs
# ---------------------------
pca_cols <- as.character(trimws(strsplit(opt$pca_col, ",")[[1]]))

output_file <- opt$out
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

# ---------------------------
# Read PCA
# ---------------------------
pca_df <- read.csv(pca_file, stringsAsFactors = FALSE, check.names = FALSE)

if (!all(c("PC1", "PC2", "Barcode") %in% colnames(pca_df))) {
  stop("PCA file must contain PC1, PC2, Barcode")
}

# ---------------------------
# Read metadata
# ---------------------------
#lines <- readLines(metadata_file)
#start_idx <- which(grepl("Pool", lines))[1]

#if (is.na(start_idx)) {
#  stop("No metadata header found (Pool missing)")
#}

metadata <- read.csv(
  metadata_file,
 # skip = start_idx,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

metadata <- metadata[!is.na(metadata$Barcode) & metadata$Barcode != "", ]
metadata$Barcode <- make.unique(trimws(metadata$Barcode))
rownames(metadata) <- metadata$Barcode
metadata$Barcode <- NULL

colnames(metadata) <- gsub(" ", "_", colnames(metadata))

# ---------------------------
# Loop over PCA columns
# ---------------------------
for (col in pca_cols) {

  col <- trimws(col)
  col_safe <- gsub(" ", "_", col)

  if (!(col_safe %in% colnames(metadata))) {
    warning(paste("Missing metadata column:", col, "→ skipping"))
    next
  }

  vals <- metadata[pca_df$Barcode, col_safe]
  vals <- as.character(trimws(vals))

  if (all(is.na(vals))) {
    warning(paste("All values NA in column:", col, "→ skipping"))
    next
  }

  # ---------------------------
  # SIMPLE FACTOR MAPPING ONLY
  # ---------------------------
  pca_df[[col_safe]] <- as.factor(vals)

  # ---------------------------
  # Output file
  # ---------------------------
  out_file <- output_file
  out_file <- gsub("\\{col\\}", col_safe, out_file)

  # safety cleanup if old templates remain
  out_file <- gsub("\\{ref\\}", "", out_file)
  out_file <- gsub("--", "-", out_file)

  # ---------------------------
  # Plot
  # ---------------------------
  p <- ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = .data[[col_safe]]), size = 3) +
  geom_text_repel(
  aes(label = Barcode),
  size = 2.5,
  max.overlaps = Inf,
  box.padding = 0.4,
  point.padding = 0.2,
  min.segment.length = 0
) +
  theme_minimal() +
  labs(
    title = paste("PCA colored by", col),
    x = "PC1",
    y = "PC2",
    color = col
  ) +
  theme(legend.position = "bottom")

  ggsave(out_file, p, width = 7, height = 8, dpi = 150, bg = "white")

  cat("✅ PCA saved:", out_file, "| column:", col, "\n")
}