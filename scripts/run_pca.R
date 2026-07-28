#!/usr/bin/env Rscript

# ---------------------------
# Load packages
# ---------------------------
suppressMessages(library(optparse))
suppressMessages(library(matrixStats))
suppressMessages(library(DESeq2))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))

# ---------------------------
# Command-line options
# ---------------------------
option_list <- list(
  make_option("--counts", type="character"),
  make_option("--metadata", type="character"),
  make_option("--annotation", type="character"),
  make_option("--out", type="character"),
  make_option("--norm", type="character", default="vst"),
  make_option("--filt", type="logical", default=TRUE)
)

opt <- parse_args(OptionParser(option_list = option_list))

counts_file   <- opt$counts
metadata_file <- opt$metadata
annot_file    <- opt$annotation
output_file   <- opt$out
norm_method   <- tolower(opt$norm)
apply_filter  <- opt$filt

if (!norm_method %in% c("vst", "cpm")) {
  stop("Invalid normalization method. Choose 'vst' or 'cpm'.")
}

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

# match samples
common_samples <- intersect(colnames(counts), rownames(metadata))

if (length(common_samples) < 2) {
  stop("❌ Not enough overlapping samples")
}

counts <- counts[, common_samples, drop = FALSE]
metadata <- metadata[common_samples, , drop = FALSE]

# ---------------------------
# 3. Read annotation (NEW)
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

protein_gene_ids <- unique(annot$ensembl_gene_id)

cat("ℹ️ Protein-coding genes in annotation:", length(protein_gene_ids), "\n")

# ---------------------------
# 4. Filter counts
# ---------------------------
counts <- counts[rownames(counts) %in% protein_gene_ids, , drop = FALSE]

if (nrow(counts) == 0) {
  stop("❌ No protein-coding genes after filtering.")
}

cat("✅ Counts filtered to", nrow(counts), "protein-coding genes.\n")

# ---------------------------
# 5. Normalize counts
# ---------------------------
if (norm_method == "vst") {

  dds <- DESeqDataSetFromMatrix(
    countData = round(counts),
    colData = metadata,
    design = ~1
  )

  dds <- estimateSizeFactors(dds)

  if (apply_filter) {
    keep_genes <- rowSums(counts(dds) >= 10) >= 3
    dds <- dds[keep_genes, ]
    cat("ℹ️ Genes after filtering:", nrow(dds), "\n")

    if (nrow(dds) < 2) stop("Too few genes after filtering.")
  }

  norm_counts <- assay(vst(dds, blind = TRUE))

} else {

  lib_sizes <- colSums(counts)
  cpm <- t(t(counts) / lib_sizes * 1e6)
  norm_counts <- log2(cpm + 1)

  if (apply_filter) {
    keep_genes <- rowSums(norm_counts >= 1) >= 3
    norm_counts <- norm_counts[keep_genes, , drop = FALSE]
  }
}

# ---------------------------
# 6. PCA
# ---------------------------
top_var <- head(order(rowVars(norm_counts), decreasing = TRUE), 500)
mat <- norm_counts[top_var, , drop = FALSE]

pca <- prcomp(t(mat), center = TRUE, scale. = FALSE)

npcs <- min(10, ncol(pca$x))
pca_df <- as.data.frame(pca$x[, 1:npcs, drop = FALSE])
pca_df$Barcode <- rownames(pca_df)

# ---------------------------
# 7. Save output
# ---------------------------
pca_out_file <- file.path(dirname(output_file), "pca.csv")
write.csv(pca_df, pca_out_file, quote = FALSE)

cat("✅ PCA saved to:", pca_out_file, "\n")