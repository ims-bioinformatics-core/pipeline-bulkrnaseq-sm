#!/usr/bin/env Rscript

suppressMessages(library(optparse))
suppressMessages(library(pheatmap))

# ---------------------------
# CLI arguments
# ---------------------------
option_list <- list(
  make_option("--counts", type="character"),
  make_option("--metadata", type="character"),
  make_option("--annotation", type="character"),
  make_option("--dge_cat", type="character"),
  make_option("--out", type="character")
)

opt <- parse_args(OptionParser(option_list = option_list))

# ---------------------------
# Load data
# ---------------------------
metadata <- read.csv(opt$metadata, stringsAsFactors = FALSE, check.names = FALSE)
counts   <- read.csv(opt$counts, stringsAsFactors = FALSE, check.names = FALSE)
annot    <- read.delim(opt$annotation, stringsAsFactors = FALSE, check.names = FALSE)

# ---------------------------
# Format metadata
# ---------------------------
rownames(metadata) <- trimws(metadata$Barcode)
metadata$Barcode <- NULL
metadata <- metadata[!is.na(rownames(metadata)) & rownames(metadata) != "", , drop = FALSE]

# ---------------------------
# Format counts
# ---------------------------
rownames(counts) <- trimws(counts[[1]])
counts <- counts[, -1, drop = FALSE]
colnames(counts) <- trimws(colnames(counts))
rownames(counts) <- sub("\\..*$", "", rownames(counts))  # remove Ensembl version

# ---------------------------
# Match samples
# ---------------------------
common <- intersect(colnames(counts), rownames(metadata))

if (length(common) < 2) {
  stop("❌ Not enough overlapping samples")
}

counts   <- counts[, common, drop = FALSE]
metadata <- metadata[common, , drop = FALSE]

# ---------------------------
# Read annotation (NEW)
# ---------------------------
cat("ℹ️ Reading annotation...\n")

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

annot$ensembl_gene_id <- trimws(annot$ensembl_gene_id)
annot$gene_symbol     <- trimws(annot$gene_symbol)

annot <- annot[annot$ensembl_gene_id != "", , drop = FALSE]

# fill missing symbols safely
annot$gene_symbol[
  is.na(annot$gene_symbol) | annot$gene_symbol == ""
] <- annot$ensembl_gene_id[
  is.na(annot$gene_symbol) | annot$gene_symbol == ""
]

# ---------------------------
# protein-coding filter
# ---------------------------
if ("protein_coding" %in% colnames(annot)) {
  annot <- annot[annot$protein_coding == TRUE, , drop = FALSE]
}

protein_genes <- unique(annot$ensembl_gene_id)

cat("ℹ️ Protein-coding genes:", length(protein_genes), "\n")

# ---------------------------
# Filter counts
# ---------------------------
counts <- counts[rownames(counts) %in% protein_genes, , drop = FALSE]

if (nrow(counts) == 0) {
  stop("❌ No protein-coding genes found in counts")
}

# ---------------------------
# Log transform
# ---------------------------
counts_log <- log2(as.matrix(counts) + 1)

# ---------------------------
# Gene symbol mapping
# ---------------------------
id_to_symbol <- setNames(annot$gene_symbol, annot$ensembl_gene_id)

symbols <- id_to_symbol[rownames(counts_log)]
symbols[is.na(symbols) | symbols == ""] <- rownames(counts_log)

rownames(counts_log) <- make.unique(symbols)

# ---------------------------
# Grouping
# ---------------------------
group_col <- opt$dge_cat

if (!group_col %in% colnames(metadata)) {
  stop(paste0("Group column not found: ", group_col))
}

groups <- metadata[[group_col]]
names(groups) <- rownames(metadata)

groups <- groups[common]

counts_log <- counts_log[, names(groups), drop = FALSE]

# ---------------------------
# TOP VARIABLE GENES PER GROUP
# ---------------------------
selected_genes <- c()

for (g in unique(groups)) {

  samples <- names(groups)[groups == g]

  if (length(samples) < 2) next

  sub <- counts_log[, samples, drop = FALSE]

  gene_var <- apply(sub, 1, var, na.rm = TRUE)
  gene_var[is.na(gene_var)] <- 0

  top_genes <- names(sort(gene_var, decreasing = TRUE))[1:min(50, length(gene_var))]

  selected_genes <- unique(c(selected_genes, top_genes))
}

if (length(selected_genes) == 0) {
  stop("❌ No genes selected")
}

# ---------------------------
# MATRIX
# ---------------------------
mat <- counts_log[selected_genes, , drop = FALSE]
mat <- t(scale(t(mat)))
mat[is.na(mat)] <- 0

# ---------------------------
# SAMPLE ORDERING
# ---------------------------
sample_order <- c()

for (g in unique(groups)) {

  samples_g <- names(groups)[groups == g]

  if (length(samples_g) > 1) {
    sub_mat <- mat[, samples_g, drop = FALSE]
    hc <- hclust(dist(t(sub_mat)))
    samples_g <- samples_g[hc$order]
  }

  sample_order <- c(sample_order, samples_g)
}

mat <- mat[, sample_order, drop = FALSE]

# ---------------------------
# ANNOTATION
# ---------------------------
annotation_col <- data.frame(Group = groups[sample_order])
rownames(annotation_col) <- sample_order

# ---------------------------
# HEATMAP
# ---------------------------
pheatmap(
  mat,
  annotation_col = annotation_col,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 6,
  fontsize_col = 6,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  filename = opt$out,
  width = 10,
  height = 6
)