#!/usr/bin/env Rscript

suppressMessages(library(optparse))
suppressMessages(library(tidyverse))
suppressMessages(library(ggrepel))

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

metadata <- metadata[
  !is.na(rownames(metadata)) & rownames(metadata) != "",
  ,
  drop = FALSE
]

# ---------------------------
# Format counts
# ---------------------------
rownames(counts) <- trimws(counts[[1]])
counts <- counts[, -1, drop = FALSE]
colnames(counts) <- trimws(colnames(counts))
rownames(counts) <- sub("\\..*$", "", rownames(counts))

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
# Annotation (STRICT protein-coding filter + intersection fix)
# ---------------------------
cat("ℹ️ Reading annotation...\n")

colnames(annot) <- trimws(colnames(annot))

required <- c("ensembl_gene_id", "gene_symbol", "protein_coding")

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

# remove invalid IDs
annot <- annot[
  !is.na(annot$ensembl_gene_id) &
  annot$ensembl_gene_id != "",
  ,
  drop = FALSE
]

# ---------------------------
# HARD FILTER: protein coding only
# ---------------------------
annot <- annot[annot$protein_coding == TRUE, , drop = FALSE]

cat("ℹ️ Protein-coding genes in annotation:", nrow(annot), "\n")

if (nrow(annot) == 0) {
  stop("❌ No protein-coding genes found in annotation")
}

# ---------------------------
# 🔥 CRITICAL FIX: intersect with counts
# ---------------------------
common_genes <- intersect(rownames(counts), annot$ensembl_gene_id)

if (length(common_genes) == 0) {
  stop("❌ No overlap between counts and protein-coding annotation")
}

counts <- counts[common_genes, , drop = FALSE]
annot  <- annot[annot$ensembl_gene_id %in% common_genes, , drop = FALSE]

cat("ℹ️ Genes after filtering:", nrow(counts), "\n")

# ---------------------------
# Fill missing symbols
# ---------------------------
annot$gene_symbol[
  is.na(annot$gene_symbol) | annot$gene_symbol == ""
] <- annot$ensembl_gene_id[
  is.na(annot$gene_symbol) | annot$gene_symbol == ""
]

id_to_symbol <- setNames(
  annot$gene_symbol,
  annot$ensembl_gene_id
)

# ---------------------------
# Log transform
# ---------------------------
counts_log <- log2(as.matrix(counts) + 1)

symbols <- id_to_symbol[rownames(counts_log)]
symbols[is.na(symbols) | symbols == ""] <- rownames(counts_log)

rownames(counts_log) <- make.unique(symbols)

# ---------------------------
# Group column
# ---------------------------
group_col <- opt$dge_cat

if (!group_col %in% colnames(metadata)) {
  stop(paste0("Group column not found: ", group_col))
}

meta_group <- metadata[[group_col]]
names(meta_group) <- rownames(metadata)

meta_group <- meta_group[!is.na(meta_group)]

metadata   <- metadata[names(meta_group), , drop = FALSE]
counts_log <- counts_log[, names(meta_group), drop = FALSE]

# ---------------------------
# Check groups
# ---------------------------
groups <- unique(meta_group)

if (length(groups) < 2) {
  stop("Need at least 2 groups for MA plot")
}

groups <- sort(groups)
group1 <- groups[1]
group2 <- groups[2]

samples1 <- names(meta_group)[meta_group == group1]
samples2 <- names(meta_group)[meta_group == group2]

# ---------------------------
# MA calculation
# ---------------------------
mean1 <- rowMeans(counts_log[, samples1, drop = FALSE], na.rm = TRUE)
mean2 <- rowMeans(counts_log[, samples2, drop = FALSE], na.rm = TRUE)

ma_df <- data.frame(
  A = (mean1 + mean2) / 2,
  M = mean2 - mean1,
  gene = rownames(counts_log),
  stringsAsFactors = FALSE
)

ma_df <- ma_df[is.finite(ma_df$A) & is.finite(ma_df$M), ]

# ---------------------------
# Top genes
# ---------------------------
top_pos <- ma_df[order(ma_df$M, decreasing = TRUE), ][1:min(20, nrow(ma_df)), ]
top_neg <- ma_df[order(ma_df$M, decreasing = FALSE), ][1:min(20, nrow(ma_df)), ]
top_genes <- rbind(top_pos, top_neg)

# ---------------------------
# Plot
# ---------------------------
p <- ggplot(ma_df, aes(x = A, y = M)) +

  geom_point(alpha = 0.15, size = 0.4) +

  geom_point(data = top_genes, color = "red", size = 1.2) +

  geom_text_repel(
    data = top_genes,
    aes(label = gene),
    size = 1,
    color = "red",
    max.overlaps = Inf
  ) +

  geom_hline(yintercept = 0, color = "red") +

  geom_smooth(method = "loess", se = FALSE, color = "blue") +

  labs(
    title = paste("MA plot:", group2, "vs", group1),
    x = "A (mean expression)",
    y = "M (log fold change)"
  ) +

  theme_classic(base_size = 8)

ggsave(opt$out, plot = p, width = 4, height = 4, dpi = 300)