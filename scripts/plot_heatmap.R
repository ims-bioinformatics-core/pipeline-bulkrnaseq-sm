
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
metadata <- read.csv(
  opt$metadata,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

counts <- read.csv(
  opt$counts,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

annot <- read.delim(
  opt$annotation,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

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

# Remove Ensembl version
rownames(counts) <- sub("\\..*$", "", rownames(counts))

# ---------------------------
# Match samples
# ---------------------------
common <- intersect(colnames(counts), rownames(metadata))

if (length(common) < 2) {
  stop("❌ Not enough overlapping samples")
}

counts <- counts[, common, drop = FALSE]
metadata <- metadata[common, , drop = FALSE]

# ---------------------------
# Read annotation
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
annot$gene_symbol <- trimws(annot$gene_symbol)

annot <- annot[
  annot$ensembl_gene_id != "",
  ,
  drop = FALSE
]

# Fill missing symbols safely
missing_symbol <- (
  is.na(annot$gene_symbol) |
  annot$gene_symbol == ""
)

annot$gene_symbol[missing_symbol] <-
  annot$ensembl_gene_id[missing_symbol]

# ---------------------------
# Protein-coding filter
# ---------------------------
if ("protein_coding" %in% colnames(annot)) {

  # Handle logical TRUE/FALSE or character values
  if (is.logical(annot$protein_coding)) {
    annot <- annot[annot$protein_coding == TRUE, , drop = FALSE]
  } else {
    annot <- annot[
      tolower(trimws(annot$protein_coding)) == "true" |
      tolower(trimws(annot$protein_coding)) == "protein_coding",
      ,
      drop = FALSE
    ]
  }
}

protein_genes <- unique(annot$ensembl_gene_id)

cat(
  "ℹ️ Protein-coding genes:",
  length(protein_genes),
  "\n"
)

# ---------------------------
# Filter counts
# ---------------------------
counts <- counts[
  rownames(counts) %in% protein_genes,
  ,
  drop = FALSE
]

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
id_to_symbol <- setNames(
  annot$gene_symbol,
  annot$ensembl_gene_id
)

symbols <- id_to_symbol[rownames(counts_log)]

symbols[
  is.na(symbols) | symbols == ""
] <- rownames(counts_log)[
  is.na(symbols) | symbols == ""
]

# Keep original Ensembl IDs before changing rownames
gene_ids <- rownames(counts_log)

# Make unique gene symbols for matrix rownames
unique_symbols <- make.unique(symbols)

rownames(counts_log) <- unique_symbols

# ---------------------------
# Grouping
# ---------------------------
group_col <- opt$dge_cat

if (!group_col %in% colnames(metadata)) {
  stop(
    paste0(
      "Group column not found: ",
      group_col
    )
  )
}

groups <- metadata[[group_col]]
names(groups) <- rownames(metadata)

groups <- groups[common]

counts_log <- counts_log[
  ,
  names(groups),
  drop = FALSE
]

# ---------------------------
# TOP VARIABLE GENES PER GROUP
# ---------------------------
selected_genes <- c()

for (g in unique(groups)) {

  samples <- names(groups)[groups == g]

  if (length(samples) < 2) {
    next
  }

  sub <- counts_log[
    ,
    samples,
    drop = FALSE
  ]

  gene_var <- apply(
    sub,
    1,
    var,
    na.rm = TRUE
  )

  gene_var[is.na(gene_var)] <- 0

  top_genes <- names(
    sort(
      gene_var,
      decreasing = TRUE
    )
  )[
    1:min(50, length(gene_var))
  ]

  selected_genes <- unique(
    c(selected_genes, top_genes)
  )
}

if (length(selected_genes) == 0) {
  stop("❌ No genes selected")
}

cat(
  "ℹ️ Selected genes:",
  length(selected_genes),
  "\n"
)

# ---------------------------
# EXPRESSION TABLE
# ---------------------------
# Expression values here are log2(count + 1)
# and are NOT z-scored.

expression_table <- as.data.frame(
  counts_log[
    selected_genes,
    ,
    drop = FALSE
  ],
  check.names = FALSE
)

# Add gene symbol
expression_table <- data.frame(
  gene_symbol = rownames(expression_table),
  expression_table,
  check.names = FALSE
)

# Recover Ensembl IDs
symbol_to_ensembl <- setNames(
  gene_ids,
  unique_symbols
)

expression_table <- cbind(
  ensembl_gene_id =
    unname(
      symbol_to_ensembl[
        expression_table$gene_symbol
      ]
    ),
  expression_table
)

# ---------------------------
# MATRIX FOR HEATMAP
# ---------------------------
mat <- counts_log[
  selected_genes,
  ,
  drop = FALSE
]

# Z-score each gene for heatmap
mat <- t(scale(t(mat)))

mat[is.na(mat)] <- 0

# ---------------------------
# SAMPLE ORDERING
# ---------------------------
sample_order <- c()

for (g in unique(groups)) {

  samples_g <- names(groups)[groups == g]

  if (length(samples_g) > 1) {

    sub_mat <- mat[
      ,
      samples_g,
      drop = FALSE
    ]

    hc <- hclust(
      dist(t(sub_mat))
    )

    samples_g <- samples_g[hc$order]
  }

  sample_order <- c(
    sample_order,
    samples_g
  )
}

mat <- mat[
  ,
  sample_order,
  drop = FALSE
]

# Reorder expression table to match heatmap
expression_table <- expression_table[
  ,
  c(
    "ensembl_gene_id",
    "gene_symbol",
    sample_order
  ),
  drop = FALSE
]

# ---------------------------
# SAVE EXPRESSION TABLE
# ---------------------------
out_prefix <- sub(
  "\\.[^.]+$",
  "",
  opt$out
)

expression_file <- paste0(
  out_prefix,
  "_gene_expression.csv"
)

write.csv(
  expression_table,
  expression_file,
  row.names = FALSE,
  quote = FALSE
)

cat(
  "✅ Gene-expression table saved to:\n",
  expression_file,
  "\n"
)

# ---------------------------
# SAVE Z-SCORE TABLE
# ---------------------------
zscore_table <- as.data.frame(
  mat,
  check.names = FALSE
)

zscore_table <- data.frame(
  gene_symbol = rownames(zscore_table),
  zscore_table,
  check.names = FALSE
)

zscore_file <- paste0(
  out_prefix,
  "_heatmap_zscores.csv"
)

write.csv(
  zscore_table,
  zscore_file,
  row.names = FALSE,
  quote = FALSE
)

cat(
  "✅ Heatmap z-score table saved to:\n",
  zscore_file,
  "\n"
)

# ---------------------------
# ANNOTATION
# ---------------------------
annotation_col <- data.frame(
  Group = groups[sample_order]
)

rownames(annotation_col) <- sample_order

# ---------------------------
# HEATMAP
# ---------------------------

# Dynamically adjust height based on number of genes
n_genes <- nrow(mat)

# Approximately 0.18 inch per gene
# Minimum 6 inches, maximum 40 inches
heatmap_height <- max(
  6,
  min(40, n_genes * 0.18)
)

# Dynamically adjust row font size
if (n_genes <= 50) {
  row_fontsize <- 8
} else if (n_genes <= 100) {
  row_fontsize <- 7
} else if (n_genes <= 150) {
  row_fontsize <- 6
} else {
  row_fontsize <- 5
}

# Slightly increase overall plot width
heatmap_width <- max(
  10,
  min(18, 8 + ncol(mat) * 0.15)
)

cat(
  "ℹ️ Heatmap:",
  n_genes,
  "genes ×",
  ncol(mat),
  "samples\n"
)

cat(
  "ℹ️ Plot size:",
  heatmap_width,
  "×",
  heatmap_height,
  "inches\n"
)

pheatmap(
  mat,
  annotation_col = annotation_col,
  cluster_rows = TRUE,
  cluster_cols = FALSE,

  show_rownames = TRUE,
  show_colnames = TRUE,

  # Improved readability
  fontsize_row = row_fontsize,
  fontsize_col = 7,

  # Slightly more space between rows
  cellheight = NA,
  cellwidth = NA,

  color = colorRampPalette(
    c("blue", "white", "red")
  )(100),

  filename = opt$out,

  width = heatmap_width,
  height = heatmap_height
)

