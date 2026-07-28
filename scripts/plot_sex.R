#!/usr/bin/env Rscript

# ---------------------------
# Load packages
# ---------------------------
suppressMessages(library(optparse))
suppressMessages(library(tidyverse))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(circlize))

# ---------------------------
# Options
# ---------------------------
option_list <- list(
  make_option("--counts", type="character"),
  make_option("--metadata", type="character"),
  make_option("--annotation", type="character"),
  make_option("--out", type="character")
)

opt <- parse_args(OptionParser(option_list = option_list))

counts_file   <- opt$counts
metadata_file <- opt$metadata
annot_file    <- opt$annotation
output_file   <- opt$out

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

# ---------------------------
# 1. Metadata
# ---------------------------
metadata <- read.csv(
  metadata_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

metadata <- metadata[
  !is.na(metadata$Barcode) &
    metadata$Barcode != "",
]

metadata$Barcode <- make.unique(trimws(metadata$Barcode))

rownames(metadata) <- metadata$Barcode
metadata$Barcode <- NULL

# ---------------------------
# 2. Counts
# ---------------------------
counts <- read.csv(
  counts_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

rownames(counts) <- counts[[1]]

counts <- counts[, -1, drop = FALSE]

colnames(counts) <- trimws(colnames(counts))

# remove Ensembl version suffix
rownames(counts) <- sub("\\..*$", "", rownames(counts))

# ---------------------------
# Sample matching
# ---------------------------
common_samples <- intersect(
  colnames(counts),
  rownames(metadata)
)

counts <- counts[, common_samples, drop = FALSE]
metadata <- metadata[common_samples, , drop = FALSE]

if (ncol(counts) < 2) {
  stop("❌ Not enough overlapping samples")
}

# ---------------------------
# 3. Read annotation
# ---------------------------
cat("ℹ️ Reading gene annotation...\n")

annot <- read.delim(
  annot_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

colnames(annot) <- trimws(colnames(annot))

cat("ℹ️ Annotation columns:\n")
print(colnames(annot))

required <- c(
  "ensembl_gene_id",
  "gene_symbol"
)

if (!all(required %in% colnames(annot))) {

  stop(
    "❌ Annotation must contain: ",
    paste(required, collapse = ", "),
    "\nFound: ",
    paste(colnames(annot), collapse = ", ")
  )
}

# ---------------------------
# Clean annotation
# ---------------------------
annot$ensembl_gene_id <- trimws(annot$ensembl_gene_id)
annot$gene_symbol     <- trimws(annot$gene_symbol)

annot <- annot[
  !is.na(annot$ensembl_gene_id) &
    annot$ensembl_gene_id != "",
  ,
  drop = FALSE
]

# fill missing symbols
missing_symbol <- is.na(annot$gene_symbol) | annot$gene_symbol == ""

annot$gene_symbol[missing_symbol] <-
  annot$ensembl_gene_id[missing_symbol]

# remove duplicate IDs
annot <- annot[
  !duplicated(annot$ensembl_gene_id),
  ,
  drop = FALSE
]

# mapping
id_to_symbol <- setNames(
  annot$gene_symbol,
  annot$ensembl_gene_id
)

# ---------------------------
# 4. Sex markers
# ---------------------------
female_markers <- c("XIST")

male_markers <- c(
  "SRY",
  "DDX3Y",
  "EIF2S3Y",
  "UTY",
  "ZFY"
)

sex_markers <- c(
  female_markers,
  male_markers
)

sex_annot <- annot[
  annot$gene_symbol %in% sex_markers,
  ,
  drop = FALSE
]

# =========================================================
# NO SEX GENES FOUND -> SAVE EMPTY FIGURE
# =========================================================
if (nrow(sex_annot) == 0) {

  cat("⚠️ No sex genes found in annotation\n")

  png(
    output_file,
    width = 900,
    height = 700,
    res = 150
  )

  plot.new()

  title("No sex genes found")

  text(
    0.5,
    0.5,
    labels = "No sex marker genes detected in annotation",
    cex = 1.2
  )

  dev.off()

  cat("✅ Empty heatmap saved to:", output_file, "\n")

  quit(save = "no", status = 0)
}

sex_ids <- unique(sex_annot$ensembl_gene_id)

# ---------------------------
# 5. Subset counts
# ---------------------------
counts_sex <- counts[
  rownames(counts) %in% sex_ids,
  ,
  drop = FALSE
]

# =========================================================
# NO SEX GENES IN COUNTS -> SAVE EMPTY FIGURE
# =========================================================
if (nrow(counts_sex) == 0) {

  cat("⚠️ No sex genes found in counts matrix\n")

  png(
    output_file,
    width = 900,
    height = 700,
    res = 150
  )

  plot.new()

  title("No sex genes found")

  text(
    0.5,
    0.5,
    labels = "No sex marker genes detected in counts matrix",
    cex = 1.2
  )

  dev.off()

  cat("✅ Empty heatmap saved to:", output_file, "\n")

  quit(save = "no", status = 0)
}

# remove zero-expression genes
counts_sex <- counts_sex[
  rowSums(counts_sex) > 0,
  ,
  drop = FALSE
]

# =========================================================
# ALL ZERO SEX GENES -> SAVE EMPTY FIGURE
# =========================================================
if (nrow(counts_sex) == 0) {

  cat("⚠️ Sex genes detected but all have zero counts\n")

  png(
    output_file,
    width = 900,
    height = 700,
    res = 150
  )

  plot.new()

  title("No expressed sex genes")

  text(
    0.5,
    0.5,
    labels = "Sex marker genes detected but expression is zero",
    cex = 1.2
  )

  dev.off()

  cat("✅ Empty heatmap saved to:", output_file, "\n")

  quit(save = "no", status = 0)
}

# ---------------------------
# 6. Convert to gene symbols
# ---------------------------
symbols <- id_to_symbol[rownames(counts_sex)]

symbols[
  is.na(symbols) |
    symbols == ""
] <- rownames(counts_sex)

rownames(counts_sex) <- symbols

# collapse duplicate symbols
counts_sex <- rowsum(
  counts_sex,
  group = rownames(counts_sex)
)

counts_sex_log <- log10(
  as.matrix(counts_sex) + 1
)

# ---------------------------
# 7. Metadata annotation
# ---------------------------
colnames(metadata)[
  colnames(metadata) == "Sex"
] <- "sex"

metadata$sex <- trimws(metadata$sex)

metadata$sex <- case_when(
  metadata$sex %in% c("Male","male","M","m") ~ "Male",
  metadata$sex %in% c("Female","female","F","f") ~ "Female",
  TRUE ~ NA_character_
)

# replace NA safely
metadata$sex[
  is.na(metadata$sex)
] <- "Unknown"

metadata$sex <- factor(metadata$sex)

sample_sex <- metadata$sex
names(sample_sex) <- rownames(metadata)

sex_colors <- c(
  Male = "#2c7fb8",
  Female = "#e41a1c",
  Unknown = "grey70"
)

sex_colors <- sex_colors[
  names(sex_colors) %in% levels(metadata$sex)
]

ha <- HeatmapAnnotation(
  Sex = sample_sex,
  col = list(Sex = sex_colors)
)

# ---------------------------
# 8. Plot heatmap
# ---------------------------
png(
  output_file,
  width = 900,
  height = 700,
  res = 150
)

Heatmap(
  counts_sex_log,

  name = "log10(count+1)",

  bottom_annotation = ha,

  cluster_rows = TRUE,
  cluster_columns = TRUE,

  show_row_names = TRUE,
  show_column_names = TRUE,

  col = colorRamp2(
    c(
      min(counts_sex_log),
      max(counts_sex_log)
    ),
    c("white", "red")
  ),

  row_names_gp = gpar(fontsize = 10),

  column_names_gp = gpar(
    fontsize = 10,
    rot = 45
  )
)

dev.off()

cat("✅ Heatmap saved to:", output_file, "\n")