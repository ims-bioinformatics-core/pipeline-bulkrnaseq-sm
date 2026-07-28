#!/usr/bin/env Rscript

# ---------------------------
# Load packages
# ---------------------------
suppressMessages(library(optparse))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))

# ---------------------------
# Command-line options
# ---------------------------
option_list <- list(
  make_option("--pca", type="character", help="PCA table CSV"),
  make_option("--metadata", type="character", help="Metadata CSV"),
  make_option("--pca_col_num", type="character", default=NULL),
  make_option("--pca_col_cat", type="character", default=NULL),
  make_option("--pca_col_ref", type="character", default=NULL),
  make_option("--out", type="character", help="Output file path")
)

opt <- parse_args(OptionParser(option_list = option_list))

pca_file <- opt$pca
metadata_file <- opt$metadata

num_cols <- if (!is.null(opt$pca_col_num))
  trimws(strsplit(opt$pca_col_num, ",")[[1]]) else character(0)

cat_cols <- if (!is.null(opt$pca_col_cat))
  trimws(strsplit(opt$pca_col_cat, ",")[[1]]) else character(0)

refs <- if (!is.null(opt$pca_col_ref))
  trimws(strsplit(opt$pca_col_ref, ",")[[1]]) else character(0)

out_file <- opt$out
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)

# ---------------------------
# Read PCA
# ---------------------------
pca_df <- read.csv(pca_file, stringsAsFactors = FALSE, check.names = FALSE)

if (!("Barcode" %in% colnames(pca_df))) {
  stop("PCA file must contain 'Barcode'")
}

pc_cols <- grep("^PC[0-9]+$", colnames(pca_df), value = TRUE)
if (length(pc_cols) < 1) stop("No PC columns found")

# ---------------------------
# Read metadata
# ---------------------------
#lines <- readLines(metadata_file)
#start_idx <- which(grepl("Pool", lines))[1]
#if (is.na(start_idx)) stop("No 'Pool' line found in metadata file")

metadata <- read.csv(
  metadata_file,
  #skip = start_idx,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

metadata <- metadata[!is.na(metadata$Barcode) & metadata$Barcode != "", ]
metadata$Barcode <- make.unique(trimws(metadata$Barcode))
rownames(metadata) <- metadata$Barcode
metadata$Barcode <- NULL

colnames(metadata) <- gsub(" ", "_", colnames(metadata))

# ---------------------------
# Merge PCA + metadata
# ---------------------------
common <- intersect(pca_df$Barcode, rownames(metadata))
if (length(common) < 2) stop("No overlap between PCA and metadata")

pca_df <- pca_df[pca_df$Barcode %in% common, ]
metadata <- metadata[common, , drop = FALSE]

# ---------------------------
# RESULT STORAGE
# ---------------------------
results <- data.frame()

# ===========================
# NUMERIC VARIABLES
# ===========================
for (col in num_cols) {

  col_safe <- gsub(" ", "_", col)
  if (!(col_safe %in% colnames(metadata))) next

  vec <- suppressWarnings(as.numeric(metadata[[col_safe]]))
  if (all(is.na(vec))) next

  for (pc in pc_cols) {

    ok <- complete.cases(vec, pca_df[[pc]])
    if (sum(ok) < 3) next

    results <- rbind(results, data.frame(
      variable = col_safe,
      PC = pc,
      correlation = cor(vec[ok], pca_df[[pc]][ok])
    ))
  }
}

# ===========================
# CATEGORICAL VARIABLES (REFERENCE ONLY FOR COMPUTATION)
# ===========================
for (i in seq_along(cat_cols)) {

  col <- cat_cols[i]
  col_safe <- gsub(" ", "_", col)

  if (!(col_safe %in% colnames(metadata))) next

  vec_raw <- as.character(metadata[[col_safe]])

  ref <- if (length(refs) >= i) refs[i] else NA

  # reference level only for analysis
  if (!is.na(ref) && ref %in% vec_raw) {
    vec <- factor(vec_raw, levels = c(ref, setdiff(unique(vec_raw), ref)))
  } else {
    vec <- factor(vec_raw)
  }

  for (pc in pc_cols) {

    df <- data.frame(y = pca_df[[pc]], x = vec)
    df <- df[complete.cases(df), ]

    if (nrow(df) < 3 || nlevels(df$x) < 2) next

    fit <- aov(y ~ x, data = df)

    ss_total <- sum((df$y - mean(df$y))^2)
    ss_res <- sum(residuals(fit)^2)

    r2 <- 1 - (ss_res / ss_total)

    results <- rbind(results, data.frame(
      variable = col_safe,
      PC = pc,
      correlation = r2
    ))
  }
}

# ---------------------------
# ORDER PCA COLUMNS
# ---------------------------
results$PC <- factor(results$PC, levels = pc_cols)

# ---------------------------
# 🔥 FINAL FIX: ORDER BY PC1 ONLY
# ---------------------------
pc1 <- pc_cols[1]

var_order <- results %>%
  filter(PC == pc1) %>%
  arrange(desc(correlation)) %>%
  pull(variable)

# reverse for ggplot (top = strongest)
results$variable <- factor(results$variable, levels = rev(var_order))

# ---------------------------
# PLOT
# ---------------------------
p <- ggplot(results, aes(x = PC, y = variable, fill = correlation)) +

  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = round(correlation, 2)), size = 3) +

  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0
  ) +

  theme_minimal(base_size = 12) +

  labs(
    title = "Metadata association with PCA axes",
    x = "",
    y = "",
    fill = "Correlation / R²"
  ) +

  theme(
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 10),
    panel.grid = element_blank(),
    legend.position = "right"
  )

ggsave(
  out_file,
  p,
  width = max(7, length(pc_cols) * 1.2),
  height = 5,
  dpi = 150,
  bg = "white"
)

cat("✅ Correlation plot saved to:", out_file, "\n")