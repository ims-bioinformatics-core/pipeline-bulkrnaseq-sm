#!/usr/bin/env Rscript

# ---------------------------
# Load packages
# ---------------------------
suppressMessages(library(optparse))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(vcd))

# ---------------------------
# Command-line options
# ---------------------------
option_list <- list(
  make_option("--pca", type="character", help="PCA table CSV"),
  make_option("--metadata", type = "character"),
  make_option("--pca_col_num", type = "character", default = NULL),
  make_option("--pca_col_cat", type = "character", default = NULL),
  make_option("--pca_col_ref", type = "character", default = NULL),
  make_option("--out", type = "character")
)

opt <- parse_args(OptionParser(option_list = option_list))

metadata_file <- opt$metadata
pca_file <- opt$pca
out_file <- opt$out

num_cols <- if (!is.null(opt$pca_col_num))
  trimws(strsplit(opt$pca_col_num, ",")[[1]]) else character(0)

cat_cols <- if (!is.null(opt$pca_col_cat))
  trimws(strsplit(opt$pca_col_cat, ",")[[1]]) else character(0)

dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)

# ---------------------------
# Read PCA
# ---------------------------
pca_df <- read.csv(pca_file, stringsAsFactors = FALSE, check.names = FALSE)

if (!("Barcode" %in% colnames(pca_df))) {
  stop("PCA file must contain 'Barcode'")
}

pca_df$Barcode <- trimws(pca_df$Barcode)

pc_cols <- grep("^PC[0-9]+$", colnames(pca_df), value = TRUE)
if (length(pc_cols) < 1) stop("No PC columns found")

# ---------------------------
# Read metadata
# ---------------------------
metadata <- read.csv(
  metadata_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (!("Barcode" %in% colnames(metadata))) {
  stop("Metadata must contain 'Barcode'")
}

metadata$Barcode <- trimws(metadata$Barcode)

# ---------------------------
# Subset metadata to PCA samples ONLY
# ---------------------------
common <- intersect(metadata$Barcode, pca_df$Barcode)

if (length(common) < 2) stop("No overlap between PCA and metadata")

metadata <- metadata[metadata$Barcode %in% common, , drop = FALSE]
pca_df <- pca_df[pca_df$Barcode %in% common, , drop = FALSE]

# Optional: enforce identical ordering
metadata <- metadata[order(metadata$Barcode), , drop = FALSE]
pca_df <- pca_df[order(pca_df$Barcode), , drop = FALSE]

# ---------------------------
# Clean metadata
# ---------------------------
colnames(metadata) <- gsub(" ", "_", colnames(metadata))

metadata[] <- lapply(metadata, function(x) {
  trimws(as.character(x))
})

# ---------------------------
# Variable lists
# ---------------------------
num_cols_safe <- gsub(" ", "_", num_cols)
cat_cols_safe <- gsub(" ", "_", cat_cols)

vars <- unique(c(cat_cols_safe, num_cols_safe))
vars <- vars[vars %in% colnames(metadata)]

# ---------------------------
# Result storage
# ---------------------------
results <- data.frame()

# =========================================================
# NUMERIC ↔ NUMERIC (Pearson correlation)
# =========================================================
for (col1 in num_cols_safe) {

  if (!(col1 %in% colnames(metadata))) next

  x <- suppressWarnings(as.numeric(metadata[[col1]]))
  if (all(is.na(x))) next

  for (col2 in num_cols_safe) {

    if (!(col2 %in% colnames(metadata))) next

    if (col1 == col2) {
      results <- rbind(results, data.frame(
        var1 = col1,
        var2 = col2,
        value = 1,
        type = "numeric"
      ))
      next
    }

    y <- suppressWarnings(as.numeric(metadata[[col2]]))
    if (all(is.na(y))) next

    ok <- complete.cases(x, y)
    if (sum(ok) < 3) next

    r <- cor(x[ok], y[ok], method = "pearson")

    results <- rbind(results, data.frame(
      var1 = col1,
      var2 = col2,
      value = r,
      type = "numeric"
    ))
  }
}

# =========================================================
# CATEGORICAL ↔ CATEGORICAL (Cramer's V)
# =========================================================
for (col1 in cat_cols_safe) {

  if (!(col1 %in% colnames(metadata))) next

  x_raw <- metadata[[col1]]

  for (col2 in cat_cols_safe) {

    if (!(col2 %in% colnames(metadata))) next

    if (col1 == col2) {
      results <- rbind(results, data.frame(
        var1 = col1,
        var2 = col2,
        value = 1,
        type = "categorical"
      ))
      next
    }

    y_raw <- metadata[[col2]]

    df <- data.frame(x = x_raw, y = y_raw)
    df <- df[complete.cases(df), ]

    if (nrow(df) < 3) next

    tbl <- table(df$x, df$y)
    if (min(dim(tbl)) < 2) next

    cramers_v <- suppressWarnings(assocstats(tbl)$cramer)

    results <- rbind(results, data.frame(
      var1 = col1,
      var2 = col2,
      value = cramers_v,
      type = "categorical"
    ))
  }
}

# =========================================================
# CATEGORICAL ↔ NUMERIC (ANOVA R²)
# =========================================================
for (cat in cat_cols_safe) {

  if (!(cat %in% colnames(metadata))) next

  x_raw <- metadata[[cat]]

  for (num in num_cols_safe) {

    if (!(num %in% colnames(metadata))) next

    y <- suppressWarnings(as.numeric(metadata[[num]]))
    if (all(is.na(y))) next

    df <- data.frame(cat = x_raw, num = y)
    df <- df[complete.cases(df), ]

    if (nrow(df) < 3) next
    if (nlevels(factor(df$cat)) < 2) next

    fit <- aov(num ~ as.factor(cat), data = df)

    ss_total <- sum((df$num - mean(df$num))^2)
    ss_res <- sum(residuals(fit)^2)

    r2 <- 1 - (ss_res / ss_total)

    results <- rbind(results, data.frame(
      var1 = cat,
      var2 = num,
      value = r2,
      type = "cat_num"
    ))

    results <- rbind(results, data.frame(
      var1 = num,
      var2 = cat,
      value = r2,
      type = "num_cat"
    ))
  }
}

# ---------------------------
# Force full matrix
# ---------------------------
full_grid <- expand.grid(
  var1 = vars,
  var2 = vars,
  stringsAsFactors = FALSE
)

results_full <- full_grid %>%
  left_join(results %>% select(var1, var2, value),
            by = c("var1", "var2")) %>%
  mutate(value = ifelse(var1 == var2, 1, value))

# ---------------------------
# Axis ordering
# ---------------------------
results_full$var1 <- factor(results_full$var1, levels = vars)
results_full$var2 <- factor(results_full$var2, levels = vars)

# ---------------------------
# Labels
# ---------------------------
results_full$label <- ifelse(
  is.na(results_full$value),
  "",
  sprintf("%.2f", results_full$value)
)

# ---------------------------
# Plot
# ---------------------------
p <- ggplot(results_full, aes(x = var1, y = var2, fill = value)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = label), size = 3) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    na.value = "grey90",
    limits = c(-1, 1)
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  labs(
    title = "Metadata Variable Correlations (PCA-subset samples)",
    x = "",
    y = "",
    fill = "Correlation / Association"
  )

# ---------------------------
# Save
# ---------------------------
ggsave(out_file, p, width = 8, height = 6, dpi = 150, bg = "white")

cat("✅ Metadata correlation plot saved to:", out_file, "\n")