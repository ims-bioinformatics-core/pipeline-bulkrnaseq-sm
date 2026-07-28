#!/usr/bin/env Rscript

# ---------------------------
# Load packages
# ---------------------------
suppressMessages(library(optparse))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))

# ---------------------------
# Command-line options
# ---------------------------
option_list <- list(
  make_option("--metadata", type="character"),
  make_option("--num_cols", type="character", default=NULL),
  make_option("--cat_cols", type="character", default=NULL),
  make_option("--out", type="character")
)

opt <- parse_args(OptionParser(option_list = option_list))

metadata_file <- opt$metadata
out_file <- opt$out

num_cols <- if (!is.null(opt$num_cols))
  trimws(strsplit(opt$num_cols, ",")[[1]]) else character(0)

cat_cols <- if (!is.null(opt$cat_cols))
  trimws(strsplit(opt$cat_cols, ",")[[1]]) else character(0)

dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)

# ---------------------------
# Read metadata
# ---------------------------
metadata <- read.csv(
  metadata_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

metadata <- metadata[!is.na(metadata$Barcode) & metadata$Barcode != "", ]
metadata$Barcode <- make.unique(trimws(metadata$Barcode))
rownames(metadata) <- metadata$Barcode

# clean names + values
colnames(metadata) <- gsub(" ", "_", colnames(metadata))
metadata[] <- lapply(metadata, function(x) trimws(as.character(x)))

# ---------------------------
# Prepare data
# ---------------------------
metadata$Sample <- rownames(metadata)

all_cols <- c(num_cols, cat_cols)
all_cols_safe <- gsub(" ", "_", all_cols)
all_cols_safe <- all_cols_safe[all_cols_safe %in% colnames(metadata)]

if (length(all_cols_safe) == 0) {
  stop("No valid columns found in metadata.")
}

df_long <- metadata %>%
  select(Sample, all_of(all_cols_safe)) %>%
  pivot_longer(
    cols = -Sample,
    names_to = "Variable",
    values_to = "Value"
  )

# ---------------------------
# Split numeric vs categorical
# ---------------------------
df_long_num <- df_long %>%
  filter(Variable %in% gsub(" ", "_", num_cols)) %>%
  mutate(Value = suppressWarnings(as.numeric(Value))) %>%
  filter(!is.na(Value))

df_long_cat <- df_long %>%
  filter(Variable %in% gsub(" ", "_", cat_cols))

# ---------------------------
# Build plots
# ---------------------------
plots <- list()

# CATEGORICAL → stacked proportion bars
if (nrow(df_long_cat) > 0) {
  p_cat <- ggplot(df_long_cat, aes(x = Sample, fill = Value)) +
    geom_bar(position = "fill") +
    facet_wrap(~ Variable, scales = "free_x") +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank()
    ) +
    labs(
      title = "Categorical variables across samples",
      x = "Sample",
      y = "Proportion",
      fill = "Category"
    )
  
  plots[["cat"]] <- p_cat
}

# NUMERIC → bar height
if (nrow(df_long_num) > 0) {
  p_num <- ggplot(df_long_num, aes(x = Sample, y = Value, fill = Variable)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap(~ Variable, scales = "free_y") +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    labs(
      title = "Numeric variables across samples",
      x = "Sample",
      y = "Value"
    )
  
  plots[["num"]] <- p_num
}

# ---------------------------
# Save output
# ---------------------------
if (length(plots) == 0) {
  stop("No plottable data found.")
}

# If both exist, combine vertically
if (length(plots) == 2) {
  suppressMessages(library(patchwork))
  final_plot <- plots$cat / plots$num
} else {
  final_plot <- plots[[1]]
}

ggsave(
  out_file,
  final_plot,
  width = 12,
  height = 8,
  dpi = 150,
  bg = "white"
)

cat("✅ Bar plot saved to:", out_file, "\n")