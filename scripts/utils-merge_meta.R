#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop("Usage: utils-generate_filepaths.R <contents.csv> <sample_info.csv> <output.csv>")
}

contents_file <- args[1]
meta_file     <- args[2]
output_file   <- args[3]

# ---------------------------
# Read input files
# ---------------------------
contents <- read_csv(contents_file, show_col_types = FALSE)
meta     <- read_csv(meta_file, show_col_types = FALSE)

#names(contents) <- gsub("\\.", "_", names(contents))
#names(meta) <- gsub("\\.", "_", names(meta))
# Check required columns
if (!"Sample name" %in% colnames(contents)) {
  stop("contents.csv must contain 'Sample name'")
}

if (!"Sample name" %in% colnames(meta)) {
  stop("sample_info.csv must contain 'Sample name'")
}

# ---------------------------
# Match metadata to contents
# ---------------------------
contents <- contents %>%
  mutate(
    `Sample name` = `Sample name` %>%
      str_remove("^GBC\\d+_") %>%   # remove GBC123_
      str_remove("_[0-9]+$")        # remove any trailing _number
  )
merged <- contents %>%
  left_join(meta, by = "Sample name")

# ---------------------------
# Save output
# ---------------------------
write_csv(merged, output_file)

cat("✅ Metadata merged and saved to:", output_file, "\n")
cat("Rows:", nrow(merged), "\n")