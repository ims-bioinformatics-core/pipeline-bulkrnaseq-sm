#!/usr/bin/env Rscript

# ---------------------------
# Load packages
# ---------------------------
suppressMessages(library(optparse))
suppressMessages(library(DESeq2))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))

# ---------------------------
# Command-line options
# ---------------------------
option_list <- list(
  make_option("--counts", type = "character"),
  make_option("--metadata", type = "character"),
  make_option("--annotation", type = "character"),
  make_option("--out", type = "character"),
  make_option("--dge_cat", type = "character"),
  make_option("--dge_cat_reference", type = "character"),
  make_option("--filt", type = "logical", default = TRUE),
  make_option("--lfc_threshold", type = "numeric", default = 1),
  make_option("--padj_threshold", type = "numeric", default = 0.05)
)

opt <- parse_args(OptionParser(option_list = option_list))

counts_file   <- opt$counts
metadata_file <- opt$metadata
annot_file    <- opt$annotation
output_file   <- opt$out
dge_cat       <- opt$dge_cat
dge_cat_ref   <- opt$dge_cat_reference
apply_filter  <- opt$filt
lfc_thresh    <- opt$lfc_threshold
padj_thresh   <- opt$padj_threshold

# ---------------------------
# Validate required args
# ---------------------------
required_args <- c("counts","metadata","annotation","out","dge_cat","dge_cat_reference")

for (arg in required_args) {
  if (is.null(opt[[arg]]) || opt[[arg]] == "") {
    stop(sprintf("❌ --%s is required.", arg))
  }
}

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

# ---------------------------
# Volcano function (SAFE FILE NAMES)
# ---------------------------
create_volcano_plot <- function(df, contrast_name, outdir, lfc_thresh, padj_thresh) {

  safe_name <- gsub("[^A-Za-z0-9_\\-]", "_", contrast_name)

  df <- df %>%
    mutate(
      neg_log10_padj = -log10(padj),
      gene_label = ifelse(is.na(gene_symbol) | gene_symbol == "",
                          ensembl_gene_id,
                          gene_symbol),
      direction = case_when(
        significant & log2FoldChange > 0 ~ "Up",
        significant & log2FoldChange < 0 ~ "Down",
        TRUE ~ "NS"
      )
    )

  top_up <- df %>%
    filter(significant, log2FoldChange > 0) %>%
    arrange(padj) %>%
    head(20)

  top_down <- df %>%
    filter(significant, log2FoldChange < 0) %>%
    arrange(padj) %>%
    head(20)

  label_df <- bind_rows(top_up, top_down)

  p <- ggplot(df, aes(x = log2FoldChange, y = neg_log10_padj)) +
    geom_point(aes(color = direction), alpha = 0.6, size = 1.2) +
    scale_color_manual(values = c("Up"="red", "Down"="blue", "NS"="grey70")) +
    geom_vline(xintercept = c(-lfc_thresh, lfc_thresh), linetype = "dashed") +
    geom_hline(yintercept = -log10(padj_thresh), linetype = "dashed") +
    geom_text(data = label_df,
              aes(label = gene_label),
              size = 3,
              check_overlap = TRUE) +
    labs(
      title = paste("Volcano:", contrast_name),
      x = "Log2 Fold Change",
      y = "-Log10 Adjusted P-value"
    ) +
    theme_bw()

  out_file <- file.path(outdir, paste0("volcano_", safe_name, ".png"))

  ggsave(out_file, p, width = 8, height = 6, dpi = 150)

  cat("✅ Volcano saved:", out_file, "\n")
}

# ---------------------------
# 1. Metadata
# ---------------------------
metadata <- read.csv(metadata_file, stringsAsFactors = FALSE, check.names = FALSE)

metadata <- metadata[!is.na(metadata$Barcode) & trimws(metadata$Barcode) != "", ]
metadata$Barcode <- make.unique(trimws(metadata$Barcode))
rownames(metadata) <- metadata$Barcode
metadata$Barcode <- NULL

colnames(metadata) <- make.names(colnames(metadata))
dge_cat <- make.names(dge_cat)

metadata[[dge_cat]] <- factor(metadata[[dge_cat]])
metadata[[dge_cat]] <- relevel(metadata[[dge_cat]], ref = dge_cat_ref)

# ---------------------------
# 2. Counts
# ---------------------------
counts <- read.csv(counts_file, stringsAsFactors = FALSE, check.names = FALSE)

rownames(counts) <- counts[[1]]
counts <- counts[, -1, drop = FALSE]
rownames(counts) <- sub("\\..*$", "", rownames(counts))

common_samples <- intersect(colnames(counts), rownames(metadata))
counts <- counts[, common_samples, drop = FALSE]
metadata <- metadata[common_samples, , drop = FALSE]

# ---------------------------
# 3. Annotation
# ---------------------------
annot <- read.delim(annot_file, header = TRUE, stringsAsFactors = FALSE)

annot <- annot %>%
  mutate(
    ensembl_gene_id = trimws(ensembl_gene_id),
    gene_symbol = trimws(gene_symbol)
  ) %>%
  filter(ensembl_gene_id != "")

protein_gene_ids <- unique(annot$ensembl_gene_id)

counts <- counts[rownames(counts) %in% protein_gene_ids, ]

# ---------------------------
# 4. DESeq2
# ---------------------------
dds <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(counts)),
  colData = metadata,
  design = as.formula(paste("~", dge_cat))
)

if (apply_filter) {
  dds <- dds[rowSums(counts(dds) >= 10) >= 3, ]
}

dds <- DESeq(dds)

result_levels <- levels(metadata[[dge_cat]])
contrast_levels <- result_levels[result_levels != dge_cat_ref]

all_results <- list()

# ---------------------------
# 5. Results + Volcano
# ---------------------------
for (lvl in contrast_levels) {

  contrast_name <- paste0(lvl, "_vs_", dge_cat_ref)

  res <- results(dds,
                 contrast = c(dge_cat, lvl, dge_cat_ref),
                 alpha = padj_thresh)

  res_df <- as.data.frame(res)
  res_df$ensembl_gene_id <- rownames(res_df)

  res_df <- res_df %>%
    left_join(annot %>% select(ensembl_gene_id, gene_symbol),
              by = "ensembl_gene_id")

  res_df$significant <- !is.na(res_df$padj) &
    res_df$padj < padj_thresh &
    abs(res_df$log2FoldChange) >= lfc_thresh

  res_df$contrast <- contrast_name

  all_results[[contrast_name]] <- res_df

  create_volcano_plot(
    res_df,
    contrast_name,
    dirname(output_file),
    lfc_thresh,
    padj_thresh
  )

  cat("✅ Completed:", contrast_name, "\n")
}

# ---------------------------
# 6. Save DGE results
# ---------------------------
combined_results <- bind_rows(all_results)
write.csv(combined_results, output_file, row.names = FALSE)

cat("✅ DGE saved:", output_file, "\n")

# ---------------------------
# 7. VST
# ---------------------------
vst_obj <- vst(dds, blind = FALSE)
vst_df <- as.data.frame(assay(vst_obj))
vst_df$ensembl_gene_id <- rownames(vst_df)

vst_out <- file.path(dirname(output_file), "vst_normalised_counts.csv")
write.csv(vst_df, vst_out, row.names = FALSE)

cat("✅ VST saved:", vst_out, "\n")