#!/usr/bin/env Rscript

# ============================================================
# Gene expression plotting
#
# ONE GENE:
#   Per-sample barplot
#   - One bar per sample
#   - Bars coloured by Experimental Group
#
# MULTIPLE GENES:
#   Per-sample violin plot
#   - One violin per sample
#   - Samples coloured by Experimental Group
#   - One facet per gene
#
# Expression:
#   log2(count + 1)
#
# Examples:
#
# Single gene:
#   --genes "Fgf21"
#
# Multiple genes:
#   --genes "Fgf21,Pparg,Adipoq"
#
# Output:
#   One PNG figure
# ============================================================


# ============================================================
# LOAD PACKAGES
# ============================================================

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})


# ============================================================
# COMMAND-LINE OPTIONS
# ============================================================

option_list <- list(

  make_option(
    "--counts",
    type = "character",
    help = "Merged counts CSV"
  ),

  make_option(
    "--metadata",
    type = "character",
    help = "Metadata CSV"
  ),

  make_option(
    "--annotation",
    type = "character",
    help = "Gene annotation TSV"
  ),

  make_option(
    "--genes",
    type = "character",
    help = "Selected gene symbols, comma-separated"
  ),

  make_option(
    "--out",
    type = "character",
    help = "Output PNG"
  )
)


opt <- parse_args(
  OptionParser(
    option_list = option_list
  )
)


# ============================================================
# CHECK REQUIRED ARGUMENTS
# ============================================================

required <- c(
  "counts",
  "metadata",
  "annotation",
  "genes",
  "out"
)

missing <- required[
  sapply(
    required,
    function(x) {
      is.null(opt[[x]]) ||
        opt[[x]] == ""
    }
  )
]

if (length(missing) > 0) {

  stop(
    "Missing required argument(s): ",
    paste(
      missing,
      collapse = ", "
    )
  )
}


# ============================================================
# SELECTED GENES
# ============================================================

gene_list <- trimws(
  unlist(
    strsplit(
      opt$genes,
      ","
    )
  )
)

gene_list <- gene_list[
  gene_list != ""
]

if (length(gene_list) == 0) {
  stop(
    "No genes supplied to --genes"
  )
}

cat(
  "Selected gene(s): ",
  paste(
    gene_list,
    collapse = ", "
  ),
  "\n",
  sep = ""
)


# ============================================================
# READ COUNTS
# ============================================================

cat(
  "Reading counts...\n"
)

counts <- read.csv(
  opt$counts,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# ============================================================
# READ METADATA
# ============================================================

cat(
  "Reading metadata...\n"
)

metadata <- read.csv(
  opt$metadata,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# ============================================================
# READ ANNOTATION
# ============================================================

cat(
  "Reading annotation...\n"
)

annotation <- read.delim(
  opt$annotation,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# ============================================================
# CLEAN METADATA COLUMN NAMES
# ============================================================

metadata_names <- colnames(
  metadata
)

bad_names <- (
  is.na(metadata_names) |
    trimws(metadata_names) == ""
)

if (any(bad_names)) {

  metadata_names[bad_names] <-
    paste0(
      "Unnamed_",
      seq_len(
        sum(bad_names)
      )
    )
}

colnames(metadata) <-
  make.unique(
    metadata_names
  )


# ============================================================
# CHECK REQUIRED METADATA COLUMNS
# ============================================================

if (
  !"Barcode" %in% colnames(metadata)
) {

  stop(
    "Metadata does not contain 'Barcode'.\n",
    "Available columns:\n",
    paste(
      colnames(metadata),
      collapse = ", "
    )
  )
}


if (
  !"Experimental Group" %in%
    colnames(metadata)
) {

  stop(
    "Metadata does not contain ",
    "'Experimental Group'.\n",
    "Available columns:\n",
    paste(
      colnames(metadata),
      collapse = ", "
    )
  )
}


# ============================================================
# CLEAN METADATA VALUES
# ============================================================

metadata$Barcode <-
  trimws(
    as.character(
      metadata$Barcode
    )
  )

metadata[["Experimental Group"]] <-
  trimws(
    as.character(
      metadata[["Experimental Group"]]
    )
  )


metadata <- metadata[
  !is.na(metadata$Barcode) &
    metadata$Barcode != "",
  ,
  drop = FALSE
]

rownames(metadata) <-
  metadata$Barcode


# ============================================================
# PROCESS COUNTS
# ============================================================

if (
  ncol(counts) < 2
) {

  stop(
    "Counts file must contain ",
    "gene IDs and sample columns."
  )
}


gene_id_col <-
  colnames(counts)[1]


gene_ids <- trimws(
  as.character(
    counts[[gene_id_col]]
  )
)


# Remove Ensembl version suffix
#
# Example:
# ENSMUSG00000012345.7
#
# becomes:
# ENSMUSG00000012345

gene_ids <- sub(
  "[.].*$",
  "",
  gene_ids
)


counts <- counts[
  ,
  -1,
  drop = FALSE
]

rownames(counts) <-
  gene_ids


colnames(counts) <-
  trimws(
    colnames(counts)
  )


# Remove duplicated gene IDs
counts <- counts[
  !duplicated(
    rownames(counts)
  ),
  ,
  drop = FALSE
]


# ============================================================
# MATCH SAMPLES
# ============================================================

common_samples <- intersect(
  colnames(counts),
  rownames(metadata)
)


if (
  length(common_samples) == 0
) {

  stop(
    "No samples overlap between ",
    "counts and metadata."
  )
}


counts <- counts[
  ,
  common_samples,
  drop = FALSE
]

metadata <- metadata[
  common_samples,
  ,
  drop = FALSE
]


cat(
  "Matched samples: ",
  length(common_samples),
  "\n",
  sep = ""
)


# ============================================================
# CLEAN ANNOTATION COLUMN NAMES
# ============================================================

annotation_names <-
  colnames(annotation)

bad_annotation_names <- (
  is.na(annotation_names) |
    trimws(annotation_names) == ""
)

if (
  any(bad_annotation_names)
) {

  annotation_names[
    bad_annotation_names
  ] <-
    paste0(
      "Unnamed_",
      seq_len(
        sum(bad_annotation_names)
      )
    )
}


colnames(annotation) <-
  make.unique(
    annotation_names
  )


# ============================================================
# CHECK ANNOTATION
# ============================================================

if (
  !"ensembl_gene_id" %in%
    colnames(annotation)
) {

  stop(
    "Annotation must contain ",
    "'ensembl_gene_id'.\n",
    "Available columns:\n",
    paste(
      colnames(annotation),
      collapse = ", "
    )
  )
}


if (
  !"gene_symbol" %in%
    colnames(annotation)
) {

  stop(
    "Annotation must contain ",
    "'gene_symbol'.\n",
    "Available columns:\n",
    paste(
      colnames(annotation),
      collapse = ", "
    )
  )
}


# ============================================================
# CLEAN ANNOTATION VALUES
# ============================================================

annotation$ensembl_gene_id <-
  trimws(
    as.character(
      annotation$ensembl_gene_id
    )
  )

annotation$ensembl_gene_id <-
  sub(
    "[.].*$",
    "",
    annotation$ensembl_gene_id
  )


annotation$gene_symbol <-
  trimws(
    as.character(
      annotation$gene_symbol
    )
  )


annotation <- annotation[
  !is.na(
    annotation$ensembl_gene_id
  ) &
    annotation$ensembl_gene_id != "",
  ,
  drop = FALSE
]


annotation <- annotation[
  !is.na(
    annotation$gene_symbol
  ) &
    annotation$gene_symbol != "",
  ,
  drop = FALSE
]


annotation <- annotation[
  !duplicated(
    annotation$ensembl_gene_id
  ),
  ,
  drop = FALSE
]


# ============================================================
# FIND SELECTED GENES
# ============================================================

selected_annotation <- annotation[
  annotation$gene_symbol %in%
    gene_list,
  c(
    "ensembl_gene_id",
    "gene_symbol"
  ),
  drop = FALSE
]


if (
  nrow(selected_annotation) == 0
) {

  stop(
    "None of the selected genes ",
    "were found in the annotation.\n",
    "Requested genes: ",
    paste(
      gene_list,
      collapse = ", "
    )
  )
}


found_genes <-
  unique(
    selected_annotation$gene_symbol
  )


missing_genes <-
  setdiff(
    gene_list,
    found_genes
  )


if (
  length(missing_genes) > 0
) {

  warning(
    "These genes were not found: ",
    paste(
      missing_genes,
      collapse = ", "
    )
  )
}


cat(
  "Found gene(s): ",
  paste(
    found_genes,
    collapse = ", "
  ),
  "\n",
  sep = ""
)


# ============================================================
# KEEP GENES PRESENT IN COUNTS
# ============================================================

selected_annotation <-
  selected_annotation[
    selected_annotation$ensembl_gene_id %in%
      rownames(counts),
    ,
    drop = FALSE
  ]


if (
  nrow(selected_annotation) == 0
) {

  stop(
    "Selected genes were found in ",
    "annotation but not in counts."
  )
}


# ============================================================
# EXTRACT SELECTED COUNTS
# ============================================================

selected_ids <-
  selected_annotation$ensembl_gene_id


selected_counts <-
  counts[
    selected_ids,
    ,
    drop = FALSE
  ]


selected_counts <-
  as.data.frame(
    lapply(
      selected_counts,
      function(x) {
        as.numeric(
          as.character(x)
        )
      }
    ),
    check.names = FALSE
  )


rownames(selected_counts) <-
  selected_ids


# ============================================================
# CONVERT TO LONG FORMAT
# ============================================================

expression_data <-
  selected_counts


expression_data$ensembl_gene_id <-
  rownames(
    expression_data
  )


expression_data <-
  expression_data %>%
  pivot_longer(
    cols = -ensembl_gene_id,
    names_to = "Barcode",
    values_to = "count"
  )


expression_data <-
  expression_data %>%
  left_join(
    selected_annotation,
    by = "ensembl_gene_id"
  )


# ============================================================
# ADD EXPERIMENTAL GROUP
# ============================================================

metadata_plot <- data.frame(

  Barcode =
    as.character(
      metadata$Barcode
    ),

  Group =
    as.character(
      metadata[[
        "Experimental Group"
      ]]
    ),

  stringsAsFactors = FALSE,

  check.names = FALSE
)


metadata_plot <- metadata_plot[
  !is.na(
    metadata_plot$Group
  ) &
    metadata_plot$Group != "",
  ,
  drop = FALSE
]


expression_data <-
  expression_data %>%
  left_join(
    metadata_plot,
    by = "Barcode"
  )


# ============================================================
# REMOVE SAMPLES WITHOUT CONDITION
# ============================================================

if (
  any(
    is.na(
      expression_data$Group
    )
  )
) {

  unmatched <- unique(
    expression_data$Barcode[
      is.na(
        expression_data$Group
      )
    ]
  )

  warning(
    "Samples missing experimental-group metadata: ",
    paste(
      unmatched,
      collapse = ", "
    )
  )
}


expression_data <-
  expression_data %>%
  filter(
    !is.na(Group),
    Group != ""
  )


if (
  nrow(expression_data) == 0
) {

  stop(
    "No samples remained after ",
    "metadata matching."
  )
}


# ============================================================
# LOG TRANSFORM
# ============================================================

expression_data <-
  expression_data %>%
  mutate(
    log2_expression =
      log2(count + 1)
  )


# ============================================================
# SET SAMPLE ORDER
# ============================================================

sample_order <-
  metadata_plot$Barcode


expression_data$Barcode <-
  factor(
    expression_data$Barcode,
    levels = sample_order
  )


# ============================================================
# SET CONDITION ORDER
# ============================================================

group_order <-
  unique(
    metadata_plot$Group
  )


group_order <-
  group_order[
    !is.na(group_order) &
      group_order != ""
  ]


expression_data$Group <-
  factor(
    expression_data$Group,
    levels = group_order
  )


# ============================================================
# SET GENE ORDER
# ============================================================

gene_order <-
  gene_list[
    gene_list %in%
      unique(
        as.character(
          expression_data$gene_symbol
        )
      )
  ]


expression_data$gene_symbol <-
  factor(
    expression_data$gene_symbol,
    levels = gene_order
  )


# ============================================================
# OUTPUT DIRECTORY
# ============================================================

output_dir <-
  dirname(
    opt$out
  )


if (
  output_dir != "." &&
  !dir.exists(output_dir)
) {

  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
}


# ============================================================
# SINGLE GENE
# ============================================================

if (
  length(gene_order) == 1
) {

  cat(
    "One gene selected: ",
    "creating per-sample barplot...\n",
    sep = ""
  )


  p <- ggplot(
    expression_data,
    aes(
      x = Barcode,
      y = log2_expression,
      fill = Group
    )
  ) +

    geom_col(
      width = 0.8,
      colour = "black",
      linewidth = 0.2
    ) +

    labs(
      x = "Sample",
      y = "log2(count + 1)",
      fill = "Condition",
      title = paste0(
        as.character(
          gene_order
        ),
        " expression by sample"
      )
    ) +

    theme_bw() +

    theme(

      plot.title =
        element_text(
          hjust = 0.5,
          face = "bold"
        ),

      axis.title =
        element_text(
          face = "bold"
        ),

      axis.text.x =
        element_text(
          angle = 90,
          hjust = 1,
          vjust = 0.5,
          size = 7
        ),

      legend.title =
        element_text(
          face = "bold"
        ),

      panel.grid.minor =
        element_blank()
    )


  # Width scales with number of samples
  plot_width <-
    max(
      10,
      min(
        20,
        0.30 *
          length(
            unique(
              expression_data$Barcode
            )
          ) +
          5
      )
    )


  plot_height <- 7


# ============================================================
# MULTIPLE GENES
# ============================================================

} else {

  cat(
    "Multiple genes selected: ",
    "creating per-sample violin plot...\n",
    sep = ""
  )


  p <- ggplot(
    expression_data,
    aes(
      x = Barcode,
      y = log2_expression,
      fill = Group
    )
  ) +

    geom_violin(
      trim = FALSE,
      scale = "width",
      colour = "black",
      linewidth = 0.2
    ) +

    geom_jitter(
      width = 0.08,
      size = 0.8,
      alpha = 0.5
    ) +

    facet_wrap(
      ~ gene_symbol,
      scales = "free_y"
    ) +

    labs(
      x = "Sample",
      y = "log2(count + 1)",
      fill = "Condition",
      title = "Expression by sample"
    ) +

    theme_bw() +

    theme(

      plot.title =
        element_text(
          hjust = 0.5,
          face = "bold"
        ),

      axis.title =
        element_text(
          face = "bold"
        ),

      axis.text.x =
        element_text(
          angle = 90,
          hjust = 1,
          vjust = 0.5,
          size = 6
        ),

      strip.text =
        element_text(
          face = "bold"
        ),

      legend.title =
        element_text(
          face = "bold"
        ),

      panel.grid.minor =
        element_blank()
    )


  plot_width <-
    max(
      12,
      min(
        24,
        0.30 *
          length(
            unique(
              expression_data$Barcode
            )
          ) +
          5
      )
    )


  plot_height <-
    max(
      7,
      4 +
        3 *
        length(gene_order)
    )
}


# ============================================================
# SAVE PNG
# ============================================================

cat(
  "Saving PNG...\n"
)

ggsave(
  filename = opt$out,
  plot = p,
  width = plot_width,
  height = plot_height,
  units = "in",
  dpi = 300,
  bg = "white"
)


# ============================================================
# CHECK OUTPUT
# ============================================================

if (
  !file.exists(opt$out)
) {

  stop(
    "PNG was not created: ",
    opt$out
  )
}


file_size <-
  file.info(
    opt$out
  )$size


if (
  is.na(file_size) ||
  file_size == 0
) {

  stop(
    "PNG was created but is empty: ",
    opt$out
  )
}


cat(
  "\nPNG saved successfully:\n",
  opt$out,
  "\n"
)

cat(
  "File size: ",
  round(
    file_size / 1024
  ),
  " KB\n",
  sep = ""
)

cat(
  "Samples plotted: ",
  length(
    unique(
      expression_data$Barcode
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Genes plotted: ",
  length(gene_order),
  "\n",
  sep = ""
)
