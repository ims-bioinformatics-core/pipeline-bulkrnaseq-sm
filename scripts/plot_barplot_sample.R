#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(optparse)
    library(ggplot2)
    library(dplyr)
    library(tidyr)
})

# ---------------------------------------------------------
# Arguments
# ---------------------------------------------------------

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
        "--group",
        type = "character",
        default = "Experimental Group",
        help = "Metadata column used for grouping"
    ),
    make_option(
        "--out",
        type = "character",
        help = "Output PNG"
    )
)

opt <- parse_args(OptionParser(option_list = option_list))

cat("\n")
cat("=============================================\n")
cat(" Mitochondrial / Ribosomal sample QC\n")
cat("=============================================\n")
cat("Counts:      ", opt$counts, "\n", sep = "")
cat("Metadata:    ", opt$metadata, "\n", sep = "")
cat("Annotation:  ", opt$annotation, "\n", sep = "")
cat("Output:      ", opt$out, "\n", sep = "")
cat("Group:       ", opt$group, "\n", sep = "")
cat("\n")


# ---------------------------------------------------------
# Read counts
# ---------------------------------------------------------

cat("Reading counts...\n")

counts <- read.csv(
    opt$counts,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

cat(
    "Counts dimensions: ",
    nrow(counts),
    " genes x ",
    ncol(counts),
    " columns\n",
    sep = ""
)

cat("First count columns:\n")
print(head(colnames(counts), 10))

gene_id_col <- colnames(counts)[1]

cat("\nGene ID column: ", gene_id_col, "\n", sep = "")

sample_names <- colnames(counts)[-1]

cat("Number of samples: ", length(sample_names), "\n", sep = "")

if (any(sample_names == "")) {
    stop("Counts contain empty sample column names.")
}

# Strip Ensembl version numbers
counts[[gene_id_col]] <- sub(
    "\\..*$",
    "",
    as.character(counts[[gene_id_col]])
)

# Convert count columns to numeric
counts[sample_names] <- lapply(
    counts[sample_names],
    function(x) as.numeric(as.character(x))
)


# ---------------------------------------------------------
# Read metadata
# ---------------------------------------------------------

cat("\nReading metadata...\n")

metadata <- read.csv(
    opt$metadata,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

cat(
    "Metadata dimensions: ",
    nrow(metadata),
    " rows x ",
    ncol(metadata),
    " columns\n",
    sep = ""
)

cat("Metadata columns:\n")
print(colnames(metadata))

# Empty metadata column names are allowed.
# They are irrelevant because we only use Barcode and Group.
if (any(colnames(metadata) == "")) {
    cat(
        "Note: metadata contains ",
        sum(colnames(metadata) == ""),
        " empty column name(s); these will be ignored.\n",
        sep = ""
    )
}

if (!"Barcode" %in% colnames(metadata)) {
    stop("Metadata must contain a 'Barcode' column.")
}

if (!opt$group %in% colnames(metadata)) {
    stop(
        "Metadata does not contain requested grouping column: ",
        opt$group
    )
}

metadata_use <- data.frame(
    Barcode = as.character(metadata[["Barcode"]]),
    Group = as.character(metadata[[opt$group]]),
    stringsAsFactors = FALSE
)

metadata_use$Barcode <- trimws(metadata_use$Barcode)

# Remove completely empty barcode rows
metadata_use <- metadata_use[
    !is.na(metadata_use$Barcode) &
        metadata_use$Barcode != "",
    ,
    drop = FALSE
]

cat(
    "Usable metadata rows: ",
    nrow(metadata_use),
    "\n",
    sep = ""
)


# ---------------------------------------------------------
# Match samples
# ---------------------------------------------------------

common_samples <- intersect(
    sample_names,
    metadata_use$Barcode
)

cat(
    "Samples matched between counts and metadata: ",
    length(common_samples),
    "\n",
    sep = ""
)

if (length(common_samples) == 0) {
    stop(
        "No samples could be matched between counts and metadata."
    )
}

if (length(common_samples) < length(sample_names)) {
    missing_metadata <- setdiff(sample_names, common_samples)

    cat("\nSamples missing from metadata:\n")
    print(missing_metadata)
}

# Keep only matched samples
counts <- counts[
    ,
    c(gene_id_col, common_samples),
    drop = FALSE
]

metadata_use <- metadata_use[
    match(common_samples, metadata_use$Barcode),
    ,
    drop = FALSE
]

rownames(metadata_use) <- NULL

# Confirm ordering
if (!identical(
    colnames(counts)[-1],
    metadata_use$Barcode
)) {
    stop("Sample ordering mismatch between counts and metadata.")
}


# ---------------------------------------------------------
# Read annotation
# ---------------------------------------------------------

cat("\nReading annotation...\n")

annotation <- read.delim(
    opt$annotation,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

required_annotation_columns <- c(
    "ensembl_gene_id",
    "gene_symbol"
)

missing_annotation_columns <- setdiff(
    required_annotation_columns,
    colnames(annotation)
)

if (length(missing_annotation_columns) > 0) {
    stop(
        "Annotation is missing required column(s): ",
        paste(missing_annotation_columns, collapse = ", ")
    )
}

annotation <- annotation[
    ,
    c("ensembl_gene_id", "gene_symbol"),
    drop = FALSE
]

annotation$ensembl_gene_id <- sub(
    "\\..*$",
    "",
    as.character(annotation$ensembl_gene_id)
)

annotation$gene_symbol <- as.character(annotation$gene_symbol)

# Remove duplicate Ensembl IDs
annotation <- annotation[
    !duplicated(annotation$ensembl_gene_id),
    ,
    drop = FALSE
]


# ---------------------------------------------------------
# Map Ensembl IDs to gene symbols
# ---------------------------------------------------------

gene_symbols <- annotation$gene_symbol[
    match(
        counts[[gene_id_col]],
        annotation$ensembl_gene_id
    )
]

cat("\nAnnotation matching:\n")
cat(
    "Genes in counts:       ",
    nrow(counts),
    "\n",
    sep = ""
)
cat(
    "Genes with annotation: ",
    sum(!is.na(gene_symbols) & gene_symbols != ""),
    "\n",
    sep = ""
)
cat(
    "Genes without symbol:  ",
    sum(is.na(gene_symbols) | gene_symbols == ""),
    "\n",
    sep = ""
)


# ---------------------------------------------------------
# Identify mitochondrial and ribosomal genes
# ---------------------------------------------------------

is_mito <- !is.na(gene_symbols) &
    grepl(
        "^mt-",
        gene_symbols,
        ignore.case = FALSE
    )

is_ribo <- !is.na(gene_symbols) &
    grepl(
        "^(Rps|Rpl)",
        gene_symbols
    )

cat("\nGene classification:\n")
cat(
    "Mitochondrial genes: ",
    sum(is_mito),
    "\n",
    sep = ""
)
cat(
    "Ribosomal genes:     ",
    sum(is_ribo),
    "\n",
    sep = ""
)

cat("\nExample mitochondrial genes:\n")
print(
    head(
        unique(gene_symbols[is_mito]),
        10
    )
)

cat("\nExample ribosomal genes:\n")
print(
    head(
        unique(gene_symbols[is_ribo]),
        10
    )
)

if (sum(is_mito) == 0) {
    stop("No mitochondrial genes were identified.")
}

if (sum(is_ribo) == 0) {
    stop("No ribosomal genes were identified.")
}


# ---------------------------------------------------------
# Calculate QC percentages
# ---------------------------------------------------------

count_matrix <- as.matrix(
    counts[, common_samples, drop = FALSE]
)

storage.mode(count_matrix) <- "numeric"

total_counts <- colSums(
    count_matrix,
    na.rm = TRUE
)

mito_counts <- colSums(
    count_matrix[is_mito, , drop = FALSE],
    na.rm = TRUE
)

ribo_counts <- colSums(
    count_matrix[is_ribo, , drop = FALSE],
    na.rm = TRUE
)

qc <- data.frame(
    Sample = common_samples,
    Group = metadata_use$Group,
    Total = total_counts,
    Mito = mito_counts,
    Ribo = ribo_counts,
    stringsAsFactors = FALSE
)

qc$MitoPercent <- 100 * qc$Mito / qc$Total
qc$RiboPercent <- 100 * qc$Ribo / qc$Total

qc_long <- qc %>%
    select(
        Sample,
        Group,
        MitoPercent,
        RiboPercent
    ) %>%
    pivot_longer(
        cols = c(MitoPercent, RiboPercent),
        names_to = "Category",
        values_to = "Percent"
    )

qc_long$Category <- recode(
    qc_long$Category,
    MitoPercent = "Mitochondrial",
    RiboPercent = "Ribosomal"
)

# Preserve sample order
qc_long$Sample <- factor(
    qc_long$Sample,
    levels = common_samples
)

# Preserve group information
qc_long$Group <- factor(
    qc_long$Group
)


# ---------------------------------------------------------
# Write QC table
# ---------------------------------------------------------

qc_out <- sub(
    "\\.[Pp][Nn][Gg]$",
    "_mito_ribo.csv",
    opt$out
)

write.csv(
    qc,
    qc_out,
    row.names = FALSE
)

cat("\nQC table written to:\n")
cat(qc_out, "\n")

# ---------------------------------------------------------
# Plot
# ---------------------------------------------------------

# Horizontal bars make sample names readable when there are
# many samples.

# Reverse sample order so the first sample appears at the top
qc_long$Sample <- factor(
    qc_long$Sample,
    levels = rev(common_samples)
)

p <- ggplot(
    qc_long,
    aes(
        x = Percent,
        y = Sample,
        fill = Group
    )
) +
    geom_col(
        width = 0.75
    ) +
    facet_wrap(
        ~ Category,
        ncol = 1,
        scales = "free_x"
    ) +
    labs(
        title = "Mitochondrial and Ribosomal RNA Levels",
        x = "Percentage of total counts",
        y = NULL,
        fill = opt$group
    ) +
    theme_bw(base_size = 11) +
    theme(
        axis.text.y = element_text(
            size = 7
        ),
        axis.text.x = element_text(
            size = 9
        ),
        strip.text = element_text(
            face = "bold",
            size = 11
        ),
        legend.position = "right",
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(
            face = "bold",
            hjust = 0.5
        )
    )

# Height scales with the number of samples so that sample
# names remain readable.
plot_height <- max(
    10,
    length(common_samples) * 0.18
)

ggsave(
    filename = opt$out,
    plot = p,
    width = 14,
    height = plot_height,
    units = "in",
    dpi = 300,
    limitsize = FALSE
)

cat("\nPlot written to:\n")
cat(opt$out, "\n")

cat("\nDone.\n")

