```r
#!/usr/bin/env Rscript

# ============================================================
# Metadata association with STAR uniquely mapped reads
# ============================================================

# ---------------------------
# Load packages
# ---------------------------

suppressMessages(library(optparse))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))


# ============================================================
# Command-line options
# ============================================================

option_list <- list(

    make_option(
        "--logs",
        type = "character",
        action = "store",
        nargs = "+",
        help = "STAR Log.final.out files"
    ),

    make_option(
        "--metadata",
        type = "character",
        help = "Metadata CSV"
    ),

    make_option(
        "--pca_col_num",
        type = "character",
        default = NULL,
        help = "Comma-separated numeric metadata columns"
    ),

    make_option(
        "--pca_col_cat",
        type = "character",
        default = NULL,
        help = "Comma-separated categorical metadata columns"
    ),

    make_option(
        "--pca_col_ref",
        type = "character",
        default = NULL,
        help = "Comma-separated reference levels for categorical variables"
    ),

    make_option(
        "--out",
        type = "character",
        help = "Output file path"
    )
)

opt <- parse_args(
    OptionParser(option_list = option_list)
)


# ============================================================
# Validate arguments
# ============================================================

if (is.null(opt$logs) || length(opt$logs) == 0) {
    stop("--logs is required")
}

if (is.null(opt$metadata)) {
    stop("--metadata is required")
}

if (is.null(opt$out)) {
    stop("--out is required")
}


# ============================================================
# Parse arguments
# ============================================================

# --logs is now passed as multiple arguments:
#
#   --logs sample1.Log.final.out sample2.Log.final.out ...
#
log_files <- opt$logs


# Numeric metadata columns
num_cols <- if (!is.null(opt$pca_col_num) &&
                nzchar(opt$pca_col_num)) {

    trimws(
        strsplit(
            opt$pca_col_num,
            ",",
            fixed = TRUE
        )[[1]]
    )

} else {
    character(0)
}


# Categorical metadata columns
cat_cols <- if (!is.null(opt$pca_col_cat) &&
                nzchar(opt$pca_col_cat)) {

    trimws(
        strsplit(
            opt$pca_col_cat,
            ",",
            fixed = TRUE
        )[[1]]
    )

} else {
    character(0)
}


# Reference levels
refs <- if (!is.null(opt$pca_col_ref) &&
            nzchar(opt$pca_col_ref)) {

    trimws(
        strsplit(
            opt$pca_col_ref,
            ",",
            fixed = TRUE
        )[[1]]
    )

} else {
    character(0)
}


out_file <- opt$out


# ============================================================
# Create output directory
# ============================================================

dir.create(
    dirname(out_file),
    recursive = TRUE,
    showWarnings = FALSE
)


# ============================================================
# STAR LOG PARSER
# ============================================================

get_star_value <- function(lines, key) {

    # Find the line containing the requested STAR statistic
    idx <- which(
        grepl(
            key,
            lines,
            fixed = TRUE
        )
    )

    if (length(idx) == 0) {
        return(NA_real_)
    }

    # STAR lines look like:
    #
    # Uniquely mapped reads number | 123456
    #
    # Split using "|" literally.
    #
    parts <- strsplit(
        lines[idx[1]],
        "|",
        fixed = TRUE
    )[[1]]

    if (length(parts) < 2) {
        return(NA_real_)
    }

    # The value is the last field
    value <- trimws(
        parts[length(parts)]
    )

    suppressWarnings(
        as.numeric(value)
    )
}


# ============================================================
# Parse STAR logs
# ============================================================

mapping_rows <- list()


for (log_file in log_files) {

    if (!file.exists(log_file)) {

        warning(
            "STAR log not found: ",
            log_file
        )

        next
    }


    lines <- readLines(
        log_file,
        warn = FALSE
    )


    # --------------------------------------------------------
    # Sample name
    # --------------------------------------------------------

    sample <- basename(log_file)

    sample <- sub(
        "\\.Log\\.final\\.out$",
        "",
        sample
    )


    # --------------------------------------------------------
    # Extract STAR mapping statistics
    # --------------------------------------------------------

    mapping_rows[[length(mapping_rows) + 1]] <- data.frame(

        Barcode = sample,

        total_reads =
            get_star_value(
                lines,
                "Number of input reads"
            ),

        unique_reads =
            get_star_value(
                lines,
                "Uniquely mapped reads number"
            ),

        multi_reads =
            get_star_value(
                lines,
                "Number of reads mapped to multiple loci"
            ),

        multi_toomany_reads =
            get_star_value(
                lines,
                "Number of reads mapped to too many loci"
            ),

        unmapped_reads =
            get_star_value(
                lines,
                "Number of reads unmapped"
            ),

        stringsAsFactors = FALSE
    )
}


# ============================================================
# Check parsed STAR logs
# ============================================================

if (length(mapping_rows) == 0) {

    stop(
        "No STAR log files could be parsed."
    )
}


mapping_df <- bind_rows(
    mapping_rows
)


# Remove logs where uniquely mapped reads could not be parsed

mapping_df <- mapping_df %>%
    filter(
        !is.na(unique_reads)
    )


if (nrow(mapping_df) == 0) {

    stop(
        "No valid 'Uniquely mapped reads number' values ",
        "could be extracted from the STAR logs."
    )
}


# ============================================================
# READ METADATA
# ============================================================

metadata <- read.csv(
    opt$metadata,
    stringsAsFactors = FALSE,
    check.names = FALSE
)


# ============================================================
# Check Barcode
# ============================================================

if (!("Barcode" %in% colnames(metadata))) {

    stop(
        "Metadata file must contain 'Barcode'.\n",
        "Available columns: ",
        paste(
            colnames(metadata),
            collapse = ", "
        )
    )
}


# ============================================================
# Clean metadata
# ============================================================

metadata <- metadata[
    !is.na(metadata$Barcode) &
        trimws(metadata$Barcode) != "",
    ,
    drop = FALSE
]


metadata$Barcode <- trimws(
    as.character(
        metadata$Barcode
    )
)


# Replace spaces in column names

colnames(metadata) <- gsub(
    " ",
    "_",
    colnames(metadata)
)


# ============================================================
# Match STAR logs to metadata
# ============================================================

common <- intersect(
    mapping_df$Barcode,
    metadata$Barcode
)


if (length(common) < 2) {

    stop(
        "Fewer than 2 samples overlap between STAR logs ",
        "and metadata.\n",
        "STAR samples: ",
        paste(
            mapping_df$Barcode,
            collapse = ", "
        ),
        "\nMetadata samples: ",
        paste(
            metadata$Barcode,
            collapse = ", "
        )
    )
}


mapping_df <- mapping_df[
    mapping_df$Barcode %in% common,
    ,
    drop = FALSE
]


metadata <- metadata[
    metadata$Barcode %in% common,
    ,
    drop = FALSE
]


# ------------------------------------------------------------
# Ensure identical ordering
# ------------------------------------------------------------

mapping_df <- mapping_df[
    match(
        common,
        mapping_df$Barcode
    ),
    ,
    drop = FALSE
]


metadata <- metadata[
    match(
        common,
        metadata$Barcode
    ),
    ,
    drop = FALSE
]


# ============================================================
# RESULT STORAGE
# ============================================================

results <- data.frame(

    variable = character(),

    statistic = character(),

    correlation = numeric(),

    stringsAsFactors = FALSE
)


# ============================================================
# NUMERIC VARIABLES
# ============================================================

for (col in num_cols) {

    # Metadata columns may contain spaces.
    # They were replaced with underscores above.

    col_safe <- gsub(
        " ",
        "_",
        col
    )


    if (!(col_safe %in% colnames(metadata))) {

        warning(
            "Numeric metadata column not found: ",
            col_safe
        )

        next
    }


    # Convert to numeric

    vec <- suppressWarnings(
        as.numeric(
            metadata[[col_safe]]
        )
    )


    if (all(is.na(vec))) {

        warning(
            "Numeric column contains no usable values: ",
            col_safe
        )

        next
    }


    # Response variable:
    # STAR uniquely mapped reads

    y <- mapping_df$unique_reads


    # Complete observations

    ok <- complete.cases(
        vec,
        y
    )


    if (sum(ok) < 3) {

        warning(
            "Not enough complete observations for: ",
            col_safe
        )

        next
    }


    # Pearson correlation

    r <- suppressWarnings(
        cor(
            vec[ok],
            y[ok],
            method = "pearson"
        )
    )


    results <- rbind(

        results,

        data.frame(

            variable = col_safe,

            statistic = "Pearson r",

            correlation = r,

            stringsAsFactors = FALSE
        )
    )
}


# ============================================================
# CATEGORICAL VARIABLES
# ============================================================

for (i in seq_along(cat_cols)) {

    col <- cat_cols[i]


    col_safe <- gsub(
        " ",
        "_",
        col
    )


    if (!(col_safe %in% colnames(metadata))) {

        warning(
            "Categorical metadata column not found: ",
            col_safe
        )

        next
    }


    # --------------------------------------------------------
    # Convert to character
    # --------------------------------------------------------

    vec_raw <- as.character(
        metadata[[col_safe]]
    )


    # Treat empty values as NA

    vec_raw[
        is.na(vec_raw) |
            trimws(vec_raw) == ""
    ] <- NA_character_


    # --------------------------------------------------------
    # Reference level
    # --------------------------------------------------------

    ref <- if (
        length(refs) >= i
    ) {
        refs[i]
    } else {
        NA_character_
    }


    # --------------------------------------------------------
    # Set factor levels
    # --------------------------------------------------------

    if (
        !is.na(ref) &&
        ref %in% vec_raw
    ) {

        other_levels <- setdiff(
            unique(
                vec_raw[
                    !is.na(vec_raw)
                ]
            ),
            ref
        )


        vec <- factor(
            vec_raw,
            levels = c(
                ref,
                other_levels
            )
        )

    } else {

        vec <- factor(
            vec_raw
        )
    }


    # --------------------------------------------------------
    # Build analysis data frame
    # --------------------------------------------------------

    df <- data.frame(

        y = mapping_df$unique_reads,

        x = vec
    )


    df <- df[
        complete.cases(df),
        ,
        drop = FALSE
    ]


    if (nrow(df) < 3) {

        warning(
            "Not enough observations for: ",
            col_safe
        )

        next
    }


    if (nlevels(df$x) < 2) {

        warning(
            "Categorical variable has fewer than 2 levels: ",
            col_safe
        )

        next
    }


    # --------------------------------------------------------
    # ANOVA
    # --------------------------------------------------------

    fit <- aov(
        y ~ x,
        data = df
    )


    # Total sum of squares

    ss_total <- sum(
        (
            df$y -
                mean(df$y)
        )^2
    )


    # Residual sum of squares

    ss_res <- sum(
        residuals(fit)^2
    )


    if (ss_total == 0) {

        warning(
            "No variance in uniquely mapped reads for: ",
            col_safe
        )

        next
    }


    # Eta-squared / ANOVA R²

    r2 <- 1 - (
        ss_res /
            ss_total
    )


    results <- rbind(

        results,

        data.frame(

            variable = col_safe,

            statistic = "ANOVA R²",

            correlation = r2,

            stringsAsFactors = FALSE
        )
    )
}


# ============================================================
# CHECK RESULTS
# ============================================================

if (nrow(results) == 0) {

    stop(
        "No valid metadata associations could be calculated."
    )
}


# ============================================================
# ORDER VARIABLES
# ============================================================

# Strongest association at the top.
#
# Numeric:
#   absolute Pearson correlation
#
# Categorical:
#   ANOVA R²
#
# ============================================================

var_order <- results %>%

    mutate(
        order_value = ifelse(
            statistic == "Pearson r",
            abs(correlation),
            correlation
        )
    ) %>%

    group_by(variable) %>%

    summarise(
        order_value = max(
            order_value,
            na.rm = TRUE
        ),
        .groups = "drop"
    ) %>%

    arrange(
        desc(order_value)
    ) %>%

    pull(variable)


results$variable <- factor(
    results$variable,
    levels = rev(var_order)
)


# ============================================================
# PLOT
# ============================================================

p <- ggplot(
    results,
    aes(
        x = statistic,
        y = variable,
        fill = correlation
    )
) +

    geom_tile(
        color = "white",
        linewidth = 0.5
    ) +

    geom_text(
        aes(
            label = round(
                correlation,
                2
            )
        ),
        size = 3
    ) +

    scale_fill_gradient2(
        low = "blue",
        mid = "white",
        high = "red",
        midpoint = 0
    ) +

    theme_minimal(
        base_size = 12
    ) +

    labs(
        title = "Metadata association with uniquely mapped reads",
        x = "",
        y = "",
        fill = "Correlation / R²"
    ) +

    theme(
        axis.text.x = element_text(
            size = 11
        ),

        axis.text.y = element_text(
            size = 10
        ),

        panel.grid = element_blank(),

        legend.position = "right"
    )


# ============================================================
# SAVE
# ============================================================

ggsave(
    out_file,
    p,
    width = 7,
    height = max(
        5,
        length(var_order) * 0.35
    ),
    dpi = 300,
    bg = "white"
)


# ============================================================
# DONE
# ============================================================

cat(
    "Metadata association plot saved to: ",
    out_file,
    "\n"
)
```
