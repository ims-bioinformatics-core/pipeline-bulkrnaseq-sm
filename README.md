# pipeline-bulkrnaseq-sm

A reproducible bulk RNA-seq pipeline implemented in Snakemake, designed for HPC execution under SLURM.

It runs end-to-end from paired-end FASTQ files to QC reports, STAR alignment, a merged gene-count matrix,
sample- and expression-level QC plots, PCA, DESeq2 differential expression (one or more named
comparisons) and GSEA pathway analysis.

---

## Table of contents

- [Workflow](#workflow)
- [Requirements](#requirements)
- [Installation](#installation)
- [Inputs](#inputs)
- [Configuration reference](#configuration-reference)
- [Running the pipeline](#running-the-pipeline)
- [Outputs](#outputs)
- [Helper scripts](#helper-scripts)
- [Known issues](#known-issues)
- [Author](#author)

---

## Workflow

Rules defined in `Snakefile`, in dependency order:

| # | Rule | Tool | Conda env | Produces |
|---|------|------|-----------|----------|
| 1 | `fastqc` | FastQC | `pipeline-bulkrnaseq-sm_env` | Per-sample QC of raw reads |
| 2 | `multiqc` | MultiQC | `pipeline-bulkrnaseq-sm_env` | Aggregated raw-read QC report |
| 3 | `trim_reads` | Cutadapt | `pipeline-bulkrnaseq-sm_env` | Trimmed FASTQ pairs |
| 4 | `fastqc_trimmed` | FastQC | `pipeline-bulkrnaseq-sm_env` | Per-sample QC of trimmed reads |
| 5 | `multiqc_trimmed` | MultiQC | `pipeline-bulkrnaseq-sm_env` | Aggregated trimmed-read QC report |
| 6 | `run_star` | STAR | `pipeline-bulkrnaseq-sm_env` | Coordinate-sorted BAM + `ReadsPerGene` counts |
| 7 | `multiqc_mapped` | MultiQC | `pipeline-bulkrnaseq-sm_env` | Aggregated alignment QC report |
| 8 | `mapping_check` | Python | `pipeline-bulkrnaseq-sm_env` | Stacked barplot of mapping rates, grouped by `mapping_check_cat` |
| 9 | `merge_counts` | Python | `pipeline-bulkrnaseq-sm_env` | `merged_counts.csv` (genes × samples) |
| 10 | `sample_check` | R | `dge_env` | Sex-marker heatmap, sample-correlation heatmap, mitochondrial/ribosomal barplot |
| 11 | `expression_check` | R | `dge_env` | MA plot, expression heatmap, barplot/violin plot of user-selected genes |
| 12 | `run_pca` | R | `dge_env` | `pca.csv` + metadata-correlation plots |
| 13 | `run_pca_category` | R | `dge_env` | One PCA plot per categorical metadata column |
| 14 | `run_pca_numeric` | R | `dge_env` | One PCA plot per numeric metadata column |
| 15 | `run_dge` | R (DESeq2) | `dge_env` | Per comparison: `dge_results.csv`, `vst_normalised_counts.csv`, volcano plot |
| 16 | `run_pathway` | R (clusterProfiler, ReactomePA) | `dge_env` | Per comparison: GO / KEGG / Reactome GSEA tables + dotplot |

The QC rules are wired into the dependency graph as gates, not just as reports:
`trim_reads` waits on the raw MultiQC report, `run_star` waits on the trimmed MultiQC
report, and `run_pca` waits on the sex and sample-correlation plots from `sample_check`.
Deleting a report therefore forces everything downstream of it to re-run.

Rules 9–16 are executed once **per group** (see [Grouping](#grouping)); rules 15–16 are
additionally executed once **per comparison** in `dge_comparisons`.

---

## Requirements

- Snakemake ≥ 7
- Conda / Mamba
- SLURM (for cluster execution)
- A pre-built STAR genome index
- A gene annotation TSV (see [Inputs](#inputs))
- Network access from the compute node for KEGG GSEA (`gseKEGG` downloads pathway data at run time)

Tools are installed via the two environment files in `envs/`:

| Environment | File | Contents |
|-------------|------|----------|
| `pipeline-bulkrnaseq-sm_env` | `envs/pipeline-bulkrnaseq-sm_env.yml` | FastQC 0.12.1, MultiQC 1.30, Cutadapt 5.1, STAR 2.7.11b, samtools 1.22, seqtk, subread, pandas, matplotlib, scikit-learn |
| `dge_env` | `envs/dge_env.yml` | R 4.3, DESeq2, tidyverse, optparse, ComplexHeatmap, circlize, pheatmap, ggrepel, patchwork, viridis, vcd, biomaRt, clusterProfiler, ReactomePA, AnnotationDbi, org.Mm.eg.db |

---

## Installation

Install conda:

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh
bash ~/miniconda.sh
source ~/.bashrc
```

Install Snakemake:

```bash
conda install -c conda-forge -c bioconda snakemake
```

Create **both** environments. The `Snakefile` references them **by name**, not by path,
so they must exist under those exact names:

```bash
conda env create -f envs/pipeline-bulkrnaseq-sm_env.yml
conda env create -f envs/dge_env.yml
```

If you created `dge_env` before the pathway packages were added, update it in place:

```bash
conda env update -n dge_env -f envs/dge_env.yml
```

Verify:

```bash
conda env list | grep -E "pipeline-bulkrnaseq-sm_env|dge_env"
conda run -n dge_env Rscript -e 'library(clusterProfiler); library(org.Mm.eg.db); library(ReactomePA)'
```

---

## Inputs

```
pipeline-bulkrnaseq-sm/
├── inputs/
│   ├── filepaths.csv   # sample ID → FASTQ paths
│   └── metadata.csv    # sample annotations
└── raw_data/           # FASTQ files (or point filepaths.csv anywhere)
```

### 1. `inputs/filepaths.csv`

Three columns: sample ID, R1 path, R2 path. Paths must be absolute.

```csv
Barcode,fastq_1,fastq_2
UDI001,/abs/path/SLX-11111.UDI001.AAA111.r_1.fq.gz,/abs/path/SLX-11111.UDI001.AAA111.r_2.fq.gz
UDI002,/abs/path/SLX-11111.UDI002.AAA111.r_1.fq.gz,/abs/path/SLX-11111.UDI002.AAA111.r_2.fq.gz
```

The `Snakefile` parses this file **positionally** and skips the first line, so the header
names are cosmetic — but column *order* matters. The first column becomes the sample
wildcard used in every output filename, and it **must match the `Barcode` column of
`metadata.csv`**.

Generate this file automatically from a FASTQ directory:

```bash
bash prerun_filepaths.sh
```

which calls `scripts/utils-generate_filepaths.sh ./raw_data ./inputs filepaths.csv`.
It recognises both `*.r_1.fq.gz` / `*.r_2.fq.gz` and `*_1.fq.gz` / `*_2.fq.gz` naming.

### 2. `inputs/metadata.csv`

```csv
Pool,Barcode,Sequence,Sample Name,Condition 1,Condition 2,Batch,Sample type,Sample age,Sample Sex
SLX-11111,UDI001,CCGCGGTT-AGCGCTAG,W1,treated,obese,1,tissue,adult,M
SLX-11111,UDI002,TTATAACC-GATATCGA,W2,untreated,obese,2,tissue,adult,F
```

Hard requirements:

- A **`Barcode`** column whose values match the first column of `filepaths.csv`. Every
  script keys on it (`merge_counts.py`, `run_pca.R`, `plot_sex.R`, `run_deseq_dge.R`, …).
- The header row must contain the strings **`Sample Name`** and **`Sample Sex`**
  (case-sensitive). `scripts/merge_counts.py` locates the header by scanning for a line
  containing both, which lets it tolerate preamble rows above the real header (as produced
  by some sequencing-facility spreadsheets).
- Any column named in `dge_cat`, `mapping_check_cat`, `sample_check_cat`, `split_by`,
  `pca_col_cat` or `pca_col_num` must exist.
- The levels of `dge_cat` must include both sides of every comparison in `dge_comparisons`.

Column names containing spaces are supported. Where a column name appears in an output
*filename*, spaces are replaced with underscores (`Condition 1` → `pca_cat_Condition_1.png`).

### 3. Gene annotation TSV (`annotation:`)

Tab-separated, with a header. Required columns:

| Column | Required | Purpose |
|--------|----------|---------|
| `ensembl_gene_id` | yes | Joined against the count-matrix row names |
| `gene_symbol` | yes | Gene labels on plots, in DGE results, and for GSEA ID mapping |
| `protein_coding` | optional | If present, only rows where this is `TRUE` are kept |

Counts are filtered to the gene IDs present in this file, so an annotation whose IDs do
not match your STAR index will silently produce an empty matrix.

---

## Configuration reference

Set in `config.yml` and passed with `--configfile`. Values can be overridden at runtime
with `--config key=value`.

### Paths and layout

| Key | Example | Description |
|-----|---------|-------------|
| `run_dir` | `/rds/project/.../pipeline-bulkrnaseq-sm` | Root for all outputs; `filepaths` and `metadata` are resolved relative to it |
| `pipeline_dir` | — | Repository root — **must be passed via `--config`**, not the YAML; used to locate `scripts/` |
| `dir_flag` | `results_raw` | Subdirectory of `run_dir` holding this run's results. Change it to keep several runs side by side |
| `filepaths` | `inputs/filepaths.csv` | Relative to `run_dir` |
| `metadata` | `inputs/metadata.csv` | Relative to `run_dir` |
| `star_index` | `/.../ensembl/homo_sapiens/release-113` | Pre-built STAR genome directory |
| `annotation` | `/.../gene_annotation.tsv` | Gene annotation TSV |

### QC and trimming

| Key | Values | Description |
|-----|--------|-------------|
| `multiqc_mode` | `interactive` \| `flat` | `flat` renders static plots; use it for large cohorts where the interactive report becomes unusable |
| `trimming_method` | `illumina` \| `nextera` \| `first` | `illumina` = Cutadapt, TruSeq adapter `AGATCGGAAGAGC`, min length 25. `nextera` = Cutadapt, Nextera/Tn5 adapter `CTGTCTCTTATACACATCT` on R1 and R2, min length 25. `first` = hard-trim a fixed number of bases (hard-coded in the `Snakefile` as 5 bp from R2) |
| `mapping_check_cat` | `"Condition 1"` | Metadata column used to group samples in the mapping barplot |
| `sample_check_cat` | `"Condition 1"` | Metadata column used to annotate the sample-correlation heatmap and mito/ribo barplot |
| `gene_list` | `"Fgf21,Pparg,Adipoq"` | Comma-separated gene symbols for the `barplot_genes.png` plot. One gene → per-sample barplot; several → per-gene violin plots. **Effectively required** — the plot script stops if it is empty |
| `fastqc_threads` | integer | Passed via `--config` at runtime |
| `star_threads` | integer | Passed via `--config` at runtime |

### Alignment (STAR)

| Key | Example | Description |
|-----|---------|-------------|
| `sam_format` | `BAM` | `--outSAMtype` format |
| `sam_sort` | `SortedByCoordinate` | `--outSAMtype` sort order |
| `sam_strand_field` | `intronMotif` | `--outSAMstrandField`; needed by downstream tools that infer strand |
| `filter_score` | `0.1` | `--outFilterScoreMinOverLread` (see [Known issues](#known-issues)) |
| `filter_match` | `0.1` | `--outFilterMatchNminOverLread` (see [Known issues](#known-issues)) |
| `quant_mode` | `GeneCounts` | `--quantMode`; `GeneCounts` is required — it produces the `ReadsPerGene.out.tab` files that `merge_counts` consumes |

### Counting

| Key | Values | Description |
|-----|--------|-------------|
| `strandness` | `1` / `unstranded`, `2` / `forward`, `3` / `reverse` | Selects which column of STAR's `ReadsPerGene.out.tab` to read. **These are 1-based STAR count columns, not the 0/1/2 featureCounts convention.** Use `3` for reverse-stranded libraries (e.g. SMARTer Stranded, TruSeq Stranded mRNA), `2` for forward-stranded, `1` for unstranded |

### Grouping

| Key | Example | Description |
|-----|---------|-------------|
| `group_mode` | `false` | When `false`, every rule from `merge_counts` onward runs once over all samples, in a directory named `all/`. When `true`, it runs once per entry in `groups` |
| `split_by` | `Sample Sex` | Metadata column used to split samples when `group_mode: true` |
| `groups` | `[F, M]` | Values of `split_by` to build groups for. Ignored when `group_mode: false`. Directory names are sanitised (non-alphanumeric characters → `_`) |

Use grouping when you want fully independent analyses per stratum — for example separate
DESeq2 models for male and female samples — rather than a single model with sex as a
covariate.

### PCA

| Key | Example | Description |
|-----|---------|-------------|
| `norm` | `cpm` \| `vst` | Normalisation applied before PCA |
| `filt` | `TRUE` | Low-count filtering before PCA **and** DESeq2 |
| `pca_col_cat` | `"Condition 1,Sample Sex"` | Comma-separated categorical metadata columns; one PCA plot per column |
| `pca_col_cat_reference` | `"untreated,F"` | Reference level for each column in `pca_col_cat`, **in the same order and of the same length** |
| `pca_col_num` | `"Batch"` | Comma-separated numeric metadata columns; one PCA plot per column |

### Differential expression

| Key | Example | Description |
|-----|---------|-------------|
| `dge_cat` | `"Condition 1"` | Metadata column used as the DESeq2 design variable (`~ dge_cat`). Also used to colour the MA plot and expression heatmap |
| `dge_comparisons` | see below | Mapping of **comparison name → reference level**. One DESeq2 run, volcano plot and GSEA run is produced per entry |
| `lfc_threshold` | `1` | Absolute log2 fold-change cutoff for calling significance |
| `padj_threshold` | `0.05` | Adjusted p-value cutoff |

Comparison names must follow `NUMERATOR_vs_REFERENCE`, and the value must equal the
`REFERENCE` part — `run_deseq_dge.R` checks this and stops on a mismatch. Both levels must
exist in the `dge_cat` column.

```yaml
dge_cat: "Condition 1"
dge_comparisons:
  treated_vs_untreated: untreated
  # KO_vs_Con: Con
```

`dge_cat_reference` from earlier versions is no longer used.

### Pathway analysis (GSEA)

All optional; defaults shown.

| Key | Default | Description |
|-----|---------|-------------|
| `pathway_organism` | `mouse` | Only `mouse` is currently supported (uses `org.Mm.eg.db`, KEGG `mmu`) |
| `pathway_databases` | `[GO, KEGG, REACTOME]` | Databases to test. All three output CSVs are still expected by the `Snakefile` |
| `pathway_min_genes` | `10` | Minimum gene-set size |
| `pathway_max_genes` | `500` | Maximum gene-set size |
| `pathway_gsea_fdr` | `0.05` | FDR cutoff for reporting / plotting enriched sets |
| `pathway_lfc_threshold`, `pathway_padj_threshold` | `lfc_threshold`, `padj_threshold` | Read by the `Snakefile` but not currently passed to the script |

Genes are ranked by the DESeq2 Wald statistic (`stat`); symbols are mapped to Entrez IDs,
keeping the gene with the largest |stat| where several map to the same ID.

### Runtime-only keys (`--config`)

| Key | Description |
|-----|-------------|
| `pipeline_dir` | Repository root; `scripts/` is resolved from it |
| `run_dir` | Overrides `run_dir` from the YAML |
| `fastqc_threads`, `star_threads` | Thread counts per job |
| `exclude` | Comma-separated sample IDs to drop, e.g. `--config exclude="UDI003,UDI007"`. Excluded samples are removed from the sample list before the DAG is built, so no outputs are produced for them |

---

## Running the pipeline

### On SLURM

Edit `run.sh` so `PIPELINE_DIR` points at your checkout. `run.sh` runs Snakemake and
prints the total runtime at the end:

```bash
PIPELINE_DIR="/rds/project/rds-O11U8YqSuCk/bioinformatics/pipelines/pipeline-bulkrnaseq-sm"
CONFIG_FILE="$PIPELINE_DIR/config.yml"

snakemake -s "$PIPELINE_DIR/Snakefile" \
          -j 2 \
          --use-conda \
          --rerun-incomplete \
          --configfile "$CONFIG_FILE" \
          --config star_threads=4 fastqc_threads=4 pipeline_dir="$PIPELINE_DIR"
```

`run.sh` no longer carries `#SBATCH` headers, so pass the resources on the command line:

```bash
sbatch -J pipeline-bulkrnaseq-sm -A <ACCOUNT> -p sapphire \
       --cpus-per-task=16 --mem=64G --time=4:00:00 --mail-type=ALL \
       run.sh
```

or run it inside an interactive allocation (`sintr` / `salloc`).

`-j` is the number of *concurrent rule jobs*; `star_threads` and `fastqc_threads` are the
threads given to each job. Their product should not exceed `--cpus-per-task`. The shipped
values (`-j 2`, 4 threads each) use 8 of 16 CPUs — raise them together to use the full
allocation. STAR needs roughly 30 GB of RAM for a human index; `--mem=64G` covers two
concurrent STAR jobs.

### Useful invocations

```bash
# Dry run — print the DAG without executing anything
snakemake -s Snakefile -n --configfile config.yml --config pipeline_dir=$PWD star_threads=4 fastqc_threads=4

# Stop after alignment
snakemake -s Snakefile --use-conda --configfile config.yml \
  --config pipeline_dir=$PWD star_threads=4 fastqc_threads=4 \
  --until run_star

# Stop after DGE (skip pathway analysis)
snakemake ... --until run_dge

# Re-run one rule for all samples
snakemake -s Snakefile --use-conda --configfile config.yml \
  --config pipeline_dir=$PWD star_threads=4 fastqc_threads=4 \
  -R merge_counts

# Drop bad samples without editing filepaths.csv
snakemake ... --config exclude="UDI003,UDI007"

# Visualise the DAG
snakemake -s Snakefile --configfile config.yml --config pipeline_dir=$PWD \
  star_threads=4 fastqc_threads=4 --dag | dot -Tpng > dag.png
```

`--until <rule>` is the supported way to run only part of the workflow.

---

## Outputs

Everything is written under `<run_dir>/<dir_flag>/`:

```
<run_dir>/<dir_flag>/
├── fastq-raw/           # FastQC html + zip, raw reads
├── multiqc-raw/         # multiqc_report.html
├── trimmed/             # {sample}.r_1.fq.gz, {sample}.r_2.fq.gz
├── fastq-trimmed/       # FastQC html + zip, trimmed reads
├── multiqc-trimmed/     # multiqc_report.html
├── mapping/             # STAR BAM, ReadsPerGene.out.tab, SJ.out.tab, Log.* per sample
├── multiqc-mapped/      # multiqc_report.html (STAR logs)
├── mapping_check/       # barplot_mapping.png
├── merged-counts/<group>/
│   └── merged_counts.csv
├── sample_check/<group>/
│   ├── heatmap_sex.png
│   ├── correlation_sample.png
│   └── barplot_mito_ribo.png
├── expression_check/<group>/
│   ├── MA_plot.png
│   ├── heatmap_plot.png
│   └── barplot_genes.png
├── pca/<group>/
│   ├── pca.csv
│   ├── correlation_metadata.png
│   ├── correlation_metadata_pca.png
│   ├── pca_cat_<column>.png
│   └── pca_num_<column>.png
├── dge/<group>/<comparison>/
│   ├── dge_results.csv
│   ├── vst_normalised_counts.csv
│   └── volcano_<comparison>.png
└── pathway/<group>/<comparison>/
    ├── GSEA_gene_ranking.csv
    ├── GO_GSEA.csv
    ├── KEGG_GSEA.csv
    ├── REACTOME_GSEA.csv
    └── pathway_dotplot.png
```

`<group>` is `all` when `group_mode: false`, otherwise one directory per entry in `groups`.
`<comparison>` is each key of `dge_comparisons`.

Key files:

| File | Contents |
|------|----------|
| `merged_counts.csv` | Raw gene counts, genes (Ensembl IDs) × samples (Barcodes) |
| `barplot_mito_ribo.png` | Fraction of counts from mitochondrial (`mt-`) and ribosomal (`Rps`/`Rpl`) genes per sample |
| `pca.csv` | Principal components per sample, with a `Barcode` column |
| `dge_results.csv` | DESeq2 results: log2 fold change, Wald statistic, p-value, adjusted p-value, gene symbol, numerator/denominator |
| `vst_normalised_counts.csv` | Variance-stabilised counts, suitable for clustering and plotting |
| `*_GSEA.csv` | Enriched gene sets per database (empty file if the database was skipped or nothing was found) |

---

## Helper scripts

### Wired into the workflow

| Script | Called by | Notes |
|--------|-----------|-------|
| `run_fastqc.sh` | `fastqc`, `fastqc_trimmed` | `run_fastqc.sh <r1> <r2> <outdir> <threads>` |
| `run_multiqc.sh` | all MultiQC rules | `run_multiqc.sh <input_dir> <out_dir> <interactive\|flat>` |
| `run_cutadapt_illumina.sh` | `trim_reads` | TruSeq adapter, min length 25 |
| `run_cutadapt_nextera.sh` | `trim_reads` | Nextera adapter, min length 25 |
| `run_cutadapt_first.sh` | `trim_reads` | `<r1> <r2> <trim_bp> <r1\|r2\|both> <out1> <out2>` |
| `run_star.sh` | `run_star` | Twelve positional arguments; sets `ulimit -n 20000` |
| `merge_counts.py` | `merge_counts` | Strips STAR's `N_*` summary rows, selects the strand column, subsets to the group |
| `plot_barplot_mapping.py` | `mapping_check` | Parses `Log.final.out` into unique / multi / too-many-loci / unmapped |
| `plot_sex.R`, `plot_correlation_sample.R`, `plot_barplot_sample.R` | `sample_check` | Sex markers cover mouse (`Xist`, `Ddx3y`, `Eif2s3y`, `Uty`, …) and human (`XIST`, `DDX3Y`, …) symbols |
| `plot_ma.R`, `plot_heatmap.R`, `plot_barplot_genes.R` | `expression_check` | |
| `run_pca.R`, `plot_correlation_metadata.R`, `plot_correlation_metadata_pca.R` | `run_pca` | |
| `plot_pca_category.R`, `plot_pca_numeric.R` | `run_pca_category`, `run_pca_numeric` | |
| `run_deseq_dge.R` | `run_dge` | DESeq2, design `~ dge_cat` (or `~ study_col + dge_cat` with `--study_col`, not exposed in the `Snakefile`) |
| `run_gsea_pathway.R` | `run_pathway` | clusterProfiler `gseGO` / `gseKEGG`, ReactomePA `gsePathway` |

### Standalone utilities (not part of the DAG)

Run these manually before or alongside the pipeline:

| Script | Usage |
|--------|-------|
| `utils-generate_filepaths.sh` | `bash utils-generate_filepaths.sh <fastq_dir> <out_dir> filepaths.csv` — build `filepaths.csv` |
| `utils-merge_lanes.sh` | `bash utils-merge_lanes.sh <fastq_dir>` — concatenate per-lane FASTQs into `<fastq_dir>/merged_lanes/` |
| `utils-merge_meta.R` | `Rscript utils-merge_meta.R <contents.csv> <sample_info.csv> <output.csv>` — join facility contents to sample metadata on `Sample name` |
| `run_seqtk.sh` | `bash run_seqtk.sh <fastq_dir> <target_reads> [out_dir]` — downsample FASTQs to a fixed depth |
| `run_samtools.sh` | `bash run_samtools.sh <in.bam> <out_prefix> <threads> <annotation.gtf> <strand>` — deduplicate, index, featureCounts, flagstat/stats |
| `utils-samtools_stats.sh` | `bash utils-samtools_stats.sh <dedup.bam> [threads]` — flagstat + stats on primary alignments |
| `run_cutadapt_polyg.sh` | `bash run_cutadapt_polyg.sh <r1> <r2> <out1> <out2>` — trim poly-G tails (NovaSeq two-colour chemistry) |
| `plot_correlation_metadata_mapping.R` | Association between metadata columns and STAR uniquely mapped reads |
| `plot_barplot_metadata.R` | `Rscript plot_barplot_metadata.R --metadata <csv> --num_cols "a,b" --cat_cols "c,d" --out <png>` — overview of metadata distributions |

---

## Known issues

- **Metadata header case.** `merge_counts.py` now looks for `Sample Name` and `Sample Sex`
  (capitalised). The example `inputs/metadata.csv` still uses `Sample name` / `Sample sex`
  and will not be recognised.
- **STAR filter thresholds are truncated to 0.** The `Snakefile` reads `filter_score` and
  `filter_match` with `int(...)`, so `0.1` becomes `0` and those STAR filters are effectively
  disabled. Change them to `float(...)` to apply the configured values.
- **Pathway analysis is mouse-only** (`org.Mm.eg.db`, KEGG `mmu`). The shipped `config.yml`
  points at a human STAR index, so `run_pathway` will stop on human data — use
  `--until run_dge` until human support is added. KEGG GSEA also needs internet access at
  run time.
- **Mito/ribo barplot uses mouse gene symbols** (`^mt-`, `^Rps|^Rpl`); on human data
  (`MT-`, `RPS`/`RPL`) it will find no genes.
- **`pathway_databases` does not change the expected outputs** — all three `*_GSEA.csv`
  files are always required by `rule all`.
- **`trimming_method: first`** always trims 5 bp from R2; the value is hard-coded in the
  `Snakefile`.

---

## Author

Maintained by mk2314
University of Cambridge
