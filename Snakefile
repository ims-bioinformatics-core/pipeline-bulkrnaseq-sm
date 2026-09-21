import os
import re


# -----------------------------
# Run directory
# -----------------------------
RUN_DIR = config.get("run_dir", ".")  # default is current directory
SCRIPTS_DIR = os.path.join(config.get("pipeline_dir"), "scripts")
DIR_FLAG = config.get("dir_flag", "")


# -----------------------------
# Exclude samples (optional)
# -----------------------------
exclude_raw = config.get("exclude", "")
exclude_samples = set(
    x.strip() for x in exclude_raw.split(",") if x.strip()
)

print("[INFO] Excluding samples:", exclude_samples)


# -----------------------------
# Read filepaths and metadata
# -----------------------------
filepaths = config.get("filepaths")
if filepaths is None:
    raise ValueError(
        "Config key 'filepaths' not found. Make sure you are passing the correct config.yml with --configfile."
    )
filepaths = os.path.join(RUN_DIR, filepaths)

samples = {}
with open(filepaths) as f:
    header = f.readline()
    for line in f:
        parts = line.strip().split(",")
        if len(parts) < 3:
            continue

        sample = parts[0]

        # 🔥 APPLY EXCLUSION HERE
        if sample in exclude_samples:
            print(f"[INFO] Skipping excluded sample: {sample}")
            continue

        fq1 = parts[1]
        fq2 = parts[2]
        samples[sample] = {"r1": fq1, "r2": fq2}
        
metadata = os.path.join(RUN_DIR, config["metadata"])
# -----------------------------
# Directories for outputs under RUN_DIR
# -----------------------------

DIR_FLAG = config["dir_flag"]

fastqc_dir        = os.path.join(RUN_DIR, DIR_FLAG, "fastq-raw")
multiqc_dir       = os.path.join(RUN_DIR, DIR_FLAG, "multiqc-raw")
trim_dir          = os.path.join(RUN_DIR, DIR_FLAG, "trimmed")
fastqc_trim_dir   = os.path.join(RUN_DIR, DIR_FLAG, "fastq-trimmed")
multiqc_trim_dir  = os.path.join(RUN_DIR, DIR_FLAG, "multiqc-trimmed")
star_out_dir      = os.path.join(RUN_DIR, DIR_FLAG, "mapping")
multiqc_mapped_dir= os.path.join(RUN_DIR, DIR_FLAG, "multiqc-mapped")
merged_counts_dir = os.path.join(RUN_DIR, DIR_FLAG, "merged-counts")
plot_pca_dir      = os.path.join(RUN_DIR, DIR_FLAG, "pca")
plots_dir    = os.path.join(RUN_DIR, DIR_FLAG, "plots")
sample_check_dir  = os.path.join(RUN_DIR, DIR_FLAG, "sample_check")
expression_check_dir  = os.path.join(RUN_DIR, DIR_FLAG, "expression_check")
mapping_check_dir = os.path.join(RUN_DIR, DIR_FLAG, "mapping_check")
dge_dir = os.path.join(RUN_DIR, DIR_FLAG, "dge")


# -----------------------------
# Parameters
# -----------------------------
fastqc_threads    = int(config["fastqc_threads"])
multiqc_mode      = config["multiqc_mode"]
trimming_method   = config["trimming_method"]
star_threads      = int(config["star_threads"])
star_index_dir    = config["star_index"]
annotation        = config["annotation"]
strandness        = config["strandness"]
sam_format        = config["sam_format"]
sam_sort          = config["sam_sort"]
sam_strand_field  = config["sam_strand_field"]
filter_score      = int(config["filter_score"])
filter_match      = int(config["filter_match"])
quant_mode        = config["quant_mode"]
norm              = config["norm"]
lfc_threshold     = config["lfc_threshold"]
padj_threshold    = config["padj_threshold"]
filt              = config["filt"]
dge_cat           = config["dge_cat"]
dge_comparisons   = config["dge_comparisons"]
dge_comparison_names = list(dge_comparisons.keys())
mapping_check_cat = config["mapping_check_cat"]


# -----------------------------
# PCA settings
# -----------------------------
pca_col_cat = [x.strip() for x in config["pca_col_cat"].split(",")]
pca_col_num = [x.strip() for x in config["pca_col_num"].split(",")]
pca_col_cat_ref = [x.strip() for x in config["pca_col_cat_reference"].split(",")]

pairs = list(zip(pca_col_cat, pca_col_cat_ref))
cat_pairs = list(zip(
    [c.replace(" ", "_") for c, r in pairs],
    [r for c, r in pairs]
))

group_mode = config.get("group_mode", False)
split_by = config.get("split_by", "")
raw_groups = config.get("groups", [])

def safe_name(x):
    return re.sub(r"[^\w.-]+", "_", x.strip()).strip("_")

if group_mode:
    group_map = {
        safe_name(group): group
        for group in raw_groups
    }
    groups = list(group_map.keys())
else:
    group_map = {"all": "all"}
    groups = ["all"]

print("[INFO] group_mode:", group_mode)
print("[INFO] split_by:", split_by)
print("[INFO] groups:", groups)
print("[INFO] group_map:", group_map)

# -----------------------------
# Rule all
# -----------------------------
rule all:
    input:
        # FastQC raw
        expand(os.path.join(fastqc_dir, "{sample}.r_1_fastqc.html"), sample=samples.keys()),
        expand(os.path.join(fastqc_dir, "{sample}.r_1_fastqc.zip"), sample=samples.keys()),
        expand(os.path.join(fastqc_dir, "{sample}.r_2_fastqc.html"), sample=samples.keys()),
        expand(os.path.join(fastqc_dir, "{sample}.r_2_fastqc.zip"), sample=samples.keys()),

        # MultiQC raw
        os.path.join(multiqc_dir, "multiqc_report.html"),

        # Trimmed
        expand(os.path.join(trim_dir, "{sample}.r_1.fq.gz"), sample=samples.keys()),
        expand(os.path.join(trim_dir, "{sample}.r_2.fq.gz"), sample=samples.keys()),

        # FastQC trimmed
        expand(os.path.join(fastqc_trim_dir, "{sample}.r_1_fastqc.html"), sample=samples.keys()),
        expand(os.path.join(fastqc_trim_dir, "{sample}.r_1_fastqc.zip"), sample=samples.keys()),
        expand(os.path.join(fastqc_trim_dir, "{sample}.r_2_fastqc.html"), sample=samples.keys()),
        expand(os.path.join(fastqc_trim_dir, "{sample}.r_2_fastqc.zip"), sample=samples.keys()),

        # MultiQC trimmed
        os.path.join(multiqc_trim_dir, "multiqc_report.html"),

        # STAR
        expand(os.path.join(star_out_dir, "{sample}.Aligned.sortedByCoord.out.bam"), sample=samples.keys()),
        expand(os.path.join(star_out_dir, "{sample}.ReadsPerGene.out.tab"), sample=samples.keys()),

        # MultiQC mapped
        os.path.join(multiqc_mapped_dir, "multiqc_report.html"),
        
        # Mapping check
        os.path.join(mapping_check_dir,"barplot_mapping.png"),

        expand(os.path.join(merged_counts_dir, "{group}", "merged_counts.csv"),group=groups),
        expand(os.path.join(plot_pca_dir, "{group}", "pca.csv"),group=groups),
        expand(os.path.join(plot_pca_dir, "{group}", "correlation_metadata_pca.png"),group=groups),    
        expand(os.path.join(plot_pca_dir, "{group}", "correlation_metadata.png"),group=groups),
        
        # Sample check
        expand(os.path.join(sample_check_dir, "{group}", "heatmap_sex.png"),group=groups),
        expand(os.path.join(sample_check_dir, "{group}", "correlation_sample.png"),group=groups),
        expand(os.path.join(sample_check_dir, "{group}", "barplot_mito_ribo.png"),group=groups),
             
        # Expression check
        expand(os.path.join(expression_check_dir, "{group}", "MA_plot.png"),group=groups),
        expand(os.path.join(expression_check_dir, "{group}", "heatmap_plot.png"),group=groups),
        

        # PCA category plots 
        expand(os.path.join(plot_pca_dir, "{group}", "pca_cat_{col_cat}.png"),group=groups,col_cat=[c.replace(" ", "_") for c in pca_col_cat]),     
        
        # PCA numeric plots
        expand(os.path.join(plot_pca_dir, "{group}", "pca_num_{col_num}.png"),group=groups,col_num=[c.replace(" ", "_") for c in pca_col_num]),
        
        # DGE
        expand(os.path.join(dge_dir, "{group}", "{comparison}", "dge_results.csv"),group=groups,comparison=dge_comparison_names),
        expand(os.path.join(dge_dir, "{group}", "{comparison}", "vst_normalised_counts.csv"),group=groups,comparison=dge_comparison_names),
        expand(os.path.join(dge_dir, "{group}", "{comparison}", "volcano_{comparison}.png"),group=groups,comparison=dge_comparison_names)

        
# -----------------------------
# FastQC raw
# -----------------------------
rule fastqc:
    input:
        r1=lambda wc: samples[wc.sample]["r1"],
        r2=lambda wc: samples[wc.sample]["r2"]
    output:
        r1_html=os.path.join(fastqc_dir, "{sample}.r_1_fastqc.html"),
        r1_zip=os.path.join(fastqc_dir, "{sample}.r_1_fastqc.zip"),
        r2_html=os.path.join(fastqc_dir, "{sample}.r_2_fastqc.html"),
        r2_zip=os.path.join(fastqc_dir, "{sample}.r_2_fastqc.zip")
    threads: 
        fastqc_threads,
    conda:
        "pipeline-bulkrnaseq-sm_env"
    shell:
        """
        mkdir -p {fastqc_dir}
        bash {SCRIPTS_DIR}/run_fastqc.sh {input.r1} {input.r2} {fastqc_dir} {fastqc_threads}

        # For R1
        r1base=$(basename "{input.r1}" .fq.gz)
        mv "{fastqc_dir}/${{r1base}}_fastqc.html" "{output.r1_html}"
        mv "{fastqc_dir}/${{r1base}}_fastqc.zip"  "{output.r1_zip}"

        # For R2
        r2base=$(basename "{input.r2}" .fq.gz)
        mv "{fastqc_dir}/${{r2base}}_fastqc.html" "{output.r2_html}"
        mv "{fastqc_dir}/${{r2base}}_fastqc.zip"  "{output.r2_zip}"
        """

# -----------------------------
# MultiQC raw
# -----------------------------
rule multiqc:
    input:
        expand(os.path.join(fastqc_dir, "{sample}.r_1_fastqc.zip"), sample=samples.keys()) +
        expand(os.path.join(fastqc_dir, "{sample}.r_2_fastqc.zip"), sample=samples.keys())
    output:
        os.path.join(multiqc_dir, "multiqc_report.html")
    conda:
        "pipeline-bulkrnaseq-sm_env"
    shell:
        """
        mkdir -p {multiqc_dir}
        bash {SCRIPTS_DIR}/run_multiqc.sh {fastqc_dir} {multiqc_dir} {multiqc_mode}
        """

# -----------------------------
# Cutadapt trimming
# -----------------------------    

rule trim_reads:
    input:
        r1=lambda wc: samples[wc.sample]["r1"],
        r2=lambda wc: samples[wc.sample]["r2"],
        multiqc_report=os.path.join(multiqc_dir, "multiqc_report.html")
    output:
        r1=os.path.join(trim_dir, "{sample}.r_1.fq.gz"),
        r2=os.path.join(trim_dir, "{sample}.r_2.fq.gz")
    conda:
        "pipeline-bulkrnaseq-sm_env"
    params:
        method=config["trimming_method"]
    shell:
        """
        mkdir -p {trim_dir}

        if [ "{params.method}" = "first" ]; then
            echo "Running fixed-base trimming"
            bash {SCRIPTS_DIR}/run_cutadapt_first.sh \
                {input.r1} {input.r2} 5 r2 \
                {output.r1} {output.r2}

        elif [ "{params.method}" = "illumina" ]; then
            echo "Running Illumina adapter trimming"
            bash {SCRIPTS_DIR}/run_cutadapt_illumina.sh \
                {input.r1} {input.r2} \
                {output.r1} {output.r2}

        elif [ "{params.method}" = "nextera" ]; then
            echo "Running Nextera adapter trimming"
            bash {SCRIPTS_DIR}/run_cutadapt_nextera.sh \
                {input.r1} {input.r2} \
                {output.r1} {output.r2}

        else
            echo "❌ Unknown trimming method: {params.method}"
            exit 1
        fi
        """

# -----------------------------
# FastQC trimmed
# -----------------------------
rule fastqc_trimmed:
    input:
        r1=lambda wc: os.path.join(trim_dir, f"{wc.sample}.r_1.fq.gz"),
        r2=lambda wc: os.path.join(trim_dir, f"{wc.sample}.r_2.fq.gz")
    output:
        r1_html=os.path.join(fastqc_trim_dir, "{sample}.r_1_fastqc.html"),
        r1_zip=os.path.join(fastqc_trim_dir, "{sample}.r_1_fastqc.zip"),
        r2_html=os.path.join(fastqc_trim_dir, "{sample}.r_2_fastqc.html"),
        r2_zip=os.path.join(fastqc_trim_dir, "{sample}.r_2_fastqc.zip")
    threads: 
        fastqc_threads,
    conda:
        "pipeline-bulkrnaseq-sm_env"
    shell:
        """
        mkdir -p {fastqc_trim_dir}
        bash {SCRIPTS_DIR}/run_fastqc.sh {input.r1} {input.r2} {fastqc_trim_dir} {fastqc_threads}
        """

# -----------------------------
# MultiQC trimmed
# -----------------------------
rule multiqc_trimmed:
    input:
        expand(os.path.join(fastqc_trim_dir, "{sample}.r_1_fastqc.zip"), sample=samples.keys()) +
        expand(os.path.join(fastqc_trim_dir, "{sample}.r_2_fastqc.zip"), sample=samples.keys())
    output:
        os.path.join(multiqc_trim_dir, "multiqc_report.html")
    conda:
        "pipeline-bulkrnaseq-sm_env"
    shell:
        """
        mkdir -p {multiqc_trim_dir}
        bash {SCRIPTS_DIR}/run_multiqc.sh {fastqc_trim_dir} {multiqc_trim_dir} {multiqc_mode} 
        """
# -----------------------------
# STAR alignment
# -----------------------------

rule run_star:
    input:
        r1=lambda wc: os.path.join(trim_dir, f"{wc.sample}.r_1.fq.gz"),
        r2=lambda wc: os.path.join(trim_dir, f"{wc.sample}.r_2.fq.gz"),
        qc=os.path.join(multiqc_trim_dir, "multiqc_report.html")
    output:
        bam=os.path.join(star_out_dir, "{sample}.Aligned.sortedByCoord.out.bam"),
        gene_counts=os.path.join(star_out_dir, "{sample}.ReadsPerGene.out.tab"),
        sj_out=os.path.join(star_out_dir, "{sample}.SJ.out.tab"),
        log_final=os.path.join(star_out_dir, "{sample}.Log.final.out"),
        log_progress=os.path.join(star_out_dir, "{sample}.Log.progress.out"),
        log_gene=os.path.join(star_out_dir, "{sample}.Log.out"),
    params:
        out_prefix=lambda wc: wc.sample,
    threads: 
        star_threads,
    conda:
        "pipeline-bulkrnaseq-sm_env"
    shell:
        """
        mkdir -p {star_out_dir}
        bash {SCRIPTS_DIR}/run_star.sh {star_threads} {star_index_dir} {input.r1} {input.r2} {star_out_dir} {params.out_prefix} {sam_format} {sam_sort} {sam_strand_field} {filter_score} {filter_match} {quant_mode}
        """
        
# -----------------------------
# MultiQC mapped
# -----------------------------

rule multiqc_mapped:
    input:
        # collect STAR output files to summarize
        expand(os.path.join(star_out_dir, "{sample}.Log.final.out"), sample=samples.keys()),
        expand(os.path.join(star_out_dir, "{sample}.ReadsPerGene.out.tab"), sample=samples.keys())
    output:
        os.path.join(multiqc_mapped_dir, "multiqc_report.html")
    conda:
        "pipeline-bulkrnaseq-sm_env"
    shell:
        """
        mkdir -p {multiqc_mapped_dir}
        bash {SCRIPTS_DIR}/run_multiqc.sh {star_out_dir} {multiqc_mapped_dir} {multiqc_mode} 
        """

# -----------------------------
# Mapping check
# -----------------------------
rule mapping_check:
    input:
        logs=expand(os.path.join(star_out_dir, "{sample}.Log.final.out"), sample=samples.keys()),
        metadata=metadata

    output:
        plot=os.path.join(mapping_check_dir, "barplot_mapping.png")

    params:
        mapping_check_cat=mapping_check_cat

    conda:
        "pipeline-bulkrnaseq-sm_env"

    shell:
        """
        mkdir -p {mapping_check_dir}

        python3 {SCRIPTS_DIR}/plot_barplot_mapping.py \
            --logs {input.logs} \
            --metadata "{input.metadata}" \
            --mapping_check_cat "{params.mapping_check_cat}" \
            --out "{output.plot}"
        """

# -----------------------------
# Merge counts (by group)
# -----------------------------
rule merge_counts:
    input:
        counts=expand(
            os.path.join(star_out_dir, "{sample}.ReadsPerGene.out.tab"),
            sample=samples.keys()
        ),
        metadata=metadata

    output:
        merged=os.path.join(
            merged_counts_dir,
            "{group}",
            "merged_counts.csv"
        )

    params:
        strandness=strandness,
        split_by=split_by,
        group_name=lambda wc: group_map[wc.group]

    conda:
        "pipeline-bulkrnaseq-sm_env"

    shell:
        """
        mkdir -p "{merged_counts_dir}/{wildcards.group}"

        echo "[DEBUG] directory group: {wildcards.group}"
        echo "[DEBUG] metadata group: {params.group_name}"

        python3 {SCRIPTS_DIR}/merge_counts.py \
            --strandness {params.strandness} \
            --split_by "{params.split_by}" \
            --group "{params.group_name}" \
            --metadata "{input.metadata}" \
            --counts {input.counts} \
            --out "{output.merged}"
        """



# -----------------------------
# Sample QC Check
# -----------------------------
rule sample_check:
    input:
        counts=os.path.join(merged_counts_dir, "{group}", "merged_counts.csv"),
        metadata=metadata,
        annotation=annotation

    output:
        sex_plot=os.path.join(sample_check_dir, "{group}", "heatmap_sex.png"),
        corr_plot=os.path.join(sample_check_dir, "{group}", "correlation_sample.png"),
        mito_ribo_plot=os.path.join(sample_check_dir, "{group}", "barplot_mito_ribo.png"),

    params:
        dge_cat = config["dge_cat"],
        dge_comparisons = config["dge_comparisons"],
        dge_comparison_names = list(dge_comparisons.keys())

    conda:
        "dge_env"

    shell:
        """
        mkdir -p {sample_check_dir}/{wildcards.group}

        # -------------------------
        # SEX PLOT
        # -------------------------
        Rscript {SCRIPTS_DIR}/plot_sex.R \
            --counts "{input.counts}" \
            --metadata "{input.metadata}" \
            --annotation "{input.annotation}" \
            --out "{output.sex_plot}"

        # -------------------------
        # SAMPLE CORRELATION PLOT
        # -------------------------
        Rscript {SCRIPTS_DIR}/plot_correlation_sample.R \
            --counts "{input.counts}" \
            --metadata "{input.metadata}" \
            --annotation "{input.annotation}" \
            --group "{params.dge_cat}" \
            --out "{output.corr_plot}"
            
        # -------------------------
        # MITOCHONDRIAL / RIBOSOMAL QC
        # -------------------------
        Rscript {SCRIPTS_DIR}/plot_barplot_sample.R \
            --counts "{input.counts}" \
            --metadata "{input.metadata}" \
            --annotation "{input.annotation}" \
            --group "{params.dge_cat}" \
            --out "{output.mito_ribo_plot}"
        """

 
# -----------------------------
# Expression QC
# -----------------------------
rule expression_check:
    input:
        counts=os.path.join(merged_counts_dir, "{group}", "merged_counts.csv"),
        metadata=metadata,
        annotation=annotation

    output:
        expr_plot=os.path.join(expression_check_dir, "{group}", "MA_plot.png"),
        heat_plot=os.path.join(expression_check_dir, "{group}", "heatmap_plot.png"),

    conda:
        "dge_env"
        
    params:
        dge_cat=dge_cat

    shell:
        """
        mkdir -p {expression_check_dir}/{wildcards.group}

        # -------------------------
        # EXPRESSION PLOT
        # -------------------------
        Rscript {SCRIPTS_DIR}/plot_ma.R \
            --counts "{input.counts}" \
            --metadata "{input.metadata}" \
            --annotation "{input.annotation}" \
            --dge_cat "{params.dge_cat}" \
            --out "{output.expr_plot}"
            
        Rscript {SCRIPTS_DIR}/plot_heatmap.R \
            --counts "{input.counts}" \
            --metadata "{input.metadata}" \
            --annotation "{input.annotation}" \
            --dge_cat "{params.dge_cat}" \
            --out "{output.heat_plot}"
        """

# -----------------------------
# PCA (by group)
# -----------------------------
rule run_pca:
    input:
        sex_plot=os.path.join(sample_check_dir, "{group}", "heatmap_sex.png"),
        corr_plot=os.path.join(sample_check_dir, "{group}", "correlation_sample.png"),
        counts=os.path.join(merged_counts_dir, "{group}", "merged_counts.csv"),
        metadata=metadata,
        annotation=annotation
    output:
        pca=os.path.join(plot_pca_dir, "{group}", "pca.csv"),
        plot_cor_pca=os.path.join(plot_pca_dir, "{group}", "correlation_metadata_pca.png"),
        plot_cor_meta=os.path.join(plot_pca_dir, "{group}", "correlation_metadata.png")
    conda:
        "dge_env"
    shell:
        """
        mkdir -p {plot_pca_dir}/{wildcards.group}

        Rscript {SCRIPTS_DIR}/run_pca.R \
            --counts "{input.counts}" \
            --metadata "{input.metadata}" \
            --annotation "{input.annotation}" \
            --out "{output.pca}" \
            --norm "{norm}" \
            --filt "{filt}"

        Rscript {SCRIPTS_DIR}/plot_correlation_metadata_pca.R \
            --pca "{output.pca}" \
            --metadata "{input.metadata}" \
            --pca_col_num "{config[pca_col_num]}" \
            --pca_col_cat "{config[pca_col_cat]}" \
            --pca_col_ref "{config[pca_col_cat_reference]}" \
            --out "{output.plot_cor_pca}"
            
       Rscript {SCRIPTS_DIR}/plot_correlation_metadata.R \
            --pca "{output.pca}" \
            --metadata "{input.metadata}" \
            --pca_col_num "{config[pca_col_num]}" \
            --pca_col_cat "{config[pca_col_cat]}" \
            --pca_col_ref "{config[pca_col_cat_reference]}" \
            --out "{output.plot_cor_meta}"    
        """
        
# -----------------------------
# PCA plot (CATEGORY)
# -----------------------------
rule run_pca_category:
    input:
        pca=os.path.join(plot_pca_dir, "{group}", "pca.csv"),
        metadata=metadata

    output:
        plot_pca_cat=os.path.join(plot_pca_dir,"{group}","pca_cat_{col_cat}.png")

    conda:
        "dge_env"

    shell:
        """
        mkdir -p {plot_pca_dir}/{wildcards.group}

        Rscript {SCRIPTS_DIR}/plot_pca_category.R \
            --pca "{input.pca}" \
            --metadata "{input.metadata}" \
            --pca_col "{wildcards.col_cat}" \
            --out "{output.plot_pca_cat}"
        """

# -----------------------------
# PCA plot (NUMERIC)
# -----------------------------
rule run_pca_numeric:
    input:
        pca=os.path.join(plot_pca_dir, "{group}", "pca.csv"),
        metadata=metadata
    output:
        plot_pca_num=os.path.join(plot_pca_dir,"{group}", "pca_num_{col_num}.png")
    conda:
        "dge_env"
    shell:
        """
        mkdir -p {plot_pca_dir}/{wildcards.group}

        Rscript {SCRIPTS_DIR}/plot_pca_numeric.R \
            --pca "{input.pca}" \
            --metadata "{input.metadata}" \
            --pca_col "{wildcards.col_num}" \
            --out "{output.plot_pca_num}"
        """


def pca_plot_inputs(wildcards):
    return [
        *[
            os.path.join(
                plot_pca_dir, wildcards.group,
                f"pca_cat_{c.replace(' ', '_')}.png"
            )
            for c in pca_col_cat
        ],
        *[
            os.path.join(
                plot_pca_dir, wildcards.group,
                f"pca_num_{c.replace(' ', '_')}.png"
            )
            for c in pca_col_num
        ]
    ]


# -----------------------------
# DGE (DESeq2, by group)
# -----------------------------
rule run_dge:
    input:
        counts=os.path.join(
            merged_counts_dir,
            "{group}",
            "merged_counts.csv"
        ),
        metadata=metadata,
        annotation=annotation

    output:
        results=os.path.join(
            dge_dir,
            "{group}",
            "{comparison}",
            "dge_results.csv"
        ),
        vst=os.path.join(
            dge_dir,
            "{group}",
            "{comparison}",
            "vst_normalised_counts.csv"
        ),
        volcano=os.path.join(
            dge_dir,
            "{group}",
            "{comparison}",
            "volcano_{comparison}.png"
        )

    params:
        dge_cat=dge_cat,
        reference=lambda wc: dge_comparisons[wc.comparison],
        filt=filt,
        lfc_threshold=lfc_threshold,
        padj_threshold=padj_threshold

    conda:
        "dge_env"

    shell:
        """
        mkdir -p "{dge_dir}/{wildcards.group}/{wildcards.comparison}"

        Rscript {SCRIPTS_DIR}/run_deseq_dge.R \
            --counts "{input.counts}" \
            --metadata "{input.metadata}" \
            --annotation "{input.annotation}" \
            --out "{output.results}" \
            --comparison "{wildcards.comparison}" \
            --reference "{params.reference}" \
            --dge_cat "{params.dge_cat}" \
            --filt "{params.filt}" \
            --lfc_threshold "{params.lfc_threshold}" \
            --padj_threshold "{params.padj_threshold}"
        """

