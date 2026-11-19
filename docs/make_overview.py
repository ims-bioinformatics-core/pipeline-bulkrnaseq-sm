"""Draw the pipeline overview. Usage: python docs/make_overview.py docs/pipeline_overview.pdf  (or .svg / .png)"""
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

out = sys.argv[1]
plt.rcParams.update({"font.family": "DejaVu Sans", "pdf.fonttype": 42})

INK = "#1f2933"; MUTED = "#5f6b7a"
C = {  # fill, edge
    "input":    ("#eef1f5", "#8a96a6"),
    "process":  ("#e3eefb", "#3b78c2"),
    "qc":       ("#fdf1df", "#d08a1e"),
    "analysis": ("#e4f4ea", "#2f8f57"),
}

fig = plt.figure(figsize=(11.69, 8.27))
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 100); ax.set_ylim(0, 70.7); ax.axis("off")


def box(x, y, w, h, kind, title, tool, lines, tool_fs=8.5):
    f, e = C[kind]
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle="round,pad=0,rounding_size=1.0",
                                fc=f, ec=e, lw=1.4))
    top = y + h / 2
    ax.text(x, top - 1.9, title, ha="center", va="center", fontsize=10.5, weight="bold", color=INK)
    ax.text(x, top - 4.0, tool, ha="center", va="center", fontsize=tool_fs, color=e, weight="bold")
    ax.text(x, top - 5.6, "\n".join(lines), ha="center", va="top", fontsize=7.6, color=MUTED, linespacing=1.45)


def arrow(p, q, dashed=False, color=INK, rad=0.0):
    ax.annotate("", xy=q, xytext=p, arrowprops=dict(
        arrowstyle="-|>", lw=1.3, color=color, linestyle=(0, (3, 2)) if dashed else "-",
        shrinkA=0, shrinkB=0, mutation_scale=11, connectionstyle=f"arc3,rad={rad}"))


# --- title ---------------------------------------------------------------
ax.text(4, 66.3, "pipeline-bulkrnaseq-sm", fontsize=20, weight="bold", color=INK, va="center")
ax.text(4, 62.9, "Paired-end bulk RNA-seq on SLURM with Snakemake: FASTQ to differential expression and pathways",
        fontsize=10.5, color=MUTED, va="center")

# --- stage bands ---------------------------------------------------------
bands = [(1.75, 16.25, "Inputs"), (19.25, 33.75, "1  Read QC + trimming"), (36.25, 50.75, "2  Alignment"),
         (53.25, 67.75, "3  Counts + sample QC"), (70.25, 98.25, "4  Differential expression")]
for x0, x1, label in bands:
    ax.add_patch(FancyBboxPatch((x0, 56.2), x1 - x0, 3.6, boxstyle="round,pad=0,rounding_size=0.8",
                                fc="#f4f6f8", ec="none"))
    ax.text((x0 + x1) / 2, 58.0, label, ha="center", va="center", fontsize=9, weight="bold", color=INK)

# --- main path (row 1) ---------------------------------------------------
Y, H, W = 43, 13.5, 14.5
xs = [9, 26.5, 43.5, 60.5, 77.5, 92]
box(xs[0], Y, W, H, "input", "Inputs", "you provide",
    ["raw_data/  FASTQ R1+R2", "inputs/filepaths.csv", "inputs/metadata.csv", "config.yml"])
box(xs[1], Y, W, H, "process", "trim_reads", "Cutadapt",
    ["remove adapters", "trimming_method:", "illumina | nextera | first"])
box(xs[2], Y, W, H, "process", "run_star", "STAR",
    ["align to genome", "gene counts per sample", "star_index"])
box(xs[3], Y, W, H, "process", "merge_counts", "Python",
    ["genes × samples matrix", "strandness: 1 | 2 | 3", "annotation"])
box(xs[4], Y, W, H, "analysis", "run_dge", "DESeq2",
    ["design ~ dge_cat", "one run per comparison", "results table, volcano"])
box(xs[5], Y, 12.5, H, "analysis", "run_pathway", "clusterProfiler",
    ["GSEA on DESeq2 stat", "GO · KEGG · Reactome", "dotplot (mouse)"])
for a, b in zip(xs[:-1], xs[1:]):
    wa = W if a != xs[-1] else 12.5
    wb = W if b != xs[-1] else 12.5
    arrow((a + wa / 2, Y), (b - wb / 2, Y))

# --- QC row (row 2) ------------------------------------------------------
QY, QH, QW = 21, 10.5, 14.5
qx = [17.75, 35, 52]
box(qx[0], QY, QW, QH, "qc", "fastqc + multiqc", "raw reads",
    ["per-sample FastQC", "combined MultiQC report"])
box(qx[1], QY, QW, QH, "qc", "fastqc + multiqc", "trimmed reads",
    ["checks trimming worked", "combined MultiQC report"])
box(qx[2], QY, QW, QH, "qc", "multiqc_mapped", "+ mapping_check",
    ["STAR log report", "mapping-rate barplot"])
box(74.5, QY, 25, QH, "qc", "Sample QC and PCA", "sample_check · expression_check · run_pca",
    ["sex (Xist / Y genes), sample correlation,", "mito/ribo fraction, MA plot, heatmap,",
     "selected genes, PCA by metadata"], tool_fs=7.6)

top = QY + QH / 2
arrow((xs[0] + 3, Y - H / 2), (qx[0] - 3, top))                       # inputs -> raw QC
arrow((qx[0] + 3, top), (xs[1] - 3, Y - H / 2), dashed=True, color=C["qc"][1])  # gate -> trim
arrow((xs[1] + 3, Y - H / 2), (qx[1] - 3, top))                       # trim -> trimmed QC
arrow((qx[1] + 3, top), (xs[2] - 3, Y - H / 2), dashed=True, color=C["qc"][1])  # gate -> STAR
arrow((xs[2] + 3, Y - H / 2), (qx[2] - 1, top))                       # STAR -> mapped QC
arrow((xs[3] + 3, Y - H / 2), (67, top))                              # counts -> sample QC

# --- scope brackets ------------------------------------------------------
def bracket(x0, x1, y, label):
    ax.plot([x0, x0, x1, x1], [y - 0.8, y, y, y - 0.8], color=MUTED, lw=1)
    ax.text((x0 + x1) / 2, y + 0.9, label, ha="center", va="bottom", fontsize=7.8, color=MUTED, style="italic")

bracket(xs[0] + 9.5, xs[2] + W / 2, Y + H / 2 + 1.2, "once per sample")
bracket(xs[3] - W / 2, xs[5] + 6.25, Y + H / 2 + 1.2, "once per group  (all samples, or per group with group_mode)")
bracket(xs[4] - W / 2, xs[5] + 6.25, Y - H / 2 - 1.6, "")
ax.text((xs[4] - W / 2 + xs[5] + 6.25) / 2, Y - H / 2 - 3.6, "once per comparison in dge_comparisons",
        ha="center", fontsize=7.8, color=MUTED, style="italic")

# --- legend + run box ----------------------------------------------------
ly = 5.8
items = [("input", "input"), ("process", "processing"), ("qc", "QC report"), ("analysis", "analysis")]
x = 4
for k, label in items:
    f, e = C[k]
    ax.add_patch(FancyBboxPatch((x, ly - 1), 2.6, 2, boxstyle="round,pad=0,rounding_size=0.4", fc=f, ec=e, lw=1.2))
    ax.text(x + 3.4, ly, label, va="center", fontsize=8.5, color=INK)
    x += 13
arrow((x, ly), (x + 5, ly)); ax.text(x + 6, ly, "data flow", va="center", fontsize=8.5, color=INK)
x += 15
arrow((x, ly), (x + 5, ly), dashed=True, color=C["qc"][1])
ax.text(x + 6, ly, "QC gate: next step waits for this report", va="center", fontsize=8.5, color=INK)

ax.add_patch(FancyBboxPatch((4, 9.6), 92, 0.08, boxstyle="square,pad=0", fc="#d5dbe2", ec="none"))
ax.text(4, 2.4, "Run:  sbatch submit_slurm.sh   (dry run: bash run.sh -n)    Results:  <run_dir>/<dir_flag>/    "
        "Test data: 8 mouse liver samples (PRJDB17560), ~12 min",
        fontsize=8, color=MUTED, va="center", family="DejaVu Sans Mono")

fig.savefig(out)
print("wrote", out)
