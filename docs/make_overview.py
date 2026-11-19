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
    "check":    ("#f1ebfb", "#7b55c7"),
    "analysis": ("#e4f4ea", "#2f8f57"),
}
GATE = C["qc"][1]

fig = plt.figure(figsize=(11.69, 8.27))
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 100); ax.set_ylim(0, 70.7); ax.axis("off")


def box(x, y, w, h, kind, title, tool, lines, title_fs=10.5):
    f, e = C[kind]
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle="round,pad=0,rounding_size=1.0",
                                fc=f, ec=e, lw=1.4))
    top = y + h / 2
    ax.text(x, top - 1.9, title, ha="center", va="center", fontsize=title_fs, weight="bold", color=INK)
    ax.text(x, top - 4.0, tool, ha="center", va="center", fontsize=8.5, color=e, weight="bold")
    ax.text(x, top - 5.6, "\n".join(lines), ha="center", va="top", fontsize=7.6, color=MUTED, linespacing=1.45)


def arrow(p, q, dashed=False, color=INK, rad=0.0):
    ax.annotate("", xy=q, xytext=p, arrowprops=dict(
        arrowstyle="-|>", lw=1.3, color=color, linestyle=(0, (3, 2)) if dashed else "-",
        shrinkA=0, shrinkB=0, mutation_scale=11, connectionstyle=f"arc3,rad={rad}"))


def elbow(points, color=INK):
    """Right-angled connector through points; arrowhead on the last segment."""
    xs_, ys_ = zip(*points[:-1])
    ax.plot(xs_, ys_, color=color, lw=1.3, solid_joinstyle="miter", solid_capstyle="butt")
    arrow(points[-2], points[-1], color=color)


# --- title ---------------------------------------------------------------
ax.text(4, 67.4, "pipeline-bulkrnaseq-sm", fontsize=20, weight="bold", color=INK, va="center")
ax.text(4, 64.2, "Paired-end bulk RNA-seq on SLURM with Snakemake: FASTQ to differential expression and pathways",
        fontsize=10.5, color=MUTED, va="center")

# --- stage bands ---------------------------------------------------------
bands = [(1.75, 16.25, "Inputs"), (19.25, 33.75, "1  Read QC + trimming"), (36.25, 50.75, "2  Alignment"),
         (53.25, 67.75, "3  Counts + sample QC"), (70.25, 98.25, "4  Differential expression")]
for x0, x1, label in bands:
    ax.add_patch(FancyBboxPatch((x0, 57.4), x1 - x0, 3.4, boxstyle="round,pad=0,rounding_size=0.8",
                                fc="#f4f6f8", ec="none"))
    ax.text((x0 + x1) / 2, 59.1, label, ha="center", va="center", fontsize=9, weight="bold", color=INK)

# --- main path -----------------------------------------------------------
Y, H, W = 45, 13.5, 14.5
xs = [9, 26.5, 43.5, 60.5, 77.5, 92]
WL = 12.5  # last box
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
box(xs[5], Y, WL, H, "analysis", "run_pathway", "clusterProfiler",
    ["GSEA per comparison", "GO · KEGG · Reactome", "dotplot (mouse)"])
for a, b in zip(xs[:-1], xs[1:]):
    arrow((a + W / 2, Y), (b - (WL if b == xs[-1] else W) / 2, Y))
MB = Y - H / 2  # bottom of main row

# --- QC boxes: row A (reports / checks) and row B (per-sample QC) -----------
AY, BY, QH, QW = 28.5, 15, 9.5, 12
qx = [17.75, 35, 50.5]          # raw, trimmed, mapped columns
A_top, A_bot = AY + QH / 2, AY - QH / 2
B_top, B_left = BY + QH / 2, [x - QW / 2 for x in qx]

box(qx[0], BY, QW, QH, "qc", "fastqc", "FastQC", ["raw reads", "per-sample report"], 9.5)
box(qx[0], AY, QW, QH, "qc", "multiqc", "MultiQC", ["raw reads", "combined report"], 9.5)
box(qx[1], BY, QW, QH, "qc", "fastqc_trimmed", "FastQC", ["trimmed reads", "per-sample report"], 9.5)
box(qx[1], AY, QW, QH, "qc", "multiqc_trimmed", "MultiQC", ["trimmed reads", "combined report"], 9.5)
box(qx[2], AY, QW, QH, "qc", "multiqc_mapped", "MultiQC", ["STAR logs + counts", "combined report"], 9.5)
box(qx[2], BY, QW, QH, "check", "mapping_check", "Python", ["mapping-rate", "barplot"], 9.5)

cx, CW = [64.75, 78.25, 91.75], 12.5
box(cx[0], AY, CW, QH, "check", "sample_check", "R", ["sex (Xist / Y genes)", "correlation, mito/ribo"], 9.5)
box(cx[1], AY, CW, QH, "qc", "run_pca", "R", ["PCA by metadata", "PCs vs metadata"], 9.5)
box(cx[2], AY, CW, QH, "check", "expression_check", "R", ["MA plot, heatmap", "gene_list plots"], 9.5)

# raw reads: inputs -> fastqc -> multiqc -> (gate) trim_reads
elbow([(5.5, MB), (5.5, BY), (B_left[0], BY)])
arrow((qx[0], B_top), (qx[0], A_bot))
arrow((qx[0] + 2, A_top), (xs[1] - 4, MB), dashed=True, color=GATE)
# trimmed reads: trim_reads -> fastqc_trimmed -> multiqc_trimmed -> (gate) run_star
elbow([(26.4, MB), (26.4, BY), (B_left[1], BY)])
arrow((qx[1], B_top), (qx[1], A_bot))
arrow((qx[1] + 2, A_top), (xs[2] - 4, MB), dashed=True, color=GATE)
# alignment: run_star -> multiqc_mapped, run_star -> mapping_check
arrow((47.5, MB), (47.5, A_top))
elbow([(42.6, MB), (42.6, BY), (B_left[2], BY)])
# counts: merge_counts -> each check; run_pca waits for sample_check
bus = 35.6
ax.plot([xs[3], xs[3], cx[-1]], [MB, bus, bus], color=INK, lw=1.3, solid_joinstyle="miter")
for c in cx:
    arrow((c, bus), (c, A_top))
arrow((cx[0], A_bot), (cx[1], A_bot), dashed=True, color=GATE, rad=0.45)

# --- scope brackets ------------------------------------------------------
def bracket(x0, x1, y, label):
    ax.plot([x0, x0, x1, x1], [y - 0.8, y, y, y - 0.8], color=MUTED, lw=1)
    ax.text((x0 + x1) / 2, y + 0.9, label, ha="center", va="bottom", fontsize=7.8, color=MUTED, style="italic")

bracket(xs[1] - W / 2, xs[2] + W / 2, Y + H / 2 + 1.2, "once per sample")
bracket(xs[3] - W / 2, xs[5] + WL / 2, Y + H / 2 + 1.2, "once per group  (all samples, or per group with group_mode)")

# --- legend + footer -----------------------------------------------------
ax.add_patch(FancyBboxPatch((4, 7.3), 92, 0.08, boxstyle="square,pad=0", fc="#d5dbe2", ec="none"))
ly = 4.6
x = 4
for k, label in [("input", "input"), ("process", "processing"), ("qc", "QC report"), ("check", "check plots"), ("analysis", "analysis")]:
    f, e = C[k]
    ax.add_patch(FancyBboxPatch((x, ly - 1), 2.6, 2, boxstyle="round,pad=0,rounding_size=0.4", fc=f, ec=e, lw=1.2))
    ax.text(x + 3.4, ly, label, va="center", fontsize=8.5, color=INK)
    x += 11.8
arrow((x, ly), (x + 5, ly)); ax.text(x + 6, ly, "data flow", va="center", fontsize=8.5, color=INK)
x += 13.5
arrow((x, ly), (x + 5, ly), dashed=True, color=GATE)
ax.text(x + 6, ly, "gate: next step waits for this", va="center", fontsize=8.5, color=INK)

ax.text(4, 1.6, "Run:  sbatch submit_slurm.sh   (dry run: bash run.sh -n)    Results:  <run_dir>/<dir_flag>/    "
        "Test data: 8 mouse liver samples (PRJDB17560), ~12 min",
        fontsize=8, color=MUTED, va="center", family="DejaVu Sans Mono")

fig.savefig(out)
print("wrote", out)
