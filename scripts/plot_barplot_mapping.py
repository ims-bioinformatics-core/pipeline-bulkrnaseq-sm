#!/usr/bin/env python3

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

# -----------------------------
# Args
# -----------------------------
parser = argparse.ArgumentParser()
parser.add_argument("--logs", nargs="+", required=True)
parser.add_argument("--metadata", required=True)
parser.add_argument("--mapping_check_cat", required=True)
parser.add_argument("--out", required=True)

args = parser.parse_args()

# -----------------------------
# STAR log parser
# -----------------------------
def get_value(lines, key):
    for l in lines:
        if key in l:
            try:
                return float(l.split("|")[-1].strip())
            except:
                return 0.0
    return 0.0

rows = []

for f in args.logs:
    with open(f) as handle:
        lines = handle.readlines()

    sample = Path(f).name.replace(".Log.final.out", "")

    total = get_value(lines, "Number of input reads")
    unique = get_value(lines, "Uniquely mapped reads number")
    multi = get_value(lines, "Number of reads mapped to multiple loci")
    multi_toomany = get_value(lines, "Number of reads mapped to too many loci")
    unmapped = get_value(lines, "Number of reads unmapped")

    rows.append({
        "sample": sample,
        "unique": unique,
        "multi": multi,
        "multi_toomany": multi_toomany,
        "unmapped": unmapped,
        "total": total
    })

df = pd.DataFrame(rows)

# -----------------------------
# Metadata
# -----------------------------
meta = pd.read_csv(args.metadata)
meta.columns = meta.columns.str.strip()

if "Barcode" not in meta.columns:
    raise ValueError(f"'Barcode' column not found. Found: {list(meta.columns)}")

meta = meta.rename(columns={"Barcode": "sample"})
meta["sample"] = meta["sample"].astype(str)

df["sample"] = df["sample"].astype(str)
df = df.merge(meta, on="sample", how="left")

# -----------------------------
# Use CLI argument as grouping column
# -----------------------------
group_col = args.mapping_check_cat

if group_col not in df.columns:
    raise ValueError(
        f"'{group_col}' not found in metadata.\nAvailable columns:\n{list(df.columns)}"
    )

df["mapping_check_cat"] = df[group_col].astype(str).fillna("Unknown")

# -----------------------------
# Group ordering
# -----------------------------
group_order = sorted(df["mapping_check_cat"].unique())

# -----------------------------
# Build plotting order with gaps
# -----------------------------
plot_order = []
group_centers = {}
y_positions = []

gap = 1.2
y = 0

for g in group_order:
    sub = df[df["mapping_check_cat"] == g].sort_values("unique", ascending=False)

    start_y = y

    for _, row in sub.iterrows():
        plot_order.append(row)
        y_positions.append(y)
        y += 1

    end_y = y
    group_centers[g] = (start_y + end_y - 1) / 2 if len(sub) > 0 else start_y

    y += gap

df_plot = pd.DataFrame(plot_order)
df_plot["y"] = y_positions

# -----------------------------
# Plot
# -----------------------------
fig_height = max(3, len(df_plot) * 0.25)
plt.figure(figsize=(10, fig_height))

bottom = np.zeros(len(df_plot))

plt.barh(df_plot["y"], df_plot["unique"], label="Uniquely mapped")
bottom += df_plot["unique"]

plt.barh(df_plot["y"], df_plot["multi"], left=bottom,
         label="Multi-mapped (multiple loci)")
bottom += df_plot["multi"]

plt.barh(df_plot["y"], df_plot["multi_toomany"], left=bottom,
         label="Multi-mapped (too many loci)")
bottom += df_plot["multi_toomany"]

plt.barh(df_plot["y"], df_plot["unmapped"], left=bottom,
         label="Unmapped")

plt.yticks(df_plot["y"], df_plot["sample"])

# -----------------------------
# group labels (RIGHT side, no overlap)
# -----------------------------
x_max = df_plot[["unique", "multi", "multi_toomany", "unmapped"]].sum(axis=1).max()

for g, yc in group_centers.items():
    plt.text(
        x_max * 1.02,
        yc,
        g,
        ha="left",
        va="center",
        fontsize=10,
        fontweight="bold"
    )

plt.xlim(0, x_max * 1.15)

# -----------------------------
# Labels
# -----------------------------
plt.xlabel("Number of reads")
plt.title("STAR alignment summary (grouped by mapping_check_cat)")

# -----------------------------
# Legend (BOTTOM)
# -----------------------------
plt.legend(
    loc="upper center",
    bbox_to_anchor=(0.5, -0.08),
    ncol=4,
    frameon=False
)

plt.tight_layout()
plt.savefig(args.out, dpi=300, bbox_inches="tight")