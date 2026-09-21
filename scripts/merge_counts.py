#!/usr/bin/env python3

import sys
import os
import pandas as pd
import argparse
from io import StringIO

# -----------------------------
# Arguments
# -----------------------------
parser = argparse.ArgumentParser()

parser.add_argument("--strandness", required=True)
parser.add_argument("--metadata", required=True)
parser.add_argument("--group", required=True)
parser.add_argument("--split_by", required=True)
parser.add_argument("--out", required=True)
parser.add_argument("--counts", nargs="+", required=True)

args = parser.parse_args()

strandness = args.strandness
metadata_file = args.metadata
group = args.group
split_by = args.split_by
count_files = args.counts
output_file = args.out

# -----------------------------
# STAR strand mapping
# -----------------------------
strand_map = {
    "unstranded": 1,
    "forward": 2,
    "reverse": 3,
    "1": 1,
    "2": 2,
    "3": 3
}

if strandness not in strand_map:
    raise ValueError(f"Invalid strandness option: {strandness}")

col = strand_map[strandness]
print(f"[INFO] Using STAR column {col} for strandness={strandness}")

# -----------------------------
# Load metadata
# -----------------------------
def load_metadata(path):
    with open(path) as f:
        lines = f.readlines()

    start_idx = None
    for i, line in enumerate(lines):
        if "Sample Name" in line and "Sample Sex" in line:
            start_idx = i
            break

    if start_idx is None:
        raise ValueError("Could not find metadata header row")

    metadata_str = "".join(lines[start_idx:])
    df = pd.read_csv(StringIO(metadata_str), sep=",", engine="python")

    df.columns = (
        df.columns.astype(str)
        .str.strip()
        .str.replace(r"\s+", " ", regex=True)
    )

    df = df.loc[:, ~df.columns.str.contains("^Unnamed")]

    print("[DEBUG] metadata columns:", df.columns.tolist())

    return df


meta = load_metadata(metadata_file)

# -----------------------------
# REQUIRED columns
# -----------------------------
if split_by not in meta.columns:
    raise ValueError(f"Column '{split_by}' not found")

sample_col = "Barcode"

if sample_col not in meta.columns:
    raise ValueError("Barcode column not found")

meta[split_by] = meta[split_by].astype(str).str.strip()
meta[sample_col] = meta[sample_col].astype(str).str.strip()

# -----------------------------
# IMPORTANT FIX: group logic
# -----------------------------
if group.lower() == "all" or split_by.strip() == "":
    group_samples = set(meta[sample_col].tolist())
    print(f"[INFO] No grouping applied → using ALL samples ({len(group_samples)})")
else:
    group_samples = set(
        meta.loc[meta[split_by] == str(group), sample_col].tolist()
    )
    print(f"[INFO] Group '{group}' has {len(group_samples)} samples")

# -----------------------------
# merge counts
# -----------------------------
dfs = []

def norm(x):
    return str(x).strip()

for f in count_files:

    sample = norm(
        os.path.basename(f)
        .replace(".ReadsPerGene.out.tab", "")
        .split(".")[0]
    )

    if sample not in group_samples:
        continue

    df = pd.read_csv(f, sep="\t", header=None, comment="#")

    # remove STAR noise
    df = df[~df[0].astype(str).str.startswith("N_")]

    df = df.iloc[:, [0, col]]
    df.columns = ["gene_id", sample]
    df = df.set_index("gene_id")

    dfs.append(df)

if len(dfs) == 0:
    raise ValueError(f"No count files found for group '{group}'")

merged = pd.concat(dfs, axis=1)

os.makedirs(os.path.dirname(output_file), exist_ok=True)
merged.to_csv(output_file, sep=",")

print(f"[DONE] Wrote: {output_file}")