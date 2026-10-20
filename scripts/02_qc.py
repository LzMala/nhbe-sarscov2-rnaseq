"""
Performs sample-level quality control on NHBE Series 1 RNA-seq samples:
    1. Library size (sequencing depth) distribution.
    2. Low-count gene filtering (>= 10 counts in >= 3 samples).
    3. Sample-to-sample Pearson correlation on log2-CPM.
    4. Principal Component Analysis (PCA) to verify mock vs. infected separation.
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

# Keep a gene with at least 10 counts in at least 3 samples, one full group.
MIN_COUNT = 10
MIN_SAMPLES = 3
MOCK_COLOR = "#4A90E2"
INFECTED_COLOR = "#D0021B"

# File paths
COUNTS_PATH = "data/processed/nhbe_series1_counts.tsv.gz"
META_PATH = "data/processed/metadata.csv"
FIG_DIR = "figures"
RES_DIR = "results"

os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(RES_DIR, exist_ok=True)

# 1. Load data and verify alignment
print(f"\nLoading counts from: {COUNTS_PATH}")
print(f"Loading metadata from: {META_PATH}")
counts = pd.read_csv(COUNTS_PATH, sep="\t", index_col=0)
meta = pd.read_csv(META_PATH)

# The count columns and the sample sheet must be in the same order.
if list(counts.columns) != list(meta["sample_id"]):
    print("Error: Count table columns do not match metadata sample_id order.")
    raise SystemExit(1)

# Clean labels for plots
labels = [
    c.replace("Series1_NHBE_", "").replace("_", " ")
    for c in counts.columns
]
# Blue = mock. Red = infected.
colors = [MOCK_COLOR if cond == "mock" else INFECTED_COLOR for cond in meta["condition"]]


### --------- Step 1: Library Size (Sequencing Depth) ---------

# Library size = add up every gene count in that sample.
library_sizes = counts.sum(axis=0)
print("\n--- Library Sizes (Total Read Counts) ---")
for label, size in zip(labels, library_sizes):
    print(f"{label}: {size:,} reads")

plt.figure(figsize=(7, 4.5))
plt.bar(labels, library_sizes / 1e6, color=colors, edgecolor="black", alpha=0.85)
plt.ylabel("Total Reads (Millions)", fontsize=11)
plt.title("Sequencing Depth per Sample", fontsize=12, fontweight="bold")
plt.xticks(rotation=25, ha="right", fontsize=10)
plt.grid(axis="y", linestyle="--", alpha=0.5)
plt.tight_layout()
plt.savefig(f"{FIG_DIR}/01_library_sizes.png", dpi=300)
plt.close()
print(f"Saved: {FIG_DIR}/01_library_sizes.png")

### --------- Step 2: Filter Lowly Expressed Genes ---------
# Keep genes with at least 10 counts in at least 3 samples (one full condition group)
keep_mask = (counts >= MIN_COUNT).sum(axis=1) >= MIN_SAMPLES
counts_filtered = counts.loc[keep_mask]
zero_genes = (counts.sum(axis=1) == 0).sum()

print("\n--- Gene Filtering Summary ---")
print(f"  Total genes annotated: {counts.shape[0]:,}")
print(f"  Genes with zero counts across all samples: {zero_genes:,}")
print(f"  Genes passing filter (>= {MIN_COUNT} counts in >= {MIN_SAMPLES} samples): {counts_filtered.shape[0]:,}\n")

filter_summary = pd.DataFrame({
    "metric": ["total_genes", "zero_genes", "min_count", "min_samples", "genes_retained"],
    "value": [counts.shape[0], zero_genes, MIN_COUNT, MIN_SAMPLES, counts_filtered.shape[0]]
})
print(filter_summary)
filter_summary.to_csv(f"{RES_DIR}/qc_gene_filter.tsv", sep="\t", index=False)
print(f"Saved: {RES_DIR}/qc_gene_filter.tsv")


### --------- Step 3: Scaling & Normalization (log2-CPM) ---------
# Counts Per Million (CPM) normalizes for sequencing depth; log2 stabilizes variance
cpm = counts_filtered.div(library_sizes, axis=1) * 1e6
log_cpm = np.log2(cpm + 1)


### --------- Step 4: Sample-to-Sample Pearson Correlation ---------
correlation = log_cpm.corr()
correlation.to_csv(f"{RES_DIR}/qc_sample_correlation.tsv", sep="\t")
print(f"Saved: {RES_DIR}/qc_sample_correlation.tsv")

# Find lowest off-diagonal correlation
arr = correlation.to_numpy(copy=True)
np.fill_diagonal(arr, np.nan)
min_corr = np.nanmin(arr)
print(f"\nLowest pairwise sample correlation: {min_corr:.4f}")

plt.figure(figsize=(6.5, 5.5))
plt.imshow(correlation, cmap="Blues", vmin=0.90, vmax=1.0)
plt.xticks(range(len(labels)), labels, rotation=30, ha="right", fontsize=10)
plt.yticks(range(len(labels)), labels, fontsize=10)
for i in range(len(labels)):
    for j in range(len(labels)):
        val = correlation.iloc[i, j]
        plt.text(j, i, f"{val:.3f}", ha="center", va="center",
                 color="white" if val > 0.96 else "black", fontsize=9)

plt.title("Sample Correlation Table (log2-CPM)", fontsize=12, fontweight="bold")
plt.colorbar(label="Pearson Correlation")
plt.tight_layout()
plt.savefig(f"{FIG_DIR}/02_sample_correlation.png", dpi=300)
plt.close()
print(f"Saved: {FIG_DIR}/02_sample_correlation.png")


### --------- Step 5: Principal Component Analysis (PCA) ---------
# Select top 500 most variable genes across samples
gene_variance = log_cpm.var(axis=1)
top_500_genes = gene_variance.sort_values(ascending=False).head(500).index
samples_for_pca = log_cpm.loc[top_500_genes].T

pca = PCA(n_components=2, random_state=0)
pca_coords = pca.fit_transform(samples_for_pca)
pc1 = pca_coords[:, 0]
pc2 = pca_coords[:, 1]
var_pc1 = pca.explained_variance_ratio_[0] * 100
var_pc2 = pca.explained_variance_ratio_[1] * 100
print("\n--- PCA variance ---")
print(f"  PC1: {var_pc1:.1f}%")
print(f"  PC2: {var_pc2:.1f}%")

plt.figure(figsize=(7, 5.5))
for i in range(len(labels)):
    plt.scatter(pc1[i], pc2[i], color=colors[i], s=120, edgecolors="black", linewidth=1)
    plt.text(pc1[i] + 0.4, pc2[i] + 0.2, labels[i], fontsize=9)

plt.xlabel(f"PC1: {var_pc1:.1f}% Variance Explained", fontsize=11)
plt.ylabel(f"PC2: {var_pc2:.1f}% Variance Explained", fontsize=11)
plt.title("PCA of NHBE Samples (Top 500 Variable Genes)", fontsize=12, fontweight="bold")
plt.grid(True, linestyle="--", alpha=0.4)
plt.tight_layout()
plt.savefig(f"{FIG_DIR}/03_pca.png", dpi=300)
plt.close()
print(f"Saved: {FIG_DIR}/03_pca.png")

### --------- Step 6: Automated Separation QC Check ---------
mock_pc1 = pc1[meta["condition"] == "mock"]
inf_pc1 = pc1[meta["condition"] == "infected"]

is_separated = (mock_pc1.max() < inf_pc1.min()) or (inf_pc1.max() < mock_pc1.min())

print("\n--- Quality Control Verdict ---")
if is_separated:
    print("PASS: On PC1, the three mock dishes sit apart from the three infected dishes.")
    print("All 6 samples stay in the differential expression test.")
else:
    print("WARNING: Mock and infected dishes overlap on PC1.")
    print("Check figures/03_pca.png before running the gene test.")

print(f"\nQC script finished successfully. Review figures in '{FIG_DIR}/'.")
