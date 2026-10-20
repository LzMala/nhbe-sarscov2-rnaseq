"""
Generates final publication-grade figures for the NHBE SARS-CoV-2 reanalysis:
1. Volcano plot highlighting the pre-specified immune gene panel.
2. Clustered expression heatmap (Row Z-scores) of top DE genes and panel.
3. Horizontal bar chart of top enriched MSigDB Hallmark pathways.
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# File paths
DE_PATH = "results/deseq2_results_full.tsv"
COUNTS_PATH = "data/processed/nhbe_series1_counts.tsv.gz"
META_PATH = "data/processed/metadata.csv"
ORA_PATH = "results/pathway_ora.tsv"
FIG_DIR = "figures"

# Same 12 names as scripts/04_gene_panel.py. Fixed before looking at results.
PANEL = [
    "IFNB1", "IFIT1", "IFIT2", "IFIT3", "ISG15", "MX1", "OAS1", "RSAD2",
    "CXCL8", "IL6", "TNF", "NFKBIA",
]
UP_COLOR = "#D0021B"
DOWN_COLOR = "#4A90E2"
FDR_COLOR = "#E6A817"
NS_COLOR = "#B0B0B0"

os.makedirs(FIG_DIR, exist_ok=True)

print(f"\nLoading DE table from: {DE_PATH}")
de = pd.read_csv(DE_PATH, sep="\t")
counts = pd.read_csv(COUNTS_PATH, sep="\t", index_col=0)
meta = pd.read_csv(META_PATH)

### --------- Step 1: Volcano Plot ---------
plot_de = de.dropna(subset=["log2FC", "padj"]).copy()
plot_de["minus_log10_padj"] = -np.log10(plot_de["padj"].clip(lower=1e-300))

plt.figure(figsize=(7.4, 5.6))

# Gray = not significant. Gold = adjusted p only. Blue = down. Red = up.
for call, color in [("not_significant", NS_COLOR), ("fdr_only", FDR_COLOR), ("down", DOWN_COLOR), ("up", UP_COLOR)]:
    part = plot_de[plot_de["call"] == call]
    plt.scatter(part["log2FC"], part["minus_log10_padj"], s=12, c=color, alpha=0.75, linewidths=0)

# Dashed vertical lines are the 2-fold cuts. The horizontal line is adjusted p = 0.05.
plt.axvline(1, color="gray", linestyle="--", linewidth=0.9, alpha=0.8)
plt.axvline(-1, color="gray", linestyle="--", linewidth=0.9, alpha=0.8)
plt.axhline(-np.log10(0.05), color="gray", linestyle="--", linewidth=0.9, alpha=0.8)

# Label only the pre-specified names that were tested.
label_shift = {
    "IFIT1": (8, 12),
    "IFIT2": (-36, 8),
    "IFIT3": (10, -12),
    "ISG15": (12, 4),
    "RSAD2": (8, 10),
    "TNF": (8, -10),
    "NFKBIA": (-48, 6),
}
labeled = plot_de[plot_de["gene"].isin(PANEL)]
for _, row in labeled.iterrows():
    plt.annotate(
        row["gene"],
        (row["log2FC"], row["minus_log10_padj"]),
        textcoords="offset points",
        xytext=label_shift.get(row["gene"], (4, 4)),
        fontsize=8,
    )

plt.xlabel("log2 Fold Change (Infected vs Mock)", fontsize=11)
plt.ylabel("-log10 Adjusted P-value", fontsize=11)
plt.title("Transcriptomic Response of NHBE to SARS-CoV-2", fontsize=12, fontweight="bold")
plt.grid(True, linestyle="--", alpha=0.4)
plt.tight_layout()
plt.savefig(f"{FIG_DIR}/05_volcano.png", dpi=300)
plt.close()
print(f"Saved: {FIG_DIR}/05_volcano.png")

absent = [gene for gene in PANEL if gene not in counts.index]
too_rare = [gene for gene in PANEL if gene in counts.index and gene not in set(plot_de["gene"])]
if absent:
    print(f"Note: Not in the authors' count table: {', '.join(absent)}")
if too_rare:
    print(f"Note: In the count table, but too rare to test: {', '.join(too_rare)}")

### --------- Step 2: Expression Heatmap (Row Z-scores) ---------
top = de[de["call"].isin(["up", "down"])].nsmallest(20, "padj")
heat_genes = []
for gene in list(top["gene"]) + PANEL:
    if gene in counts.index and gene not in heat_genes and gene in set(de["gene"]):
        heat_genes.append(gene)

sample_order = list(meta["sample_id"])
library_size = counts[sample_order].sum(axis=0)

# log2 counts per million, then each gene is colored against its own average.
log_cpm = np.log2(counts.loc[heat_genes, sample_order].div(library_size, axis=1) * 1e6 + 1)
zscore = log_cpm.sub(log_cpm.mean(axis=1), axis=0).div(log_cpm.std(axis=1).replace(0, np.nan), axis=0)

labels = [name.replace("Series1_NHBE_", "").replace("_", " ") for name in sample_order]

plt.figure(figsize=(8, 8))
image = plt.imshow(zscore.to_numpy(), aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
plt.xticks(range(len(labels)), labels, rotation=30, ha="right")
plt.yticks(range(len(heat_genes)), heat_genes, fontsize=8)
plt.title("Top changed genes and the pre-specified list", fontsize=12, fontweight="bold")
plt.colorbar(image, label="Row z-score of log2 counts per million", fraction=0.03)
plt.tight_layout()
plt.savefig(f"{FIG_DIR}/06_heatmap.png", dpi=300)
plt.close()
print(f"Saved: {FIG_DIR}/06_heatmap.png")

### --------- Step 3: Top Enriched Hallmark Pathways ---------
ora = pd.read_csv(ORA_PATH, sep="\t")
hallmark = ora[ora["Gene_set"] == "MSigDB_Hallmark_2020"].sort_values("Adjusted P-value").head(10)

# A horizontal bar chart draws the first name at the bottom, so flip the strongest to the top.
hallmark = hallmark.iloc[::-1]
scores = -np.log10(hallmark["Adjusted P-value"].clip(lower=1e-300))

plt.figure(figsize=(8, 5.5))
plt.barh(hallmark["Term"], scores, color="#D0021B", edgecolor="black")
plt.xlabel("-log10 adjusted p-value", fontsize=11)
plt.title("Hallmark groups enriched in genes that rose", fontsize=12, fontweight="bold")
plt.grid(axis="x", linestyle="--", alpha=0.4)
plt.tight_layout()
plt.savefig(f"{FIG_DIR}/07_pathways.png", dpi=300)
plt.close()
print(f"Saved: {FIG_DIR}/07_pathways.png")
print("\nFigure script finished successfully. ")
