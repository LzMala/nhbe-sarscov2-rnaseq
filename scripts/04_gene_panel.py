"""
Evaluates the pre-specified biological hypothesis:
    Comparing Type-I Interferon / ISGs against pro-inflammatory genes
    in primary NHBE cells infected with SARS-CoV-2.
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# Fixed before looking at the results. scripts/06_figures.py repeats this list.
TYPE_I_ISG = ["IFNB1", "IFIT1", "IFIT2", "IFIT3", "ISG15", "MX1", "OAS1", "RSAD2"]
INFLAMMATORY = ["CXCL8", "IL6", "TNF", "NFKBIA"]

# File paths
COUNTS_PATH = "data/processed/nhbe_series1_counts.tsv.gz"
DE_PATH = "results/deseq2_results_full.tsv"
META_PATH = "data/processed/metadata.csv"
RES_DIR = "results"
FIG_DIR = "figures"

PANEL_PATH = f"{RES_DIR}/gene_panel.tsv"
SUMMARY_PATH = f"{RES_DIR}/gene_panel_summary.tsv"
FIG_PATH = f"{FIG_DIR}/04_gene_panel.png"

# On this figure, blue is the interferon list and red is the inflammation list.
ISG_COLOR = "#4A90E2"
INFLAM_COLOR = "#D0021B"

os.makedirs(RES_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)

### --------- Step 1: Load Data ---------
print(f"\nLoading DE table from: {DE_PATH}")
print(f"Loading counts from: {COUNTS_PATH}")
de = pd.read_csv(DE_PATH, sep="\t").set_index("gene")
counts = pd.read_csv(COUNTS_PATH, sep="\t", index_col=0)
meta = pd.read_csv(META_PATH)

mock_samples = meta.loc[meta["condition"] == "mock", "sample_id"]
infected_samples = meta.loc[meta["condition"] == "infected", "sample_id"]

panel_genes = TYPE_I_ISG + INFLAMMATORY
groups = ["type_I_isg"] * len(TYPE_I_ISG) + ["inflammatory"] * len(INFLAMMATORY)


### --------- Step 2: Score Each Pre-Specified Gene ---------
rows = []
for gene, group in zip(panel_genes, groups):
    in_counts = gene in counts.index
    in_model = gene in de.index

    if in_counts:
        mean_mock = counts.loc[gene, mock_samples].mean()
        mean_infected = counts.loc[gene, infected_samples].mean()
    else:
        mean_mock = float("nan")
        mean_infected = float("nan")

    if not in_counts:
        call = "not_in_count_table"
        log2fc = float("nan")
        padj = float("nan")
        base_mean = float("nan")
    elif not in_model:
        # Present, but it failed rule of 10 counts in 3 samples.
        call = "filtered_low_count"
        log2fc = float("nan")
        padj = float("nan")
        base_mean = float("nan")
    else:
        call = de.loc[gene, "call"]
        log2fc = de.loc[gene, "log2FC"]
        padj = de.loc[gene, "padj"]
        base_mean = de.loc[gene, "baseMean"]

    rows.append({
        "gene": gene,
        "group": group,
        "in_count_table": in_counts,
        "in_model": in_model,
        "call": call,
        "log2FC": log2fc,
        "padj": padj,
        "baseMean": base_mean,
        "mean_mock": mean_mock,
        "mean_infected": mean_infected,
    })

panel = pd.DataFrame(rows)
panel.to_csv(PANEL_PATH, sep="\t", index=False)
print(f"\nSaved: {PANEL_PATH}")

print("\n--- Pre-specified Gene Panel Results ---")
print(panel[["gene", "group", "call", "log2FC",
             "padj","mean_mock", "mean_infected"]].to_string(index=False))

### --------- Step 3: Summary by gene group ---------
summary = pd.crosstab(panel["group"], panel["call"])
summary.to_csv(SUMMARY_PATH, sep="\t")
print(f"\nSaved: {SUMMARY_PATH}")

print("\n--- Summary of Calls by Group ---")
print(summary)

### --------- Step 4: Gene panel bar chart ---------
# Only genes that entered the test have a fold change to draw.
tested = panel[panel["in_model"]].copy()
tested = tested.sort_values(by="log2FC", ascending=True)

bar_colors = [ISG_COLOR if grp == "type_I_isg" else INFLAM_COLOR for grp in tested["group"]]

plt.figure(figsize=(7.5, 5.5))
plt.barh(tested["gene"], tested["log2FC"], color=bar_colors, edgecolor="black", alpha=0.9)
plt.axvline(0, color="black", linewidth=0.8)
plt.axvline(1, color="gray", linestyle="--", linewidth=1)
plt.axvline(-1, color="gray", linestyle="--", linewidth=1)

plt.xlabel("log2 Fold Change (Infected vs Mock)", fontsize=11)
plt.title("Response of Pre-Specified Immune Genes (SARS-CoV-2 in NHBE)", fontsize=12, fontweight="bold")
plt.legend(
    handles=[
        Patch(facecolor=ISG_COLOR, edgecolor="black", label="Type-I / ISGs (Antiviral)"),
        Patch(facecolor=INFLAM_COLOR, edgecolor="black", label="Inflammatory (NF-κB)"),
    ],
    frameon=False,
    loc="lower right"
)
plt.grid(axis="x", linestyle="--", alpha=0.5)
plt.tight_layout()
plt.savefig(FIG_PATH, dpi=300)
plt.close()
print(f"\nSaved: {FIG_PATH}")

absent = panel.loc[panel["call"] == "not_in_count_table", "gene"].tolist()
too_rare = panel.loc[panel["call"] == "filtered_low_count", "gene"].tolist()
if absent:
    print(f"\nNote: Not in the authors' count table: {', '.join(absent)}")
if too_rare:
    print(f"Note: In the count table, but too rare to test: {', '.join(too_rare)}")
print("\nGene panel script finished successfully.")

