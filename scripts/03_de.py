"""
Tests which genes change in primary NHBE cells after SARS-CoV-2 versus mock.
PyDESeq2, design ~ condition, mock as the reference. A positive log2 fold
change means higher after infection.
"""

import os
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

MIN_COUNT = 10
MIN_SAMPLES = 3

# Primary call: adjusted p < 0.05 and at least a 2-fold change.
PADJ_CUTOFF = 0.05
LFC_CUTOFF = 1.0

# File paths
COUNTS_PATH = "data/processed/nhbe_series1_counts.tsv.gz"
META_PATH = "data/processed/metadata.csv"
RES_DIR = "results"

FULL_PATH = f"{RES_DIR}/deseq2_results_full.tsv"
SIG_PATH = f"{RES_DIR}/deseq2_results_significant.tsv"
LOOSE_PATH = f"{RES_DIR}/deseq2_results_fdr05.tsv"
SUMMARY_PATH = f"{RES_DIR}/deseq2_summary.tsv"

os.makedirs(RES_DIR, exist_ok=True)

### --------- Step 1: Load the six samples ---------
print(f"\nLoading counts from: {COUNTS_PATH}")
print(f"Loading metadata from: {META_PATH}")
counts = pd.read_csv(COUNTS_PATH, sep="\t", index_col=0)
meta = pd.read_csv(META_PATH)

if list(counts.columns) != list(meta["sample_id"]):
    print("Error: Count table columns do not match metadata sample_id order.")
    raise SystemExit(1)

if counts.index.duplicated().any():
    print("Error: The count table has duplicate gene names.")
    raise SystemExit(1)

### --------- Step 2: Filter Lowly Expressed Genes ---------
keep_mask = (counts >= MIN_COUNT).sum(axis=1) >= MIN_SAMPLES
counts_filtered = counts.loc[keep_mask]
print("\n--- Genes entering the model ---")
print(f"  Genes in the six-sample table: {counts.shape[0]:,}")
print(f"  Genes tested (>= {MIN_COUNT} counts in >= {MIN_SAMPLES} samples): {counts_filtered.shape[0]:,}")

# PyDESeq2 wants one row per sample and whole-number counts.
counts_for_model = counts_filtered.T.loc[meta["sample_id"]]
if ((counts_for_model % 1) != 0).any().any() or (counts_for_model < 0).any().any():
    print("Error: Counts must be whole numbers, zero or higher.")
    raise SystemExit(1)
counts_for_model = counts_for_model.astype(int)

model_meta = meta.set_index("sample_id")[["condition"]]
print(model_meta)


### --------- Step 3: Fit the model ---------
# Mock is the baseline. Infected is the group we compare to it.
# With only 3 dishes per group, a wild count is flagged and left as it is.
print("\nFitting PyDESeq2.")
print("This step can take a few minutes.")
dds = DeseqDataSet(
    counts=counts_for_model,
    metadata=model_meta,
    design="~condition",
    refit_cooks=False,
    n_cpus=1,
    quiet=True,
)
dds.deseq2()

### --------- Step 4: Wald test, infected versus mock ---------
stats = DeseqStats(
    dds,
    contrast=["condition", "infected", "mock"],
    n_cpus=1,
    quiet=True,
)
stats.summary()
results = stats.results_df.copy()
results.index.name = "gene"

results = results.rename(columns={"log2FoldChange": "log2FC"})
results = results[["baseMean", "log2FC", "lfcSE", "stat", "pvalue", "padj"]]
results = results.sort_values(["padj", "pvalue"], ascending=True)

# up / down = passes both cuts. fdr_only = passes the adjusted p, but the fold change is smaller than 2.
results["call"] = "not_significant"
results.loc[results["padj"] < PADJ_CUTOFF, "call"] = "fdr_only"
results.loc[(results["padj"] < PADJ_CUTOFF) & (results["log2FC"] >= LFC_CUTOFF), "call"] = "up"
results.loc[(results["padj"] < PADJ_CUTOFF) & (results["log2FC"] <= -LFC_CUTOFF), "call"] = "down"
results.loc[results["padj"].isna(), "call"] = "no_padj"

### --------- Step 5: Save the tables ---------
results.to_csv(FULL_PATH, sep="\t")
print(f"\nSaved: {FULL_PATH}")

significant = results[results["call"].isin(["up", "down"])].copy()
significant.to_csv(SIG_PATH, sep="\t")
print(f"Saved: {SIG_PATH}")


# Looser table: adjusted p < 0.05, no fold-change cut. Used to see how the cutoff changes the list.
fdr_only_table = results[results["padj"] < PADJ_CUTOFF].copy()
fdr_only_table.to_csv(LOOSE_PATH, sep="\t")
print(f"Saved: {LOOSE_PATH}")

n_up = int((results["call"] == "up").sum())
n_down = int((results["call"] == "down").sum())
n_fdr = int((results["call"] == "fdr_only").sum())
n_tested = int(results.shape[0])

summary = pd.DataFrame({
    "metric": [
        "genes_tested",
        "padj_cutoff",
        "log2fc_cutoff",
        "genes_up",
        "genes_down",
        "genes_fdr_only",
        "genes_significant_both_cuts",
    ],
    "value": [
        n_tested,
        PADJ_CUTOFF,
        LFC_CUTOFF,
        n_up,
        n_down,
        n_fdr,
        n_up + n_down,
    ],
})
summary.to_csv(SUMMARY_PATH, sep="\t", index=False)
print(f"Saved: {SUMMARY_PATH}")

print("\n--- Differential expression summary ---")
print(f"  Contrast: infected versus mock")
print(f"  Genes tested: {n_tested:,}")
print(f"  Up  (padj < {PADJ_CUTOFF} and log2FC >= {LFC_CUTOFF}): {n_up:,}")
print(f"  Down (padj < {PADJ_CUTOFF} and log2FC <= -{LFC_CUTOFF}): {n_down:,}")
print(f"  Adjusted p < {PADJ_CUTOFF} but |log2FC| < {LFC_CUTOFF}: {n_fdr:,}")

print("\nTop genes by adjusted p:")
print(significant.head(10)[["baseMean", "log2FC", "padj", "call"]])

print("\nDE script finished successfully.")
