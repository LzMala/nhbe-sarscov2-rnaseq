"""
Refits the gene test after dropping the deepest infected dish,
Series1_NHBE_SARS-CoV-2_3. Five dishes remain: 3 mock and 2 infected.
Asks whether the 12 pre-specified gene calls stay the same.
"""

import os
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

# File paths
COUNTS_PATH = "data/processed/nhbe_series1_counts.tsv.gz"
META_PATH = "data/processed/metadata.csv"
PANEL_PATH = "results/gene_panel.tsv"
OUT_PATH = "results/sensitivity_without_deep_library.tsv"

# Same count rule as the main test: at least 10 counts in at least 3 samples.
# Three is still the full mock group. The infected group now has two dishes.
MIN_COUNT = 10
MIN_SAMPLES = 3
PADJ_CUTOFF = 0.05
LFC_CUTOFF = 1.0
DROP_SAMPLE = "Series1_NHBE_SARS-CoV-2_3"

### --------- Step 1: Drop the deepest infected dish ---------
print(f"\nLoading counts from: {COUNTS_PATH}")
counts = pd.read_csv(COUNTS_PATH, sep="\t", index_col=0)
meta = pd.read_csv(META_PATH)
panel = pd.read_csv(PANEL_PATH, sep="\t")

meta = meta[meta["sample_id"] != DROP_SAMPLE].copy()
counts = counts[list(meta["sample_id"])]
print(f"Dropped {DROP_SAMPLE}. Samples left: {list(meta['sample_id'])}")

keep = (counts >= MIN_COUNT).sum(axis=1) >= MIN_SAMPLES
counts_for_model = counts.loc[keep].T.loc[meta["sample_id"]]
if ((counts_for_model % 1) != 0).any().any() or (counts_for_model < 0).any().any():
    print("Error: Counts must be whole numbers, zero or higher.")
    raise SystemExit(1)
counts_for_model = counts_for_model.astype(int)
model_meta = meta.set_index("sample_id")[["condition"]]

### --------- Step 2: Fit the same model on the five remaining dishes ---------
print("\nFitting PyDESeq2 on five samples (3 mock, 2 infected).")
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

stats = DeseqStats(dds, contrast=["condition", "infected", "mock"], n_cpus=1, quiet=True)
stats.summary()
results = stats.results_df.copy()
results.index.name = "gene"

### --------- Step 3: Compare the 12 gene calls ---------
rows = []
for _, gene_row in panel.iterrows():
    gene = gene_row["gene"]
    old_call = gene_row["call"]
    if gene not in results.index:
        # Absent from the file, or still too rare. Keep the original reason.
        if old_call == "not_in_count_table":
            new_call = "not_in_count_table"
        else:
            new_call = "filtered_low_count"
        log2fc = ""
        padj = ""
    else:
        log2fc = float(results.loc[gene, "log2FoldChange"])
        padj = results.loc[gene, "padj"]
        if pd.isna(padj):
            new_call = "no_padj"
        elif padj < PADJ_CUTOFF and log2fc >= LFC_CUTOFF:
            new_call = "up"
        elif padj < PADJ_CUTOFF and log2fc <= -LFC_CUTOFF:
            new_call = "down"
        elif padj < PADJ_CUTOFF:
            new_call = "fdr_only"
        else:
            new_call = "not_significant"
        log2fc = round(log2fc, 3)
    same = "same" if new_call == old_call else "changed"

    rows.append({
        "gene": gene,
        "call_with_all_six": old_call,
        "call_without_deep_library": new_call,
        "log2FC_without_deep_library": log2fc,
        "padj_without_deep_library": padj,
        "comparison": same,
    })

table = pd.DataFrame(rows)
table.to_csv(OUT_PATH, sep="\t", index=False)
print(f"\nSaved: {OUT_PATH}")

print("\n--- Sensitivity Analysis Summary ---")
print(table[["gene", "call_with_all_six", "call_without_deep_library",
             "comparison"]].to_string(index=False))
print("\nSensitivity script finished successfully.")
