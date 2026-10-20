"""
Benchmarks our Series 1 reanalysis against the qualitative claims
made by Blanco-Melo et al. (Cell 2020) for primary NHBE cells:
1. Validates paper-highlighted chemokines and cytokines.
2. Validates lack of Type-I / Type-III interferon induction.
3. Compares our top 20 differentially expressed genes.
"""

import os
import pandas as pd

# File paths
DE_PATH = "results/deseq2_results_full.tsv"
COUNTS_PATH = "data/processed/nhbe_series1_counts.tsv.gz"
SIG_PATH = "results/deseq2_results_significant.tsv"
OUT_PATH = "results/paper_overlap.tsv"

# Genes explicitly named in the NHBE text of Blanco-Melo et al., Cell 2020.
PAPER_UP = [
    "CCL20", "CXCL1", "IL1B", "IL6", "CXCL3", "CXCL5",
    "CXCL6", "CXCL2", "CXCL16", "TNF", "IRF9", "ICAM1",
]
PAPER_QUIET = ["IFNB1", "IFNL1", "IFNL2", "IFNL3"]

os.makedirs("results", exist_ok=True)

### --------- Step 1: Score the genes named in the NHBE paragraph ---------
print(f"\nLoading DE table from: {DE_PATH}")
de = pd.read_csv(DE_PATH, sep="\t").set_index("gene")
counts = pd.read_csv(COUNTS_PATH, sep="\t", index_col=0)

rows = []
for gene in PAPER_UP + PAPER_QUIET:
    if gene in PAPER_UP:
        paper_claim = "paper reports this gene induced in NHBE"
    else:
        paper_claim = "paper reports type I or III interferon undetectable in NHBE"

    if gene not in de.index:
        if gene in counts.index:
            our_call = "filtered_low_count"
            agreement = "agree, too rare to test"
        else:
            our_call = "not_in_count_table"
            agreement = "cannot score"
        rows.append({
            "gene": gene,
            "paper_claim": paper_claim,
            "our_call": our_call,
            "log2FC": "",
            "padj": "",
            "agreement": agreement,
        })
        continue

    call = de.loc[gene, "call"]
    log2fc = de.loc[gene, "log2FC"]
    padj = de.loc[gene, "padj"]

    if gene in PAPER_UP and call == "up":
        agreement = "agree, up"
    elif gene in PAPER_UP and call == "fdr_only":
        agreement = "same direction, below the 2-fold line"
    elif gene in PAPER_UP and call == "down":
        agreement = "disagree, we call it down"
    elif gene in PAPER_UP and call == "not_significant":
        agreement = "disagree, we do not call it changed"
    elif gene in PAPER_QUIET and call == "up":
        agreement = "disagree, we call it up"
    elif gene in PAPER_QUIET:
        agreement = "agree, not a clear induced interferon"
    else:
        agreement = "partial"

    rows.append({
        "gene": gene,
        "paper_claim": paper_claim,
        "our_call": call,
        "log2FC": round(float(log2fc), 3),
        "padj": padj,
        "agreement": agreement,
    })

### --------- Step 2: Add our top 20 genes the paragraph did not name ---------
top20 = pd.read_csv(SIG_PATH, sep="\t").head(20)
paper_names = set(PAPER_UP + PAPER_QUIET)
for gene in top20["gene"]:
    if gene in paper_names:
        continue
    row = top20.loc[top20["gene"] == gene].iloc[0]
    rows.append({
        "gene": gene,
        "paper_claim": "in our top 20, not one of the names in the NHBE paragraph",
        "our_call": row["call"],
        "log2FC": round(float(row["log2FC"]), 3),
        "padj": row["padj"],
        "agreement": "our result, paper did not highlight this name in the NHBE text",
    })

table = pd.DataFrame(rows)
table.to_csv(OUT_PATH, sep="\t", index=False)
print(f"\nSaved: {OUT_PATH}")

print(table[["gene", "paper_claim", "our_call", "log2FC", "agreement"]].to_string(index=False))
print("\nPaper comparison script finished successfully.")
