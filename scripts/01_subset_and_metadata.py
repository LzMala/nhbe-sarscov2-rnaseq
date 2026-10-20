"""
Subsets the GSE147507 raw counts table to the six Series 1 samples
  (primary NHBE cells: 3 mock vs 3 SARS-CoV-2) and generates metadata.csv.
"""

import os
import pandas as pd

# File paths
COUNT_FILE = "data/raw/GSE147507_RawReadCounts_Human.tsv.gz"
OUT_COUNTS = "data/processed/nhbe_series1_counts.tsv.gz"
OUT_META = "data/processed/metadata.csv"

# Series 1 sample identifiers and annotations (verified from GEO GSE147507)
SAMPLES = [
    "Series1_NHBE_Mock_1",
    "Series1_NHBE_Mock_2",
    "Series1_NHBE_Mock_3",
    "Series1_NHBE_SARS-CoV-2_1",
    "Series1_NHBE_SARS-CoV-2_2",
    "Series1_NHBE_SARS-CoV-2_3",
]
GEO_IDS = [
    "GSM4432378",
    "GSM4432379",
    "GSM4432380",
    "GSM4432381",
    "GSM4432382",
    "GSM4432383",
]
CONDITIONS = ["mock", "mock", "mock", "infected", "infected", "infected"]

### --------- Step 1: Load the full count table ---------
print(f"\nLoading counts from: {COUNT_FILE}")
counts = pd.read_csv(COUNT_FILE, sep="\t", index_col=0)
print(f"Full table shape (genes, samples): {counts.shape}")

# Stop if one of the six Series 1 names is missing. Do not guess a nearby column.
missing = [s for s in SAMPLES if s not in counts.columns]
if missing:
    print(f"\nError: Missing expected samples: {missing}")
    print(f"Available columns: {list(counts.columns)}")
    raise SystemExit(1)

### --------- Step 2: Keep the six Series 1 dishes ---------
subset_counts = counts[SAMPLES]
print(f"\nSubset table shape (genes, samples): {subset_counts.shape}")

print("\n--- Library Sizes (Total reads per sample) ---")
print(subset_counts.sum(axis=0))

zero_genes = (subset_counts.sum(axis=1) == 0).sum()
print(f"\nGenes with zero reads across all 6 samples: {zero_genes}")

### --------- Step 3: Save the six-dish table and the sample sheet ---------
os.makedirs("data/processed", exist_ok=True)
subset_counts.to_csv(OUT_COUNTS, sep="\t", compression="gzip")
print(f"\nSaved: {OUT_COUNTS}")

metadata = pd.DataFrame({
    "sample_id": SAMPLES,
    "geo_accession": GEO_IDS,
    "condition": CONDITIONS,
    "cell_type": ["NHBE"] * 6,
    "series": [1] * 6,
})

metadata.to_csv(OUT_META, index=False)
print(f"Saved: {OUT_META}")
print("\n--- Metadata ---")
print(metadata)
print("\nSubset script finished successfully.")