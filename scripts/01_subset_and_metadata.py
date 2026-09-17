"""Keep Series1 NHBE samples only and write metadata.csv."""

from pathlib import Path
import pandas as pd

RAW = Path("data/raw/GSE147507_RawReadCounts_Human.tsv.gz")
OUT_COUNTS = Path("data/processed/nhbe_series1_counts.tsv.gz")
OUT_META = Path("data/processed/metadata.csv")

KEEP = [
    "Series1_NHBE_Mock_1",
    "Series1_NHBE_Mock_2",
    "Series1_NHBE_Mock_3",
    "Series1_NHBE_SARS-CoV-2_1",
    "Series1_NHBE_SARS-CoV-2_2",
    "Series1_NHBE_SARS-CoV-2_3",
]

META_ROWS = [
    {"sample": "Series1_NHBE_Mock_1", "geo": "GSM4432378", "condition": "mock", "cell_type": "NHBE", "series": 1},
    {"sample": "Series1_NHBE_Mock_2", "geo": "GSM4432379", "condition": "mock", "cell_type": "NHBE", "series": 1},
    {"sample": "Series1_NHBE_Mock_3", "geo": "GSM4432380", "condition": "mock", "cell_type": "NHBE", "series": 1},
    {"sample": "Series1_NHBE_SARS-CoV-2_1", "geo": "GSM4432381", "condition": "infected", "cell_type": "NHBE", "series": 1},
    {"sample": "Series1_NHBE_SARS-CoV-2_2", "geo": "GSM4432382", "condition": "infected", "cell_type": "NHBE", "series": 1},
    {"sample": "Series1_NHBE_SARS-CoV-2_3", "geo": "GSM4432383", "condition": "infected", "cell_type": "NHBE", "series": 1},
]


def main() -> None:
    OUT_COUNTS.parent.mkdir(parents=True, exist_ok=True)

    counts = pd.read_csv(RAW, sep="\t", index_col=0)
    print("Full table shape (genes x samples):", counts.shape)
    print("Column names in the file:")
    for name in counts.columns:
        print(" ", name)

    missing = [c for c in KEEP if c not in counts.columns]
    if missing:
        raise SystemExit(
            "These sample names were not in the count file:\n"
            + "\n".join(missing)
            + "\nLook at the column list above and stop."
        )

    small = counts[KEEP]
    print("Subset shape:", small.shape)
    print("Library sizes (total counts per sample):")
    print(small.sum(axis=0).to_string())

    zeros = int((small.sum(axis=1) == 0).sum())
    print("Genes with zero counts in all 6 samples:", zeros)

    small.to_csv(OUT_COUNTS, sep="\t", compression="gzip")
    pd.DataFrame(META_ROWS).to_csv(OUT_META, index=False)
    print("Wrote", OUT_COUNTS)
    print("Wrote", OUT_META)


if __name__ == "__main__":
    main()