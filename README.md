# NHBE SARS-CoV-2 bulk RNA-seq reanalysis

In primary human bronchial epithelial cells, which genes change 24 hours after SARS-CoV-2 compared with mock, and is the type I interferon response actually on?

This is a reanalysis of public data. I did not generate these samples.

## Data

GEO series [GSE147507](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE147507), Series 1 only.

| Sample | GEO | Group |
| --- | --- | --- |
| Series1_NHBE_Mock_1 | GSM4432378 | mock |
| Series1_NHBE_Mock_2 | GSM4432379 | mock |
| Series1_NHBE_Mock_3 | GSM4432380 | mock |
| Series1_NHBE_SARS-CoV-2_1 | GSM4432381 | infected |
| Series1_NHBE_SARS-CoV-2_2 | GSM4432382 | infected |
| Series1_NHBE_SARS-CoV-2_3 | GSM4432383 | infected |

Paper: Blanco-Melo et al., Cell 2020 (PMID 32416070).

Counts: the author file `GSE147507_RawReadCounts_Human.tsv.gz`. Six columns only. The other cell types, flu, and the later NHBE experiment stay out.

## What this reanalysis did

1. Downloaded the author count file and recorded the link, date, and fingerprint.
2. Kept the six Series 1 samples and wrote `data/processed/metadata.csv`.
3. Checked depth, sample similarity, and a PCA. A gene entered later steps only with at least 10 counts in at least 3 samples (12,291 genes).
4. Tested infected versus mock with PyDESeq2. Mock is the reference. A positive log2 fold change means higher after infection. Primary cutoff: adjusted p below 0.05 and at least a 2-fold change.
5. Scored 12 genes chosen before looking at the top of the list.
6. Tested whether known gene groups, especially NF-kB, inflammation, and interferon, are enriched.
7. Compared the genes named for NHBE in the paper with this table.

## Main result

95 genes went up and 16 went down at adjusted p below 0.05 and at least a 2-fold change. IL6 and TNF are clearly up. IFN-beta (IFNB1) is too rare to test. MX1 and OAS1 are up, and several other interferon-stimulated genes are not. The strongest gene group among genes that rose is TNF-alpha signaling through NF-kB. Interferon-alpha response is enriched, and it is not the top group. This agrees with the paper's NHBE claim that inflammation is the louder program. It does not show that the response is weak compared with flu, because those samples were not included.

## How to reproduce

Start in this folder, not inside `scripts`.

```powershell
cd C:\Users\lzmah\GROK\nhbe-sarscov2-rnaseq
conda activate nhbe-sarscov2
python scripts/run_all.py
```

A run worked when the terminal prints lines that start with `Saved`. The same command ends with the five-dish check. That check can also be repeated alone with `python scripts/08_sensitivity.py`.

The four pictures to look at:

- `figures/03_pca.png` — mock versus infected
- `figures/05_volcano.png` — the pre-specified genes on the full test
- `figures/06_heatmap.png` — top changed genes plus that list
- `figures/07_pathways.png` — strongest Hallmark groups

A plain-language walkthrough is in `docs/NHBE_case_study.txt`. Decisions are in `docs/decisions.md`.

## Limitations

See `LIMITATIONS.md`. Short version: n = 3, cells in a dish, author-made counts, no flu comparison, no new biomarkers.

## How to cite the original data

See `CITATION.md`. Cite Blanco-Melo et al., Cell 2020, and GEO GSE147507. This repository is a reanalysis of that public data.
