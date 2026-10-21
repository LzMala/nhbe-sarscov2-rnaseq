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

## Reproduce and run

The commands below use the Miniconda environment named in `environment.yaml`. Start in this project folder. The scripts look for `data\`, `figures\`, and `results\` from here. Starting inside `scripts` writes a second copy of those folders.

Create the environment once:

```powershell
conda env create -f environment.yaml
conda activate nhbe-sarscov2
```

If the environment is already on the computer, skip the create line and activate it:

```powershell
conda activate nhbe-sarscov2
python scripts/run_all.py
```

That one command runs the scripts in this order:

1. `scripts/00_download.py` downloads the authors' count table. This step needs an internet connection.
2. `scripts/01_subset_and_metadata.py` keeps the six Series 1 dishes and writes the sample sheet.
3. `scripts/02_qc.py` checks sequencing depth, whether the dishes agree, and whether mock sits apart from infected.
4. `scripts/03_de.py` tests each gene. This step takes a few minutes.
5. `scripts/04_gene_panel.py` scores the 12 genes chosen before looking at the winners.
6. `scripts/05_pathways.py` tests the named gene groups. This step needs an internet connection and can take several minutes.
7. `scripts/06_figures.py` draws the volcano plot, the heatmap, and the pathway bars.
8. `scripts/07_paper_comparison.py` scores the gene names from the paper's NHBE paragraph.
9. `scripts/08_sensitivity.py` repeats the gene test without the deepest infected dish. This step takes a few minutes. It is a check, not a second main result.

A run worked when each script prints lines that start with `Saved` and ends with `finished successfully`.

To repeat one step, activate the same environment and run that file alone. Example:

```powershell
python scripts/04_gene_panel.py
```

The new files land in `data\processed\`, `figures\`, and `results\`.

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
