# Decisions

## Samples and counts (Weeks 1–2)

- Used the author file `GSE147507_RawReadCounts_Human.tsv.gz`, not the NCBI GSM-ID matrix.
- Kept only Series 1 NHBE: 3 mock + 3 SARS-CoV-2.
- Dropped A549, Calu-3, Series 9 NHBE, biopsies, and other viruses so the contrast stays one question.
- Library sizes are about 10–16 million for five samples. SARS-CoV-2 replicate 3 is about 29 million. Keep all six. The depth difference is checked again in Week 3.
- 5,437 genes have zero counts in all six samples. Week 3 removes those, and other very low genes, before PCA.

## Week 3 — which genes enter the plots

A gene is kept when it has at least 10 counts in at least 3 samples. Three is the number of replicates in each group, so a gene has to be present across a whole condition, not in a single library. Ten removes genes that are almost all zeros, because their fold changes are unstable. The same rule will be used again before the differential-expression model.

The plots use log2(counts per million + 1). Dividing by library size stops the deepest sample from looking high for every gene. The log makes a doubled gene equal to +1. PCA uses the 500 genes that vary the most across the six samples. The sample labels are not used to build the plot.

The script does not drop a sample. On this run, 12,291 of 21,797 genes pass the filter. PC1 (63.2% of the variation) puts all three mock samples on the left and all three infected samples on the right. SARS-CoV-2 replicate 3 has 29.3 million counts, about 2.7 times the median library, and it sits with the other infected samples. The lowest correlation between two samples is 0.981. All six samples stay in. Mock 3 is higher on PC2 than the other mocks, and it is still on the mock side of PC1, so it is not removed.

Figures from this step: `figures/01_library_sizes.png`, `figures/02_sample_correlation.png`, `figures/03_pca.png`.

The sample sheet columns are `sample_id`, `geo_accession`, `condition`, `cell_type`, and `series`. Later scripts should use those names.

## Week 4 — the differential expression model

PyDESeq2, design `~ condition`. The contrast is infected versus mock, so a positive log2 fold change means higher after infection. The same Week 3 filter is used first: at least 10 counts in at least 3 samples. That leaves 12,291 genes.

The test is a Wald test. Adjusted p-values are Benjamini–Hochberg. With 3 replicates, outlier counts are flagged and not replaced.

Two cutoffs are saved on purpose:

- Primary list: adjusted p < 0.05 and absolute log2 fold change at least 1 (a 2-fold change). 95 genes up, 16 genes down.
- Looser list: adjusted p < 0.05 with no fold-change cut. Another 314 genes pass the adjusted p only. Those stay out of the headline list.

Tables: `results/deseq2_results_full.tsv`, `results/deseq2_results_significant.tsv`, `results/deseq2_results_fdr05.tsv`. The volcano is Week 7.

## Week 5 — the pre-specified gene panel

The author count file already uses gene symbols, so there is no GeneID-to-symbol step.

The list was fixed before looking at the top genes from Week 4:

- Type I interferon and interferon-stimulated genes: IFNB1, IFIT1, IFIT2, IFIT3, ISG15, MX1, OAS1, RSAD2.
- Inflammatory genes: CXCL8, IL6, TNF, NFKBIA.

CXCL8 is not a row in this count table, and neither is the older symbol IL8. It is left as missing. It is not replaced with CXCL1 or another chemokine after seeing the results.

IFNB1 is in the table but fails the Week 3 count filter (about 1 count in mock and below 1 in infected). It is not tested.

Of the interferon-stimulated genes that were tested, MX1 and OAS1 are up (more than 2-fold). IFIT1 and IFIT3 pass the adjusted p-value with a smaller fold change. IFIT2, ISG15, and RSAD2 do not pass. IL6 and TNF are up. NFKBIA passes the adjusted p-value with a fold change just under 2.

That is the result for this slice: IL6 and TNF are clearly induced, and the type I interferon genes are uneven. The list is not expanded to make the interferon story look stronger.

## Week 6 — gene groups

Two tests, both on the Week 4 table.

The list test asks whether the 95 genes that rose at least 2-fold are crowded in known groups, compared with all 12,291 tested genes. The ranked test puts every tested gene in order by the Week 4 statistic, so a gene that missed the 2-fold line can still contribute.

The strongest Hallmark group in the list test is TNF-alpha signaling via NF-kB (adjusted p about 1e-26). Inflammatory response and IL-6/JAK/STAT3 signaling are also near the top. Interferon-gamma response ranks above interferon-alpha response. Interferon-alpha response is enriched (adjusted p about 3e-10) and is not the top group.

## Week 7 — four pictures

The pictures for a reader are `figures/03_pca.png`, `figures/05_volcano.png`, `figures/06_heatmap.png`, and `figures/07_pathways.png`. The volcano labels the pre-specified genes. IFNB1 and CXCL8 are not on it because they were not tested. The heatmap uses the same log2 counts-per-million scale as Week 3, with each gene scaled across the six dishes so a high-count gene does not paint the whole row one color.

## Week 8 — the paper

The NHBE paragraph in Blanco-Melo et al. names CCL20, CXCL1, IL1B, IL6, CXCL3, CXCL5, CXCL6, CXCL2, CXCL16, and TNF as induced, and also IL6, IRF9, ICAM1, and TNF among genes shared with flu. In this table those genes are up, except CXCL16, which passes the adjusted p-value and stays under the 2-fold line.

The paper says type I and type III interferon are undetectable in these cells. IFNB1, IFNL1, IFNL2, and IFNL3 are in the count file and are too rare to test. That agrees. The comparison is in `results/paper_overlap.tsv`. Differences that remain can come from the genome build, the filters, and the software. The authors' full supplemental table was not copied into this repository.

## Weeks 11–12 — one extra check

The model was fit again without Series1_NHBE_SARS-CoV-2_3, the deepest infected dish. IL6, TNF, MX1, and OAS1 stay up. NFKBIA stays under the 2-fold line. IFIT1 and IFIT3 lose the adjusted-p pass. The headline does not depend on that dish. The borderline interferon genes do. With only two infected dishes left, the uncertainty estimate is shaky, so this check is a note, not a second main result. Table: `results/sensitivity_without_deep_library.tsv`.
