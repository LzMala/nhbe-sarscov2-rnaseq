"""
Performs pathway enrichment analysis on SARS-CoV-2 infected NHBE cells:
1. Over-Representation Analysis (ORA / Enrichr) on significant upregulated genes.
2. Gene Set Enrichment Analysis (GSEA / PreRank) on all tested genes ranked by Wald stat.
Gene sets: MSigDB Hallmark 2020 and GO Biological Process 2023.
"""

import os
import pandas as pd
import gseapy as gp

# File paths
DE_PATH = "results/deseq2_results_full.tsv"
RES_DIR = "results"
ORA_PATH = f"{RES_DIR}/pathway_ora.tsv"
RANK_PATH = f"{RES_DIR}/pathway_rank.tsv"
GSEA_PATH = f"{RES_DIR}/pathway_gsea.tsv"

# call == "up" already means adjusted p < 0.05 and at least a 2-fold rise,
# assigned in scripts/03_de.py. Changing a cutoff here would not change that label.
GENE_SETS = ["MSigDB_Hallmark_2020", "GO_Biological_Process_2023"]
SEED = 0

os.makedirs(RES_DIR, exist_ok=True)

### --------- Step 1: Load DE Results ---------
print(f"\nLoading DE table from: {DE_PATH}")
de = pd.read_csv(DE_PATH, sep="\t")

up_genes = de.loc[de["call"] == "up", "gene"].dropna().tolist()
background = de["gene"].dropna().tolist()

print(f"Upregulated genes for the list test: {len(up_genes)}")
print(f"Background, every tested gene: {len(background)}")

### --------- Step 2: Over-Representation Analysis (ORA) ---------
# Among the genes that rose, are inflammation or interferon groups
#   more common than they are in the full tested list?
print("\nRunning Over-Representation Analysis (ORA via Enrichr).")
print("This step needs a network connection.")
ora = gp.enrichr(
    gene_list=up_genes,
    gene_sets=GENE_SETS,
    organism="human",
    background=background,
    outdir=None,
    cutoff=1,
)
ora_table = ora.results.copy()
ora_table.to_csv(ORA_PATH, sep="\t", index=False)
print(f"Saved: {ORA_PATH}")


### --------- Step 3: Ranked GSEA (PreRank) ---------
# Rank all genes by Wald statistic (stat = log2FC / lfcSE)
#   A positive statistic means higher after infection.
#   This can see a weak interferon signal that missed the 2-fold line.
ranked = de.dropna(subset=["gene", "stat"]).copy()
ranked["abs_stat"] = ranked["stat"].abs()
ranked = ranked.sort_values("abs_stat", ascending=False).drop_duplicates("gene")
rank = ranked.set_index("gene")["stat"].sort_values(ascending=False)
rank.to_csv(RANK_PATH, sep="\t", header=False)
print(f"Saved: {RANK_PATH}")

print("\nRunning GSEA PreRank (200 shuffles).")
print("This step can take several minutes.")
gsea = gp.prerank(
    rnk=rank,
    gene_sets=GENE_SETS,
    outdir=f"{RES_DIR}/gsea_prerank",
    permutation_num=200,
    seed=SEED,
    threads=1,
    min_size=15,
    max_size=500,
    no_plot=True,
)
gsea_table = gsea.res2d.copy()
gsea_table.to_csv(GSEA_PATH, sep="\t", index=False)
print(f"Saved: {GSEA_PATH}")

### --------- Step 4: Summary of Top Hallmark Pathways ---------
print("\n--- Top Hallmark Pathways (ORA / Enrichr) ---")
hallmark_ora = ora_table[ora_table["Gene_set"] == "MSigDB_Hallmark_2020"].copy()
hallmark_ora = hallmark_ora.sort_values("Adjusted P-value")
print(hallmark_ora[["Term", "Adjusted P-value"]].head(8).to_string(index=False))

print("\n--- Top Pathways (GSEA PreRank) ---")
top_gsea = gsea_table.sort_values("FDR q-val").head(8)
print(top_gsea[["Term", "NES", "FDR q-val"]].to_string(index=False))

print("\nPathway script finished successfully.")