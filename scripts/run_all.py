"""
Runs the analysis from the download through the five-sample check.
Start in the project folder, not inside scripts.
Script 08 refits the model after dropping the deepest infected dish.
"""

import subprocess
import sys

SCRIPTS = [
    "scripts/00_download.py",
    "scripts/01_subset_and_metadata.py",
    "scripts/02_qc.py",
    "scripts/03_de.py",
    "scripts/04_gene_panel.py",
    "scripts/05_pathways.py",
    "scripts/06_figures.py",
    "scripts/07_paper_comparison.py",
    "scripts/08_sensitivity.py",
]

for script in SCRIPTS:
    print(f"\n========== {script} ==========")
    subprocess.run([sys.executable, script], check=True)

print("\nFull run finished.")
