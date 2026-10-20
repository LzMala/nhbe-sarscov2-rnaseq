"""
Downloads the authors' human count table for GSE147507 from NCBI GEO
and writes a receipt with the time, the file size, and a SHA-256 fingerprint.
"""

from datetime import datetime, timezone
import hashlib
import os
import urllib.request

# Authors' supplementary count table, Blanco-Melo et al., Cell 2020.
URL = ( "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE147nnn/"
        "GSE147507/suppl/GSE147507_RawReadCounts_Human.tsv.gz")

# File paths
OUT = "data/raw/GSE147507_RawReadCounts_Human.tsv.gz"
LOG = "data/raw/DOWNLOAD_LOG.txt"

### --------- Step 1: Download the count table ---------
os.makedirs("data/raw", exist_ok=True)
print(f"\nDownloading GSE147507 counts from: {URL}")
urllib.request.urlretrieve(URL, OUT)

### --------- Step 2: Write the download receipt ---------
# The fingerprint changes if any byte in the file changes.
with open(OUT, "rb") as f:
    file_bytes = f.read()

fingerprint = hashlib.sha256(file_bytes).hexdigest()
size = len(file_bytes)
timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")

log_content = (
    f"url={URL}\n"
    f"downloaded on (UTC)={timestamp}\n"
    f"path={OUT}\n"
    f"bytes={size}\n"
    f"sha256={fingerprint}\n")

with open(LOG, "w", encoding="utf-8") as f:
    f.write(log_content)

print(f"Saved: {OUT}")
print(f"Size (bytes): {size}")
print(f"Sha256: {fingerprint}")
print(f"Saved: {LOG}")
print("\nDownload script finished successfully.")
