"""Download GSE147507 human raw counts. Run from the project root."""

from datetime import datetime, timezone
from pathlib import Path
import hashlib
import urllib.request

URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE147nnn/"
    "GSE147507/suppl/GSE147507_RawReadCounts_Human.tsv.gz"
)
OUT = Path("data/raw/GSE147507_RawReadCounts_Human.tsv.gz")
LOG = Path("data/raw/DOWNLOAD_LOG.txt")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    OUT.parent.mkdir(parents=True, exist_ok=True)
    print("Downloading:", URL)
    urllib.request.urlretrieve(URL, OUT)
    digest = sha256(OUT)
    size = OUT.stat().st_size
    now = datetime.now(timezone.utc).isoformat()
    LOG.write_text(
        f"url={URL}\n"
        f"downloaded_utc={now}\n"
        f"path={OUT}\n"
        f"bytes={size}\n"
        f"sha256={digest}\n",
        encoding="utf-8",
    )
    print("Saved", OUT, "bytes=", size)
    print("sha256", digest)
    print("Log written to", LOG)


if __name__ == "__main__":
    main()