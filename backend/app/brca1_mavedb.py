"""Optional BRCA1 MaveDB evidence package loader."""
from functools import lru_cache
import gzip
import json
import os
from pathlib import Path

from src.brca1_mavedb import assess_variant

DEFAULT_PATH = Path(__file__).resolve().parents[2] / "data/processed/brca1_mavedb_v1.json.gz"

@lru_cache(maxsize=2)
def load_package(path, modified_ns):
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return json.load(handle)

def evidence(record):
    if not record or record.get("gene_symbol") != "BRCA1":
        return {"status": "out_of_scope", "reason": "This evidence pilot covers BRCA1 only."}
    path = Path(os.environ.get("GENEVISTA_BRCA1_MAVEDB", DEFAULT_PATH))
    try:
        package = load_package(str(path), path.stat().st_mtime_ns)
        return assess_variant(record, package)
    except FileNotFoundError:
        return {"status": "not_configured", "reason": "The BRCA1 MaveDB evidence package is not installed."}
    except (OSError, EOFError, ValueError, KeyError, TypeError):
        return {"status": "unavailable", "reason": "The BRCA1 MaveDB evidence package needs rebuilding."}
