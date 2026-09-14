"""Small, optional TP53 evidence package for local and hosted workspaces."""
from functools import lru_cache
import gzip
import json
import logging
import os
from pathlib import Path

from src.tp53_evidence import assess_variant

DEFAULT_PATH = Path(__file__).resolve().parents[2] / 'data/processed/tp53_evidence_v1.json.gz'


@lru_cache(maxsize=2)
def load_package(path, modified_ns):
    # The modification key allows an atomically rebuilt artifact to be reloaded.
    with gzip.open(path, 'rt', encoding='utf-8') as handle:
        return json.load(handle)


def evidence(record):
    if not record or record.get('gene_symbol') != 'TP53':
        return {'status': 'out_of_scope', 'reason': 'This evidence pilot covers TP53 only.'}
    path = Path(os.environ.get('GENEVISTA_TP53_EVIDENCE', DEFAULT_PATH))
    try:
        package = load_package(str(path), path.stat().st_mtime_ns)
        return assess_variant(record, package)
    except FileNotFoundError:
        return {'status': 'not_configured', 'reason': 'The TP53 evidence package is not installed.'}
    except (OSError, EOFError, ValueError, KeyError, TypeError):
        logging.exception('TP53 evidence package could not be read')
        return {'status': 'unavailable', 'reason': 'The TP53 evidence package needs rebuilding.'}
