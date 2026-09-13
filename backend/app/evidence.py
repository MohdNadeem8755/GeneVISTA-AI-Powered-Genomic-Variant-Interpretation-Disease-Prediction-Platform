"""Prepared gene-level associations; these are not variant or patient predictions."""
import csv
import gzip
import re
from functools import lru_cache
from pathlib import Path

SOURCE = Path(__file__).resolve().parents[2] / 'data/processed/clingen_gene_evidence_mapped.csv.gz'
ORDER = {'Definitive': 0, 'Strong': 1, 'Moderate': 2, 'Limited': 3,
         'Disputed': 4, 'Refuted': 5, 'No Known Disease Relationship': 6,
         'No Reported Evidence': 6}

@lru_cache(maxsize=2)
def _load(stamp):
    index = {}
    with gzip.open(SOURCE, 'rt', encoding='utf-8-sig', newline='') as handle:
        for row in csv.DictReader(handle):
            index.setdefault(row['gene_symbol'].upper(), []).append(row)
    return index

def disease_evidence(genes):
    if not SOURCE.exists():
        return {'status': 'unavailable', 'items': [], 'reason': 'Prepared ClinGen evidence file is unavailable.'}
    index = _load(SOURCE.stat().st_mtime_ns)
    rows = [dict(row) for gene in sorted(set(genes)) for row in index.get(gene.upper(), [])]
    rows.sort(key=lambda row: (ORDER.get(row['classification'], 7), row['disease_name'], row['gene_symbol'], row['mode_of_inheritance']))
    for row in rows:
        row['probability'] = None
        row['confidence_score'] = None
        row['evidence_scope'] = 'gene_disease'
        # Only expose links to the trusted curation host.
        if not row.get('report_url', '').startswith('https://search.clinicalgenome.org/'):
            row['report_url'] = None
    return {'status': 'available', 'items': rows, 'source': 'ClinGen prepared gene-disease curations',
            'ordering': 'Evidence strength, then disease name; not likelihood of disease',
            'probabilities_supported': False}
