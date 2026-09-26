"""Read-only dataset preflight; sampled checks are never called full validation."""
import csv
import gzip
from itertools import islice
from pathlib import Path

SNAPSHOTS = ('2022-01', '2023-01', '2024-01', '2025-01', '2026-01', '2026-09')
CLINVAR_COLUMNS = ('VariationID', '#AlleleID', 'Type', 'Name', 'GeneID', 'GeneSymbol',
    'HGNC_ID', 'ClinicalSignificance', 'ClinSigSimple', 'LastEvaluated', 'PhenotypeIDS',
    'PhenotypeList', 'Origin', 'OriginSimple', 'Assembly', 'Chromosome', 'Start', 'Stop',
    'ReferenceAllele', 'AlternateAllele', 'ReviewStatus', 'NumberSubmitters',
    'SubmitterCategories', 'RCVaccession')

def table_check(path, required, delimiter='\t', limit=10000, header_first=None):
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError('Missing or empty file: ' + str(path))
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt', encoding='utf-8-sig', newline='') as stream:
        reader = csv.reader(stream, delimiter=delimiter)
        header = None
        for row in reader:
            if not row: continue
            if header_first:
                if row[0] != header_first: continue
            elif row[0].startswith('#') and row[0] != '#AlleleID':
                continue
            header = row
            break
        if header is None or len(header) != len(set(header)):
            raise ValueError('Missing or duplicate header')
        missing = sorted(set(required) - set(header))
        if missing: raise ValueError('Missing columns: ' + ', '.join(missing))
        count = 0
        missing_fields = dict.fromkeys(required, 0)
        rows = reader if limit is None else islice(reader, limit)
        indices = {name: header.index(name) for name in required}
        for row in rows:
            if not row or all(not v.strip() or set(v.strip()) == {'+'} for v in row): continue
            if len(row) != len(header):
                raise ValueError(f'Malformed data row {count + 1}: {len(row)} fields, expected {len(header)}')
            count += 1
            for name, index in indices.items():
                if row[index].strip().lower() in ('', '-', '.', 'na', 'n/a', 'nan', 'none', 'null'):
                    missing_fields[name] += 1
        if count == 0: raise ValueError('No data records found')
    return {'rows_checked': count, 'scope': 'full file' if limit is None else f'first up to {limit} records only',
            'columns': len(header), 'missing_values_observed': missing_fields}

def ontology_check(path):
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0: raise ValueError('Missing or empty ontology')
    ids, term, name, current = set(), None, False, False
    def finish():
        if current and (term is None or not name): raise ValueError('Term missing ID or name')
    with path.open(encoding='utf-8-sig') as stream:
        if not stream.readline().startswith('format-version:'): raise ValueError('Invalid OBO format header')
        for line in stream:
            line = line.strip()
            if line.startswith('['):
                finish(); current = line == '[Term]'; term = None; name = False
            elif current and line.startswith('id: '):
                term = line[4:]
                if not term.startswith('HP:') or not term[3:].isdigit() or term in ids:
                    raise ValueError('Invalid or duplicate HPO term ID: ' + term)
                ids.add(term)
            elif current and line.startswith('name: '): name = bool(line[6:])
        finish()
    if not ids: raise ValueError('No ontology terms')
    return {'terms_checked': len(ids), 'scope': 'full file structural checks'}

def validate_datasets(root, clinvar_limit=10000):
    root = Path(root)
    results = []
    tasks = [(f'data/raw/clinvar/variant_summary_{s}.txt.gz', lambda p: table_check(p, CLINVAR_COLUMNS, limit=clinvar_limit)) for s in SNAPSHOTS]
    tasks += [
        ('data/raw/hpo/hp.obo', ontology_check),
        ('data/raw/hpo/genes_to_phenotype.txt', lambda p: table_check(p, ('ncbi_gene_id','gene_symbol','hpo_id','hpo_name','frequency','disease_id'), limit=None)),
        ('data/raw/hpo/phenotype.hpoa', lambda p: table_check(p, ('database_id','disease_name','qualifier','hpo_id','reference','evidence','onset','frequency','sex','modifier','aspect','biocuration'), limit=None)),
        ('data/raw/clingen/gene_disease_validity.csv', lambda p: table_check(p, ('GENE SYMBOL','GENE ID (HGNC)','DISEASE LABEL','DISEASE ID (MONDO)','MOI','SOP','CLASSIFICATION','ONLINE REPORT','CLASSIFICATION DATE','GCEP'), delimiter=',', limit=None, header_first='GENE SYMBOL'))]
    for relative, check in tasks:
        try:
            result = {'file': relative, 'status': 'passed', **check(root / relative)}
        except (OSError, EOFError, UnicodeError, ValueError, csv.Error) as error:
            result = {'file': relative, 'status': 'failed', 'error': str(error)}
        results.append(result)
        print(f"{result['status'].upper()}: {relative}", flush=True)
    return {'clinvar_scan_limit': clinvar_limit, 'results': results,
            'passed': sum(r['status'] == 'passed' for r in results), 'expected': len(tasks)}
