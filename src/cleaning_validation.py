"""Helpers for notebook 02; exact targets and grouped historical states differ."""
from pathlib import Path
import pandas as pd

EXACT_CLASSES = ('Pathogenic', 'Likely pathogenic', 'Uncertain significance', 'Likely benign', 'Benign')

def exact_class(value):
    if pd.isna(value): return pd.NA
    return {x.lower(): x for x in EXACT_CLASSES}.get(str(value).strip().lower(), pd.NA)

def parsed_dates(values):
    result = pd.to_datetime(values, format='mixed', errors='coerce', utc=True)
    invalid = values.notna() & result.isna()
    if invalid.any():
        raise ValueError('Unrecognized classification dates: ' + str(values[invalid].unique().tolist()))
    return result

def source_date(path):
    import csv
    with Path(path).open(encoding='utf-8-sig', newline='') as stream:
        for row in csv.reader(stream):
            if row and row[0].startswith('FILE CREATED:'):
                return row[0].split(':', 1)[1].strip()
    raise ValueError('ClinGen source creation date missing')

def export_exact_classes(source, destination, chunksize=100000):
    """Preserve every retained record; mixed assertions never become exact labels."""
    import gzip
    source, destination = Path(source), Path(destination)
    temporary = destination.with_suffix('.partial.gz')
    count, eligible = 0, 0
    try:
        with gzip.open(temporary, 'wt', encoding='utf-8', newline='') as handle:
            for chunk in pd.read_csv(source, dtype='string', keep_default_na=False,
                    usecols=['VariationID','snapshot','ClinicalSignificance'], chunksize=chunksize):
                chunk['exact_clinical_class'] = chunk['ClinicalSignificance'].map(exact_class).astype('string')
                chunk['five_class_eligible'] = chunk['exact_clinical_class'].notna()
                chunk.to_csv(handle,index=False,header=count==0)
                count += len(chunk); eligible += int(chunk['five_class_eligible'].sum())
        if count == 0: raise ValueError('No retained records for exact-class export')
        from itertools import zip_longest
        checked = 0
        original = pd.read_csv(source, dtype='string', keep_default_na=False,
            usecols=['VariationID','snapshot','ClinicalSignificance'], chunksize=chunksize)
        saved = pd.read_csv(temporary,dtype='string',keep_default_na=False,chunksize=chunksize)
        for before, after in zip_longest(original, saved):
            if before is None or after is None: raise ValueError('Exact-class export length mismatch')
            pd.testing.assert_frame_equal(before.reset_index(drop=True),after[before.columns].reset_index(drop=True))
            expected = before['ClinicalSignificance'].map(exact_class).astype('string').fillna('')
            if not expected.reset_index(drop=True).equals(after['exact_clinical_class'].reset_index(drop=True)):
                raise ValueError('Exact class does not match source label')
            flags = expected.ne('').map({True:'True',False:'False'}).astype('string')
            if not flags.reset_index(drop=True).equals(after['five_class_eligible'].reset_index(drop=True)):
                raise ValueError('Eligibility mismatch')
            checked += len(after)
        if checked != count: raise ValueError('Exact export row count mismatch')
        temporary.replace(destination)
    finally:
        if temporary.exists(): temporary.unlink()
    return {'file':destination.name,'records':count,'five_class_eligible':eligible,'unresolved':count-eligible}
