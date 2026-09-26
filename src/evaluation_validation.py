"""Input and artifact checks for frozen temporal-model evaluation."""
from pathlib import Path
from types import SimpleNamespace
import hashlib
import json
import numpy as np
import pandas as pd
from src.feature_validation import load_label_manifest, load_labels
from src.training_validation import validate_config, validate_manifest, KEYS


def safe_file(folder, name):
    if not isinstance(name,str) or not name or Path(name).name != name or '/' in name or '\\' in name:
        raise ValueError('Expected a local filename')
    return Path(folder)/name


def resolve_bundle(root):
    root=Path(root).resolve()
    report=json.loads((root/'reports/training04_validation.json').read_text(encoding='utf-8'))
    if report.get('status')!='passed': raise ValueError('No successfully validated training run')
    bundle=(root/report['evaluation_bundle']).resolve()
    if not bundle.is_relative_to(root/'models/temporal_v1'):
        raise ValueError('Evaluation bundle must be inside models/temporal_v1')
    if not (bundle/'evaluation_plan.json').is_file(): raise FileNotFoundError(bundle)
    return bundle


def sha256(path):
    digest=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):digest.update(block)
    return digest.hexdigest()


def verify_development(root,bundle):
    root=Path(root).resolve()
    entries=json.loads((Path(bundle)/'development_inputs.json').read_text())
    if not entries: raise ValueError('Missing development provenance')
    for entry in entries:
        path=(root/entry['file']).resolve()
        if not path.is_relative_to(root/'data/processed') or sha256(path)!=entry['sha256']:
            raise ValueError('Development inputs changed since model freezing')
    return entries


def load_evaluation_data(processed,plan):
    processed=Path(processed);config=plan['feature_config'];validate_config(config,temporal=True)
    labels_manifest=load_label_manifest(processed/'model_labels',processed)
    baseline_manifest=pd.read_csv(processed/'model_features/feature_manifest.csv',dtype='string',keep_default_na=False)
    temporal_manifest=pd.read_csv(processed/'temporal_features_v1/temporal_manifest.csv',dtype='string',keep_default_na=False)
    validate_manifest(baseline_manifest,labels_manifest)
    validate_manifest(temporal_manifest,labels_manifest,baseline_manifest)
    test=temporal_manifest.loc[temporal_manifest['split'].eq('test')]
    if len(test)!=1: raise ValueError('Evaluation requires exactly one frozen test interval')
    row=test.iloc[0];meta=labels_manifest.loc[labels_manifest['feature_snapshot'].eq(row.feature_snapshot)].iloc[0]
    if plan['test_interval']!=f'{meta.feature_snapshot}_to_{meta.outcome_snapshot}' or int(plan['expected_test_rows'])!=int(row.rows):
        raise ValueError('Test interval or row count differs from frozen plan')
    expected={'baseline_features':row.feature_file,'temporal_features':row.temporal_file,'labels':row.label_file}
    if plan['test_files']!=expected:raise ValueError('Test files differ from manifests')
    if not plan['review_sizes'] or any(type(v)!=int or v<=0 for v in plan['review_sizes']):raise ValueError('Review sizes must be positive integers')
    combined=load_labels(processed/'model_labels',SimpleNamespace(**meta.to_dict()))
    for folder,name in [('model_features',row.feature_file),('temporal_features_v1',row.temporal_file)]:
        table=pd.read_csv(safe_file(processed/folder,name),dtype='string',keep_default_na=False)
        if len(table)!=int(row.rows) or table.duplicated(KEYS).any() or not table['feature_snapshot'].eq(row.feature_snapshot).all():raise ValueError('Invalid feature keys/count/snapshot')
        if set(table.columns)&(set(combined.columns)-set(KEYS)):raise ValueError('Unexpected overlapping feature columns')
        combined=combined.merge(table,on=KEYS,how='outer',validate='one_to_one',indicator=True)
        if not combined.pop('_merge').eq('both').all():raise ValueError('Unmatched feature/label IDs')
    X=combined[config['categorical_features']+config['numeric_features']].copy()
    for c in config['numeric_features']:
        X[c]=pd.to_numeric(X[c].replace('',np.nan),errors='raise').astype('float64')
        if np.isinf(X[c]).any():raise ValueError('Infinite numeric feature')
    for c in config['categorical_features']:
        values=X[c].astype(object);X[c]=values.where(values.notna()&values.ne(''),np.nan)
    totals={c:int(pd.to_numeric(labels_manifest.loc[labels_manifest['split'].eq('train'),c]).sum()) for c in config['target_classes']}
    majority=max(sorted(totals),key=totals.get)
    return X,combined['target'].copy(),combined[KEYS].copy(),majority


def verify_predictions(path,expected,chunksize=100000):
    offset=0
    for chunk in pd.read_csv(path,dtype='string',keep_default_na=False,chunksize=chunksize):
        part=expected.iloc[offset:offset+len(chunk)].reset_index(drop=True)
        if len(part)!=len(chunk) or list(part.columns)!=list(chunk.columns):raise ValueError('Prediction artifact schema/length mismatch')
        for c in chunk:
            if c.startswith('score_'):
                np.testing.assert_allclose(pd.to_numeric(chunk[c],errors='raise').to_numpy(dtype=float),part[c].to_numpy(dtype=float),rtol=1e-12,atol=1e-12)
            else:
                pd.testing.assert_series_equal(chunk[c].reset_index(drop=True),part[c].astype('string'),check_names=False)
        offset+=len(chunk)
    if offset!=len(expected):raise ValueError('Truncated prediction artifact')


def fixed_balanced_accuracy(actual, predictions):
    from sklearn.metrics import recall_score
    classes = ['stayed_vus', 'became_benign', 'became_pathogenic']
    if set(actual) != set(classes):
        return float('nan')
    return float(recall_score(actual, predictions, labels=classes, average='macro', zero_division=0))
