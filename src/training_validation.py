"""Checked data loading and evaluation for temporal VUS training."""
from pathlib import Path
from types import SimpleNamespace
import warnings
import numpy as np
import pandas as pd
from sklearn.exceptions import ConvergenceWarning
from sklearn.metrics import average_precision_score
from src.feature_validation import load_labels

KEYS = ['VariationID', 'feature_snapshot']
BASE_CATEGORICAL = ['variant_type', 'review_status']
BASE_NUMERIC = ['number_submitters', 'submitter_count_missing', 'days_since_evaluation', 'evaluation_date_missing', 'evaluation_date_after_reference']
HISTORY_NUMERIC = ['previous_number_submitters', 'previous_submitter_count_missing', 'submitter_count_change', 'submitter_change_missing']


def validate_config(config, temporal=False):
    expected_cat = BASE_CATEGORICAL + (['previous_state'] if temporal else [])
    expected_num = BASE_NUMERIC + (HISTORY_NUMERIC if temporal else [])
    if config['categorical_features'] != expected_cat or config['numeric_features'] != expected_num:
        raise ValueError('Unexpected model features: identifiers, labels and future data are not allowed')
    if config['join_keys'] != KEYS or set(config['target_classes']) != {'stayed_vus', 'became_benign', 'became_pathogenic'} or len(config['target_classes']) != 3:
        raise ValueError('Invalid keys or target classes')


def validate_manifest(manifest, labels, baseline=None):
    if manifest.empty or not manifest['feature_snapshot'].is_unique:
        raise ValueError('Empty or duplicated feature manifest')
    expected = labels.rename(columns={'file':'label_file'})
    columns = ['label_file','feature_snapshot','split','rows']
    pd.testing.assert_frame_equal(manifest[columns].sort_values('feature_snapshot').reset_index(drop=True), expected[columns].sort_values('feature_snapshot').reset_index(drop=True))
    if baseline is not None:
        columns = ['feature_file','label_file','feature_snapshot','split','rows']
        pd.testing.assert_frame_equal(manifest[columns].sort_values('feature_snapshot').reset_index(drop=True), baseline[columns].sort_values('feature_snapshot').reset_index(drop=True))
    for column in ('feature_file','label_file','temporal_file'):
        if column in manifest:
            for name in manifest[column]:
                if not name or Path(name).name != name or '/' in name or '\\' in name:
                    raise ValueError('Manifest filename must stay inside its data folder')


def load_split(processed, manifest, label_manifest, config, split, temporal=False):
    if split not in ('train','validation'):
        raise ValueError('Training notebook cannot load the held-out test split')
    validate_config(config, temporal)
    validate_manifest(manifest,label_manifest)
    processed=Path(processed)
    parts=[]
    selected=manifest.loc[manifest['split'].eq(split)]
    if selected.empty: raise ValueError('Requested development split is empty')
    for row in selected.itertuples(index=False):
        metadata=label_manifest.loc[label_manifest['feature_snapshot'].eq(row.feature_snapshot)].iloc[0]
        labels=load_labels(processed/'model_labels',SimpleNamespace(**metadata.to_dict()))
        tables=[pd.read_csv(processed/'model_features'/row.feature_file,dtype='string',keep_default_na=False)]
        if temporal:
            tables.append(pd.read_csv(processed/'temporal_features_v1'/row.temporal_file,dtype='string',keep_default_na=False))
        combined=labels
        for table in tables:
            if len(table)!=int(row.rows) or table.duplicated(KEYS).any() or not table['feature_snapshot'].eq(row.feature_snapshot).all():
                raise ValueError('Feature row count, keys, or snapshot mismatch')
            if set(table.columns)&(set(combined.columns)-set(KEYS)):
                raise ValueError('Overlapping feature/label columns could leak targets')
            combined=combined.merge(table,on=KEYS,how='outer',validate='one_to_one',indicator=True)
            if not combined.pop('_merge').eq('both').all(): raise ValueError('Feature/label IDs do not match')
        parts.append(combined)
    data=pd.concat(parts,ignore_index=True)
    if data.duplicated(KEYS).any(): raise ValueError('Duplicate examples across source files')
    columns=config['categorical_features']+config['numeric_features']
    X=data[columns].copy()
    for column in config['numeric_features']:
        X[column]=pd.to_numeric(X[column].replace('',np.nan),errors='raise').astype('float64')
        if np.isinf(X[column]).any(): raise ValueError(f'Infinite numeric feature: {column}')
    for column in config['categorical_features']:
        values = X[column].astype(object)
        X[column] = values.where(values.notna() & values.ne(''), np.nan)
    y=data['target'].copy()
    if split=='train' and set(y)!=set(config['target_classes']):
        raise ValueError('Training split must contain every target class')
    return X,y,data[KEYS].copy()


def fit_converged(model,X,y):
    with warnings.catch_warnings():
        warnings.simplefilter('error',ConvergenceWarning)
        try: model.fit(X,y)
        except ConvergenceWarning as exc:
            raise RuntimeError('Model did not converge; do not evaluate or save it. Adjust solver settings and rerun.') from exc
    classifier=model.named_steps['classifier']
    if hasattr(classifier,'n_iter_') and isinstance(classifier.n_iter_,np.ndarray) and classifier.n_iter_.max()>=classifier.max_iter:
        raise RuntimeError('Solver reached its iteration limit; model is not approved for saving')
    return model


def safe_ratio(numerator,denominator):
    return float(numerator/denominator) if denominator else float('nan')


def ranking_metrics(actual,scores):
    actual=np.asarray(actual,dtype=bool);scores=np.asarray(scores,dtype=float)
    if actual.ndim!=1 or scores.shape!=actual.shape or not len(actual) or not np.isfinite(scores).all():
        raise ValueError('Invalid ranking arrays')
    prevalence=float(actual.mean())
    ap=float(average_precision_score(actual,scores)) if actual.any() else float('nan')
    return prevalence,ap,safe_ratio(ap,prevalence)


def select_primary(ranking,models):
    targets={'became_benign','became_pathogenic'}
    rows=[]
    for name in models:
        metrics=ranking.loc[ranking['model'].eq(name)]
        if len(metrics)!=2 or set(metrics['target'])!=targets or not np.isfinite(metrics['average_precision']).all():
            raise ValueError(f'Incomplete validation ranking for {name}')
        rows.append({'model':name,'mean_average_precision':float(metrics['average_precision'].mean())})
    selection=pd.DataFrame(rows).sort_values(['mean_average_precision','model'],ascending=[False,True]).reset_index(drop=True)
    return selection.loc[0,'model'],selection
