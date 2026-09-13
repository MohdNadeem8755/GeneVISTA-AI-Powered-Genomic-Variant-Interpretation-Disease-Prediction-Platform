"""Inference for future annotation ranking, using the exact training features."""
from functools import lru_cache
from pathlib import Path
import csv
import joblib
from src.disease_ranking import feature_dict, csr32

MODEL_DIR = Path(__file__).resolve().parents[2] / 'models/disease_ranking_v1'
@lru_cache(maxsize=1)
def bundle():
    return joblib.load(MODEL_DIR / 'disease_ranker.joblib')

def evaluation():
    with (MODEL_DIR / 'temporal_evaluation.csv').open(encoding='utf-8') as f:
        rows = list(csv.DictReader(f))
    for row in rows:
        for key in ('average_precision','prevalence','precision_at_1000'):
            row[key] = float(row[key]) if row[key] else None
        for key in ('positive_events','train_positive_events','rows','review_size','iterations'):
            row[key] = int(row[key])
    return {'metadata': bundle()['metadata'], 'items': rows}

def rank(details):
    if not details or details.get('clinical_state') != 'vus':
        return {'status': 'ineligible', 'items': [], 'reason': 'This model ranks starting variants of uncertain significance only.'}
    if not (MODEL_DIR / 'disease_ranker.joblib').exists():
        return {'status': 'unavailable', 'items': [], 'reason': 'Disease ranking model is not installed.'}
    model = bundle()
    row = dict(details, number_submitters=details.get('number_submitters_raw',0))
    features = feature_dict(row)
    vocabulary = model['vectorizer'].vocabulary_
    if not any(key.startswith('gene=') and key in vocabulary for key in features):
        return {'status':'unsupported_gene','items':[], 'reason':'No input gene was seen during model training.'}
    x = csr32(model['vectorizer'].transform([features]))
    meta = model['metadata']
    items = [{'disease_id': disease, 'disease_name': meta['names'][disease],
              'score': float(classifier.decision_function(x)[0])}
             for disease, classifier in zip(meta['classes'],model['classifiers'])]
    items.sort(key=lambda item: (-item['score'],item['disease_id']))
    return {'status':'available', 'items':items, 'score_type':meta['score_type'],
            'target':meta['target'], 'source_snapshot':details.get('snapshot'),
            'limitations':meta['limitations']}
