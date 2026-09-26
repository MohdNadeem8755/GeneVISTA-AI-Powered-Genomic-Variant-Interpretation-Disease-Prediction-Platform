"""Train-only preprocessing and validation-only ANN experiments for annotation labels.

Does not read the already-used test set or change the deployed model.
"""
import argparse
from contextlib import closing
import copy
import json
from pathlib import Path
import sqlite3
import time

import joblib
import numpy as np
from sklearn.compose import ColumnTransformer
from sklearn.impute import SimpleImputer
from sklearn.neural_network import MLPClassifier
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler

from src.annotation_pipeline import meta
from src.balanced_annotation import evaluate, read_split
from src.variant_features import CLASSES

ROOT = Path(__file__).resolve().parents[1]
CONFIGS = [('ann64_sqrt', (64, 32), .5), ('ann128_balanced', (128, 64), 1.),
           ('ann128_sqrt', (128, 64), .5)]


def normalized_probabilities(p):
    """Correct float32 softmax roundoff before mixing with float64 outputs."""
    p = np.asarray(p, dtype=np.float64)
    if (p.ndim != 2 or not np.isfinite(p).all() or (p < 0).any()
            or (p > 1).any() or not np.allclose(p.sum(axis=1), 1, atol=1e-6, rtol=0)):
        raise ValueError('Invalid class probabilities')
    return p / p.sum(axis=1, keepdims=True)


def transform_inputs(x, features):
    """Log frequencies/counts without treating an absent annotation as zero."""
    x = np.asarray(x, dtype=np.float32).copy()
    if x.ndim != 2 or x.shape[1] != len(features) + 1:
        raise ValueError('Unexpected feature matrix shape')
    if np.isinf(x).any():
        raise ValueError('Infinite feature values')
    for i, name in enumerate(features):
        if name.startswith('gnomAD'):
            present = ~np.isnan(x[:, i])
            values = x[present, i]
            if (values < 0).any():
                raise ValueError('Negative population annotation')
            if name.endswith('_AF'):
                if (values > 1).any():
                    raise ValueError('Allele frequency above one')
                x[present, i] = np.log10(np.maximum(values, 1e-9))
            else:
                x[present, i] = np.log1p(values)
    return x


def preprocessor(n_features):
    numeric = make_pipeline(SimpleImputer(strategy='median', add_indicator=True,
                                          keep_empty_features=True), StandardScaler())
    gene = make_pipeline(SimpleImputer(strategy='constant', fill_value=-1,
                                       keep_empty_features=True),
                         OneHotEncoder(handle_unknown='ignore', sparse_output=False,
                                       dtype=np.float32))
    return ColumnTransformer([('numeric', numeric, list(range(n_features))),
                              ('gene', gene, [n_features])])


def class_weights(y, exponent):
    if not 0 <= exponent <= 1 or set(y) != set(CLASSES):
        raise ValueError('Five classes and a weight exponent in [0, 1] required')
    labels, counts = np.unique(y, return_counts=True)
    lookup = dict(zip(labels, counts))
    weights = np.asarray([(len(y) / (len(labels) * lookup[v])) ** exponent
                          for v in y], dtype=np.float32)
    return weights / weights.mean()


def objective(result):
    """Optimize the weaker of accuracy and balanced accuracy, then macro F1."""
    return (min(result['accuracy'], result['balanced_accuracy']), result['macro_f1'])


def development_split(con, split, features, vocab, modulus):
    if split not in ('train', 'validation'):
        raise ValueError('ANN development can only read train and validation')
    return read_split(con, split, features, vocab, modulus)


def run(database, incumbent_path, output, epochs=30, train_modulus=10, candidate=None):
    if epochs < 1 or train_modulus < 1:
        raise ValueError('Positive epoch limit and sampling modulus required')
    configs = [config for config in CONFIGS if candidate is None or config[0] == candidate]
    if not configs:
        raise ValueError('Unknown candidate')
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    with (output / 'started.json').open('x') as stream:
        json.dump({'test_rows_read': 0, 'objective': 'min(accuracy, balanced_accuracy)'}, stream)
    started = time.monotonic()
    incumbent = joblib.load(incumbent_path)
    features = tuple(incumbent['features'])
    vocab = incumbent['gene_vocab']
    with closing(sqlite3.connect(Path(database).resolve().as_uri() + '?mode=ro', uri=True)) as con:
        if meta(con, 'annotation_status') != 'ready':
            raise ValueError('Verified annotation import required')
        overlap = con.execute("SELECT COUNT(*) FROM (SELECT chrom,pos FROM variants WHERE split IN ('train','validation') GROUP BY chrom,pos HAVING COUNT(DISTINCT split)>1)").fetchone()[0]
        if overlap:
            raise ValueError('Train/validation locus overlap')
        print('Loading locus-sampled ANN train/validation; test remains unread', flush=True)
        x, y = development_split(con, 'train', features, vocab, train_modulus)
        xv, yv = development_split(con, 'validation', features, vocab, 5)
    classifier = incumbent['classifier']
    classes = np.asarray(classifier.classes_)
    prior = np.asarray(incumbent['class_priors'])
    incumbent_p = normalized_probabilities(classifier.predict_proba(xv))
    baseline = evaluate(yv, incumbent_p, classes, prior, incumbent['decision_exponent'])
    natural_baseline = evaluate(yv, incumbent_p, classes)
    prep = preprocessor(len(features))
    x = prep.fit_transform(transform_inputs(x, features)).astype(np.float32)
    xv = prep.transform(transform_inputs(xv, features)).astype(np.float32)
    report = {'status': 'training', 'train_rows': len(y), 'validation_rows': len(yv),
              'train_modulus': train_modulus, 'validation_modulus': 5,
              'test_rows_read': 0, 'production_model_changed': False,
              'position_overlap': overlap, 'incumbent_same_validation': baseline,
              'incumbent_natural_same_validation': natural_baseline,
              'objective': 'min(accuracy, balanced_accuracy), then macro F1',
              'scope': 'Internal development comparison; validation used for early stopping and selection, not an independent test. Upstream training overlap unresolved.',
              'features': list(features), 'candidates': []}

    def save():
        report['elapsed_seconds'] = time.monotonic() - started
        (output / 'report.json').write_text(json.dumps(report, indent=2, allow_nan=False), encoding='utf-8')

    print(json.dumps({'baseline_accuracy': baseline['accuracy'],
                      'baseline_balanced_accuracy': baseline['balanced_accuracy'],
                      'train': len(y), 'validation': len(yv)}), flush=True)
    save()
    best = None
    for name, layers, exponent in configs:
        model = MLPClassifier(hidden_layer_sizes=layers, alpha=.01, batch_size=1024,
                              learning_rate_init=.001, random_state=42, shuffle=True,
                              early_stopping=False)
        weights = class_weights(y, exponent)
        record = {'name': name, 'layers': list(layers), 'weight_exponent': exponent, 'epochs': []}
        report['candidates'].append(record)
        best_local = None
        stale = 0
        for epoch in range(1, epochs + 1):
            model.partial_fit(x, y, classes=classes, sample_weight=weights)
            p = normalized_probabilities(model.predict_proba(xv))
            np.testing.assert_array_equal(model.classes_, classes)
            result = {'epoch': epoch, **evaluate(yv, p, classes)}
            record['epochs'].append(result)
            if best_local is None or objective(result) > objective(best_local):
                best_local = result
                best_model = copy.deepcopy(model)
                best_p = p.copy()
                stale = 0
            else:
                stale += 1
            print(json.dumps({'name': name, 'epoch': epoch, 'accuracy': result['accuracy'],
                              'balanced_accuracy': result['balanced_accuracy']}), flush=True)
            record['selected'] = best_local
            save()
            if stale >= 7:
                break
        record['stopping_reason'] = 'validation_patience' if stale >= 7 else 'epoch_budget'
        joblib.dump({'classifier': best_model, 'preprocessor': prep, 'features': features,
                     'gene_vocab': vocab, 'selection': best_local, 'research_only': True},
                    output / (name + '.joblib'), compress=3)
        # Explicitly exploratory ensembles: these weighted ANN outputs are not calibrated.
        blends = []
        for ann_fraction in (.25, .5, .75, 1.):
            blend = normalized_probabilities(ann_fraction * best_p + (1 - ann_fraction) * incumbent_p)
            result = {'ann_fraction': ann_fraction, **evaluate(yv, blend, classes)}
            blends.append(result)
        record['ensembles'] = blends
        choice = max(blends, key=objective)
        candidate = {'name': name, 'epoch': best_local['epoch'], **choice}
        if best is None or objective(candidate) > objective(best):
            best = candidate
        save()
    report['status'] = 'validation_complete'
    report['selected'] = best
    report['both_above_90_on_validation'] = min(best['accuracy'], best['balanced_accuracy']) > .9
    report['beats_incumbent_objective'] = objective(best) > objective(baseline)
    # Verify serialization and preprocessing before declaring completion.
    restored = joblib.load(output / (best['name'] + '.joblib'))
    restored_p = normalized_probabilities(restored['classifier'].predict_proba(xv))
    fraction = best['ann_fraction']
    reproduced = evaluate(yv, fraction * restored_p + (1 - fraction) * incumbent_p, classes)
    for metric in ('accuracy', 'balanced_accuracy', 'macro_f1'):
        np.testing.assert_allclose(reproduced[metric], best[metric], rtol=0, atol=1e-12)
    report['artifact_roundtrip_verified'] = True
    save()
    print(json.dumps({'selected': best['name'], 'accuracy': best['accuracy'],
                      'balanced_accuracy': best['balanced_accuracy'], 'test_rows_read': 0}), flush=True)
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--database', type=Path, default=ROOT / 'data/processed/annotation_v3/full_cohort.sqlite')
    parser.add_argument('--incumbent', type=Path, default=ROOT / 'models/annotation_prediction_v4/predictor.joblib')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--epochs', type=int, default=30)
    parser.add_argument('--train-modulus', type=int, default=10)
    parser.add_argument('--candidate', choices=[config[0] for config in CONFIGS])
    args = parser.parse_args()
    run(args.database, args.incumbent, args.output, args.epochs, args.train_modulus, args.candidate)
