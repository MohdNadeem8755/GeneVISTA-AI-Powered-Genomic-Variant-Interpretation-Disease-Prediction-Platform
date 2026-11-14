"""Retrain biological-feature ensemble; reused internal evaluation, never external validation."""
import argparse
import sqlite3
from pathlib import Path
import joblib
import numpy as np
from catboost import CatBoostClassifier
from sklearn.base import BaseEstimator, ClassifierMixin
from sklearn.calibration import CalibratedClassifierCV
from sklearn.frozen import FrozenEstimator
from src.balanced_annotation import read_split, estimator, fit_weighted, evaluate
from src.robust_variant_training import metrics, save_json

ROOT = Path(__file__).resolve().parents[1]


class AnnotationVote(ClassifierMixin, BaseEstimator):
    def __init__(self, models, weights):
        self.models = models
        self.weights = weights
        self.classes_ = models[0].classes_

    def fit(self, X, y):
        return self

    def predict_proba(self, X):
        result = []
        for model in self.models:
            order = [list(model.classes_).index(c) for c in self.classes_]
            result.append(model.predict_proba(X)[:, order])
        return np.average(result, axis=0, weights=self.weights)

    def predict(self, X):
        return self.classes_[self.predict_proba(X).argmax(axis=1)]


def run(output, modulus=5, iterations=300):
    out = Path(output)
    out.mkdir(parents=True, exist_ok=False)
    incumbent = joblib.load(ROOT / 'models/annotation_prediction_v4/predictor.joblib')
    features, vocab = incumbent['features'], incumbent['gene_vocab']
    database = ROOT / 'data/processed/annotation_v3/full_cohort.sqlite'
    with sqlite3.connect(database.as_uri() + '?mode=ro', uri=True) as con:
        overlap = con.execute("SELECT COUNT(*) FROM (SELECT chrom,pos FROM variants WHERE split IN ('train','validation','calibration','test') GROUP BY chrom,pos HAVING COUNT(DISTINCT split)>1)").fetchone()[0]
        if overlap:
            raise ValueError('Locus overlap detected')
        print('Loading annotated training and validation records', flush=True)
        x, y = read_split(con, 'train', features, vocab, modulus)
        xv, yv = read_split(con, 'validation', features, vocab, 5)
        classes, counts = np.unique(y, return_counts=True)
        prior = counts / counts.sum()
        save_json(out / 'audit.json', {'train_rows': len(y), 'validation_rows': len(yv),
                                     'modulus': modulus, 'position_overlap': overlap,
                                     'features': list(features), 'evaluation_status': 'reused_internal_test',
                                     'objective': 'balanced_accuracy_then_macro_f1'})
        models, reports = [], []
        for power, depth in [(.5, 6), (1., 8)]:
            weights = np.array([(len(y) / (len(classes) * counts[list(classes).index(c)])) ** power for c in y])
            model = CatBoostClassifier(iterations=iterations, depth=depth, learning_rate=.05,
                                       l2_leaf_reg=10, thread_count=4, random_seed=42,
                                       verbose=False, allow_writing_files=False)
            print(f'Training CatBoost depth={depth}, weight power={power}', flush=True)
            model.fit(x, y, sample_weight=weights)
            models.append(model)
            report = metrics(yv, model.predict(xv).ravel(), classes)
            reports.append({'model': f'catboost_{depth}', **report})
            save_json(out / 'development.json', reports)
        print('Training histogram gradient boosting', flush=True)
        model = fit_weighted(estimator(31, iterations, features), x, y, 1.)
        models.append(model)
        reports.append({'model': 'histogram_boosting', **metrics(yv, model.predict(xv), classes)})
        del x
        # Tune mixture on validation only; calibration is a separate fixed partition.
        xc, yc = read_split(con, 'calibration', features, vocab, 5)
        calibrated = [CalibratedClassifierCV(FrozenEstimator(m), method='sigmoid').fit(xc, yc) for m in models]
        candidates = []
        best = None
        for weights in [(1,1,1), (2,1,1), (1,2,1), (1,1,2)]:
            vote = AnnotationVote(calibrated, weights)
            p = vote.predict_proba(xv)
            for exponent in (0., .25, .5, .75, 1.):
                result = evaluate(yv, p, vote.classes_, prior, exponent)
                candidate = {'weights': list(weights), 'exponent': exponent, **result}
                candidates.append(candidate)
                key = (result['balanced_accuracy'], result['macro_f1'])
                if best is None or key > best[0]:
                    best = (key, vote, candidate)
        old = evaluate(yv, incumbent['classifier'].predict_proba(xv), incumbent['classifier'].classes_,
                       incumbent['class_priors'], incumbent['decision_exponent'])
        save_json(out / 'development.json', {'individual': reports, 'ensembles': candidates, 'incumbent': old})
        _, vote, selected = best
        save_json(out / 'selection_frozen.json', selected)
        joblib.dump({'classifier': vote, 'features': features, 'gene_vocab': vocab,
                     'class_priors': prior, 'decision_exponent': selected['exponent']}, out / 'candidate.joblib')
        print('Selection frozen; evaluating reused internal test', flush=True)
        xt, yt = read_split(con, 'test', features, vocab, test_authorized=True)
        p = vote.predict_proba(xt)
        predicted = vote.classes_[(p / prior[None, :] ** selected['exponent']).argmax(axis=1)]
        new_report = metrics(yt, predicted, vote.classes_)
        previous = evaluate(yt, incumbent['classifier'].predict_proba(xt), incumbent['classifier'].classes_,
                            incumbent['class_priors'], incumbent['decision_exponent'])
        result = {'candidate': new_report, 'incumbent_same_test': previous,
                  'validation_selected': selected, 'test_status': 'previously_inspected_internal_test',
                  'target_met': bool(new_report['accuracy'] > .95 and new_report['balanced_accuracy'] > .95),
                  'production_model_changed': False}
        save_json(out / 'results.json', result)
        print({k: new_report[k] for k in ('accuracy','balanced_accuracy','macro_f1')}, flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True)
    parser.add_argument('--modulus', type=int, default=5)
    parser.add_argument('--iterations', type=int, default=300)
    args = parser.parse_args()
    if args.modulus < 1 or args.iterations < 1:
        parser.error('Positive modulus and iterations required')
    run(args.output, args.modulus, args.iterations)


if __name__ == '__main__':
    from src.retrain_annotation_ensemble import main as entrypoint
    entrypoint()
