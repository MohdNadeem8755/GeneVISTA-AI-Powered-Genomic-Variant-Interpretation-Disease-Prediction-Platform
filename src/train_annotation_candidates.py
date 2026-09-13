"""Full matched-cohort validation experiments; never read locked test labels."""
import json
from pathlib import Path
import argparse
import joblib
import numpy as np
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import accuracy_score, balanced_accuracy_score, classification_report, log_loss
from src.annotation_pipeline import connect, meta, FEATURES
from src.variant_features import CLASSES, variant_features

ROOT = Path(__file__).resolve().parents[1]
# Restrict sequence inputs to shared numeric features, never names or labels.
SEQUENCE = ('transition', 'intronic', 'utr', 'log_cdna_position', 'synonymous',
            'stop_gain', 'log_protein_position', 'has_nucleotide_change',
            'has_protein_change', 'event=del', 'event=dup', 'event=ins', 'event=fs', 'event=inv')
INDEPENDENT = tuple(n for n in FEATURES if n.startswith('gnomAD') or n.startswith(('GERP', 'phyloP')))

def matrix(rows, features):
    result = []
    for label, scores, gene, name, kind in rows:
        values = json.loads(scores)
        sequence = variant_features({'gene_symbol': gene, 'name': name, 'variant_type': kind})
        result.append([values.get(n) if n in FEATURES else sequence.get(n) for n in features])
    return np.asarray(result, dtype=np.float32)

def metrics(y, probabilities, classes):
    predicted = np.asarray(classes)[np.argmax(probabilities, axis=1)]
    report = classification_report(y, predicted, labels=CLASSES, output_dict=True, zero_division=0)
    return {'accuracy': float(accuracy_score(y, predicted)),
            'balanced_accuracy': float(balanced_accuracy_score(y, predicted)),
            'macro_f1': report['macro avg']['f1-score'],
            'log_loss': float(log_loss(y, probabilities, labels=classes)), 'per_class': report}

def train(database, output):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    report_path = output / 'full_validation.json'
    if report_path.exists():
        raise ValueError('Output already evaluated; use a new experiment directory')
    with connect(database) as con:
        if meta(con, 'annotation_status') != 'ready':
            raise ValueError('Verified annotation import required')
        provenance = meta(con, 'annotation')
        splits = {}
        coverage = {}
        for split in ('train', 'validation'):
            rows = con.execute('''SELECT v.label,a.scores,v.gene,v.name,v.variant_type
                FROM variants v JOIN annotations a USING(chrom,pos,ref,alt)
                WHERE v.split=? ORDER BY v.chrom,v.pos,v.ref,v.alt''', (split,)).fetchall()
            if set(r[0] for r in rows) != set(CLASSES):
                raise ValueError('All five classes required')
            splits[split] = rows
            coverage[split] = {'matched': len(rows), 'total': con.execute(
                'SELECT COUNT(*) FROM variants WHERE split=?', (split,)).fetchone()[0]}
    experiments = []
    best = None
    for scope, features in [('population_conservation', INDEPENDENT + SEQUENCE),
                             ('all_annotations', FEATURES + SEQUENCE)]:
        x = {s: matrix(rows, features) for s, rows in splits.items()}
        y = {s: np.asarray([r[0] for r in rows]) for s, rows in splits.items()}
        for weight in (0.0, 0.5):
            for leaves in (15, 31):
                counts = dict(zip(*np.unique(y['train'], return_counts=True)))
                weights = np.asarray([(len(y['train']) / (5 * counts[label])) ** weight
                                      for label in y['train']])
                model = HistGradientBoostingClassifier(max_iter=250, max_leaf_nodes=leaves,
                    min_samples_leaf=30, l2_regularization=10, learning_rate=.07,
                    early_stopping=False, random_state=42)
                model.fit(x['train'], y['train'], sample_weight=weights)
                result = {'scope': scope, 'weight_exponent': weight, 'leaves': leaves,
                          'features': list(features), **metrics(y['validation'],
                              model.predict_proba(x['validation']), model.classes_)}
                experiments.append(result)
                print(json.dumps({k: result[k] for k in ('scope','weight_exponent','leaves','accuracy','macro_f1')}), flush=True)
                if best is None or result['macro_f1'] > best['macro_f1']:
                    best = result
                    joblib.dump({'classifier': model, 'features': features,
                                 'selection': result, 'annotation': provenance},
                                output / 'validation_candidate.joblib', compress=3)
                report_path.write_text(json.dumps({'status': 'validation_in_progress',
                    'candidates': experiments, 'coverage': coverage}, indent=2))
    report = {'status': 'validation_complete', 'candidates': experiments, 'selected': best,
        'coverage': coverage, 'test_rows_read': 0, 'calibration_rows_read': 0,
        'production_model_changed': False,
        'scope': 'Full matched internal same-snapshot SNVs; not an external or clinical validation.',
        'provenance': provenance,
        'upstream_overlap': 'All-annotation candidates include supervised pathogenicity scores with possible ClinVar training overlap. Population/conservation ablation excludes those scores. Neither establishes patient disease risk.'}
    report_path.write_text(json.dumps(report, indent=2))
    return report

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--database', type=Path, default=ROOT/'data/processed/annotation_v2/download_verified_cohort.sqlite')
    parser.add_argument('--output', type=Path, default=ROOT/'reports/annotation_v2/full_candidates')
    args = parser.parse_args()
    train(args.database, args.output)
