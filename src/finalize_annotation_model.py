"""Calibrate the validation-selected candidate and evaluate the locked test once."""
import hashlib
import json
from pathlib import Path
import joblib
import numpy as np
from sklearn.calibration import CalibratedClassifierCV
from sklearn.frozen import FrozenEstimator
from src.annotation_pipeline import connect, meta
from src.train_annotation_candidates import matrix, metrics, ROOT
from src.variant_features import CLASSES


def finalize(database, experiments, output):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    lock = output / 'evaluation_started.json'
    selection = json.loads((experiments / 'full_validation.json').read_text())
    if selection['status'] != 'validation_complete':
        raise ValueError('Validation selection not complete')
    candidate_path = experiments / 'validation_candidate.joblib'
    candidate = joblib.load(candidate_path)
    if candidate['selection'] != selection['selected']:
        raise ValueError('Saved model does not match validation selection')
    with connect(database) as con:
        if meta(con, 'annotation_status') != 'ready':
            raise ValueError('Annotation import incomplete')
        if meta(con, 'annotation')['md5'] != candidate['annotation']['md5']:
            raise ValueError('Annotation provenance mismatch')
        overlap = con.execute("""SELECT COUNT(*) FROM
            (SELECT chrom,pos FROM variants WHERE split NOT IN ('excluded','reserved')
             GROUP BY chrom,pos HAVING COUNT(DISTINCT split)>1)""").fetchone()[0]
        if overlap:
            raise ValueError('Position overlap across partitions')
        # An exclusive marker prevents silent repeated test-based model selection.
        with lock.open('x', encoding='utf-8') as stream:
            json.dump({'candidate_sha256': hashlib.sha256(candidate_path.read_bytes()).hexdigest(),
                       'selection_metric': 'validation_macro_f1',
                       'calibration_method': 'sigmoid', 'test_used_for_selection': False}, stream)
        rows = con.execute("""SELECT v.label,a.scores,v.gene,v.name,v.variant_type
            FROM variants v JOIN annotations a USING(chrom,pos,ref,alt)
            WHERE v.split='calibration' ORDER BY v.chrom,v.pos,v.ref,v.alt""").fetchall()
        y = np.asarray([r[0] for r in rows])
        if set(y) != set(CLASSES):
            raise ValueError('Calibration requires all five classes')
        model = CalibratedClassifierCV(FrozenEstimator(candidate['classifier']), method='sigmoid')
        model.fit(matrix(rows, candidate['features']), y)
        calibration_count = len(rows)
        print(f'Calibrated on {calibration_count} rows', flush=True)
        rows = con.execute("""SELECT v.label,a.scores,v.gene,v.name,v.variant_type
            FROM variants v JOIN annotations a USING(chrom,pos,ref,alt)
            WHERE v.split='test' ORDER BY v.chrom,v.pos,v.ref,v.alt""").fetchall()
        y = np.asarray([r[0] for r in rows])
        if set(y) != set(CLASSES):
            raise ValueError('Test requires all five classes')
        probabilities = model.predict_proba(matrix(rows, candidate['features']))
        if not np.isfinite(probabilities).all() or not np.allclose(probabilities.sum(axis=1), 1):
            raise ValueError('Invalid model probabilities')
        test = metrics(y, probabilities, model.classes_)
        report = {'model_id': 'annotation_prediction_v2', 'status': 'research_evaluated',
            'validation_selection': selection['selected'], 'test': test,
            'test_count': len(rows), 'calibration_count': calibration_count,
            'test_total_eligible': con.execute("SELECT COUNT(*) FROM variants WHERE split='test'").fetchone()[0],
            'majority_class_accuracy': float(max(np.unique(y, return_counts=True)[1]) / len(y)),
            'features': list(candidate['features']), 'annotation': candidate['annotation'],
            'scope': selection['scope'], 'upstream_overlap': selection['upstream_overlap'],
            'calibration_method': 'sigmoid on separate calibration positions',
            'test_used_for_selection': False, 'split_position_overlap': overlap,
            'deployment_scope': 'Only exact GRCh38 allele matches with imported annotations. Other variants require the sequence model or abstention.',
            'production_model_changed': False}
        joblib.dump({'classifier': model, 'features': candidate['features'], 'metadata': report},
                    output / 'predictor.joblib', compress=3)
        (output / 'model_card.json').write_text(json.dumps(report, indent=2))
        print(json.dumps({k: test[k] for k in ('accuracy','balanced_accuracy','macro_f1','log_loss')}), flush=True)
    return report


if __name__ == '__main__':
    finalize(ROOT/'data/processed/annotation_v2/download_verified_cohort.sqlite',
             ROOT/'reports/annotation_v2/full_candidates', ROOT/'models/annotation_prediction_v2')
