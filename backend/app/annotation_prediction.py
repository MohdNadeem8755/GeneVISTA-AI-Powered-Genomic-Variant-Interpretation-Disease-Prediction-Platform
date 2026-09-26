"""Optional local annotation inference with exact cached-record matching."""
from functools import lru_cache
from contextlib import closing
import os
from pathlib import Path
import sqlite3
import json
from datetime import datetime, timezone
import joblib
from src.train_annotation_candidates import matrix
from src.variant_features import CLASSES
from src.balanced_annotation import enhanced_matrix, decision_indices
from src.annotation_pipeline import bucket

ROOT = Path(__file__).resolve().parents[2]
MODEL = ROOT / 'models/annotation_prediction_v2/predictor.joblib'

def model_path():
    candidate = ROOT / 'models/annotation_prediction_v4/predictor.joblib'
    return candidate if candidate.is_file() else MODEL

@lru_cache(maxsize=2)
def _load_bundle(path, modified):
    return joblib.load(path)

def bundle():
    path = model_path()
    return _load_bundle(str(path), path.stat().st_mtime_ns)

def model_card():
    path = model_path()
    if not path.is_file():
        return {'status': 'not_configured'}
    return {**bundle()['metadata'],
            'model_revision': str(path.stat().st_mtime_ns),
            'model_saved_at': datetime.fromtimestamp(path.stat().st_mtime, timezone.utc).isoformat()}

def predict(details, variant_id):
    default_db = ROOT / ('data/processed/annotation_v3/full_cohort.sqlite' if model_path() != MODEL else 'data/processed/annotation_v2/download_verified_cohort.sqlite')
    database = Path(os.environ.get('GENEVISTA_ANNOTATION_DB', default_db)).resolve()
    if not model_path().is_file() or not database.is_file():
        return {'status': 'not_configured', 'reason': 'The optional annotation model and lookup are not installed on this server.'}
    if not details:
        return {'status': 'insufficient_evidence', 'reason': 'Current variant details are required.'}
    with closing(sqlite3.connect(database.as_uri() + '?mode=ro', uri=True)) as con:
        state = con.execute("SELECT value FROM metadata WHERE key='annotation_status'").fetchone()
        if not state or json.loads(state[0]) != 'ready':
            return {'status': 'not_ready', 'reason': 'Annotation import is incomplete.'}
        provenance = json.loads(con.execute("SELECT value FROM metadata WHERE key='annotation'").fetchone()[0])
        rows = con.execute('''SELECT a.scores,v.gene,v.name,v.variant_type,v.split,v.chrom,v.pos,v.ref,v.alt
            FROM variants v JOIN annotations a USING(chrom,pos,ref,alt)
            WHERE v.variation_id=? AND v.split<>'excluded' ''', (str(variant_id),)).fetchall()
    if len(rows) != 1:
        return {'status': 'insufficient_evidence', 'reason': 'No unique annotated GRCh38 allele is available in this sampled lookup. Refer to the sequence model result.'}
    scores, gene, name, kind, split, chrom, pos, ref, alt = rows[0]
    if name != details.get('name') or gene != details.get('gene_symbol'):
        return {'status': 'insufficient_evidence', 'reason': 'Annotation does not match current variant details.'}
    model = bundle()
    if provenance['md5'] != model['metadata']['annotation']['md5']:
        return {'status': 'not_ready', 'reason': 'Annotation version does not match the evaluated model.'}
    rows = [(None,scores,gene,name,kind)]
    inputs = enhanced_matrix(rows, model['features'], model['gene_vocab']) if 'gene_vocab' in model else matrix(rows, model['features'])
    probabilities = model['classifier'].predict_proba(inputs)[0]
    values = dict(zip(model['classifier'].classes_, map(float, probabilities)))
    highest = max(values, key=values.get)
    winner = highest
    if 'class_priors' in model:
        index = decision_indices([probabilities], model['class_priors'], model['decision_exponent'])[0]
        winner = str(model['classifier'].classes_[index])
    if model['metadata'].get('model_id') == 'annotation_prediction_v4' and split == 'test' and bucket(f'genevista-annotation-v2-sample|{chrom}|{pos}',10)==0:
        split = 'legacy_test_not_in_current_evaluation'
    return {'status': 'available', 'model_id': model['metadata']['model_id'],
        'most_probable_class': highest,
        'decision_policy': 'The decision rule balances recall across all five classes; its selected class can differ from the largest probability.' if 'class_priors' in model else None,
        'predicted_class': winner, 'top_probability': values[winner],
        'probabilities': [{'class_name': c, 'probability': values[c]} for c in CLASSES],
        'cohort_membership': split, 'allele': f'GRCh38 {chrom}:{pos} {ref}>{alt}',
        'disagrees_with_recorded': winner != details.get('clinical_significance'),
        'evaluation': {k:model['metadata']['test'][k] for k in ('accuracy','balanced_accuracy','macro_f1')},
        'test_count': model['metadata']['test_count'], 'scope': model['metadata']['scope'],
        'caution': model['metadata']['upstream_overlap']}
