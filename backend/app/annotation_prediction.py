"""Optional local annotation inference with exact cached-record matching."""
from functools import lru_cache
import os
from pathlib import Path
import sqlite3
import json
import joblib
from src.train_annotation_candidates import matrix
from src.variant_features import CLASSES

ROOT = Path(__file__).resolve().parents[2]
MODEL = ROOT / 'models/annotation_prediction_v2/predictor.joblib'

@lru_cache(maxsize=1)
def bundle():
    return joblib.load(MODEL)

def predict(details, variant_id):
    database = Path(os.environ.get('GENEVISTA_ANNOTATION_DB', ROOT / 'data/processed/annotation_v2/download_verified_cohort.sqlite')).resolve()
    if not MODEL.is_file() or not database.is_file():
        return {'status': 'not_configured', 'reason': 'The optional annotation model and lookup are not installed on this server.'}
    if not details:
        return {'status': 'insufficient_evidence', 'reason': 'Current variant details are required.'}
    with sqlite3.connect(database.as_uri() + '?mode=ro', uri=True) as con:
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
    probabilities = model['classifier'].predict_proba(matrix([(None,scores,gene,name,kind)], model['features']))[0]
    values = dict(zip(model['classifier'].classes_, map(float, probabilities)))
    winner = max(values, key=values.get)
    return {'status': 'available', 'model_id': 'annotation_prediction_v2',
        'predicted_class': winner, 'top_probability': values[winner],
        'probabilities': [{'class_name': c, 'probability': values[c]} for c in CLASSES],
        'cohort_membership': split, 'allele': f'GRCh38 {chrom}:{pos} {ref}>{alt}',
        'disagrees_with_recorded': winner != details.get('clinical_significance'),
        'evaluation': {k:model['metadata']['test'][k] for k in ('accuracy','balanced_accuracy','macro_f1')},
        'test_count': model['metadata']['test_count'], 'scope': model['metadata']['scope'],
        'caution': model['metadata']['upstream_overlap']}
