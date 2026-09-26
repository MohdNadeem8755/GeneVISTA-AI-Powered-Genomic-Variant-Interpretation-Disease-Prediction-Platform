"""Predict with a saved research candidate using exact cached allele matches."""
import argparse
from contextlib import closing
import json
from pathlib import Path
import sqlite3
import joblib
import numpy as np
import pandas as pd
from src.ann_annotation import normalized_probabilities
from src.annotation_pipeline import meta
from src.protein_lm_features import FEATURES as PROTEIN_FEATURES
from src.train_annotation_candidates import matrix
from src.variant_features import CLASSES

ROOT=Path(__file__).resolve().parents[1]


def predict(variation_id,directory,database,protein_database):
    from catboost import CatBoostClassifier
    directory=Path(directory)
    config=joblib.load(directory/'candidate_metadata.joblib')
    if config['feature_scope']!='protein_lm':raise ValueError('Protein candidate required')
    with closing(sqlite3.connect(Path(database).resolve().as_uri()+'?mode=ro',uri=True)) as con:
        con.execute('ATTACH DATABASE ? AS protein',(Path(protein_database).resolve().as_uri()+'?mode=ro',))
        manifest=json.loads(con.execute("SELECT value FROM protein.metadata WHERE key='protein_features'").fetchone()[0])
        if (meta(con,'annotation_status')!='ready' or manifest['status']!='ready'
            or manifest['md5']!=meta(con,'annotation')['md5'] or manifest['md5']!=config['annotation_md5']):
            raise ValueError('Candidate annotation provenance mismatch')
        rows=con.execute('''SELECT a.scores,v.gene,v.name,v.variant_type,p.scores,v.split,v.chrom,v.pos,v.ref,v.alt
            FROM variants v JOIN annotations a USING(chrom,pos,ref,alt)
            JOIN protein.features p USING(chrom,pos,ref,alt)
            WHERE v.variation_id=? AND v.split IN ('train','validation')''',(str(variation_id),)).fetchall()
    if len(rows)!=1:
        return {'status':'outside_research_cache','reason':'No unique development-cohort allele available; no prediction inferred from a different allele.'}
    scores,gene,name,kind,protein,split,chrom,pos,ref,alt=rows[0]
    values=matrix([(None,scores,gene,name,kind)],config['base_features'])[0].tolist()
    protein=json.loads(protein)
    values.extend(protein[name] for name in PROTEIN_FEATURES)
    expected=list(config['base_features'])+list(PROTEIN_FEATURES)+['gene']
    if config['features']!=expected:raise ValueError('Candidate feature order mismatch')
    frame=pd.DataFrame([values],columns=expected[:-1],dtype=np.float32)
    frame['gene']=gene or '__missing__'
    model=CatBoostClassifier();model.load_model(str(directory/'candidate.cbm'))
    if set(model.classes_)!=set(CLASSES) or list(model.classes_)!=config['classes']:raise ValueError('Candidate class mismatch')
    p=normalized_probabilities(model.predict_proba(frame))[0]
    return {'status':'research_prediction','model':'protein_lm_catboost_candidate',
            'variation_id':str(variation_id),'allele':f'GRCh38 {chrom}:{pos} {ref}>{alt}',
            'predicted_class':str(model.classes_[p.argmax()]),
            'class_estimates':dict(zip(map(str,model.classes_),map(float,p))),
            'cohort_membership':split,'evaluation_split':'development validation',
            'caution':'Uncalibrated weighted model outputs, not patient disease probabilities. Training-record predictions are not independent validation.'}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--variant',required=True)
    p.add_argument('--model-directory',type=Path,default=ROOT/'reports/catboost_protein_v6_full')
    p.add_argument('--database',type=Path,default=ROOT/'data/processed/annotation_v3/full_cohort.sqlite')
    p.add_argument('--protein-features',type=Path,default=ROOT/'data/processed/protein_lm_v1/features.sqlite')
    a=p.parse_args();print(json.dumps(predict(a.variant,a.model_directory,a.database,a.protein_features),indent=2))
