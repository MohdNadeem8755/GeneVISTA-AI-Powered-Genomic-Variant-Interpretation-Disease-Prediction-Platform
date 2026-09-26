"""Validation-only CatBoost ablations using cached protein language-model scores."""
import argparse
from contextlib import closing
import json
from pathlib import Path
import sqlite3
import time
import joblib
import numpy as np
import pandas as pd
from src.annotation_pipeline import bucket, meta
from src.train_annotation_candidates import matrix
from src.protein_lm_features import FEATURES as PROTEIN_FEATURES
from src.ann_annotation import class_weights, objective, normalized_probabilities
from src.balanced_annotation import evaluate
from src.variant_features import CLASSES

ROOT=Path(__file__).resolve().parents[1]


def load_development(database,feature_database,split,features,modulus):
    if split not in ('train','validation'):raise ValueError('Only development splits may be read')
    arrays=[];labels=[];genes=[]
    with closing(sqlite3.connect(Path(database).resolve().as_uri()+'?mode=ro',uri=True)) as con:
        if meta(con,'annotation_status')!='ready':raise ValueError('Annotations are not ready')
        con.create_function('pilot',2,lambda chrom,pos:bucket(f'balanced-v4-pilot|{chrom}|{pos}',modulus))
        protein=bool(feature_database)
        if protein:
            con.execute('ATTACH DATABASE ? AS protein',(Path(feature_database).resolve().as_uri()+'?mode=ro',))
            metadata=con.execute("SELECT value FROM protein.metadata WHERE key='protein_features'").fetchone()
            if not metadata:raise ValueError('Protein feature provenance missing')
            metadata=json.loads(metadata[0])
            if metadata['status']!='ready' or metadata['md5']!=meta(con,'annotation')['md5']:raise ValueError('Protein feature provenance mismatch')
        extra=',p.scores' if protein else ''
        join=' LEFT JOIN protein.features p USING(chrom,pos,ref,alt)' if protein else ''
        cursor=con.execute('SELECT v.label,a.scores,v.gene,v.name,v.variant_type'+extra+' FROM variants v JOIN annotations a USING(chrom,pos,ref,alt)'+join+" WHERE v.split=? AND pilot(v.chrom,v.pos)=0 ORDER BY v.chrom,v.pos,v.ref,v.alt",(split,))
        while rows:=cursor.fetchmany(10000):
            base=matrix([r[:5] for r in rows],features)
            if protein:
                if any(r[5] is None for r in rows):raise ValueError('Missing matched protein feature row')
                parsed=[json.loads(r[5]) for r in rows]
                extra_values=np.asarray([[scores[name] for name in PROTEIN_FEATURES] for scores in parsed],dtype=np.float32)
                base=np.column_stack((base,extra_values))
            arrays.append(base);labels.extend(r[0] for r in rows);genes.extend(r[2] or '__missing__' for r in rows)
    if set(labels)!=set(CLASSES):raise ValueError('All five classes required')
    names=list(features)+(list(PROTEIN_FEATURES) if protein else [])
    frame=pd.DataFrame(np.concatenate(arrays),columns=names)
    frame['gene']=genes
    return frame,np.asarray(labels)


def run(database,feature_database,output,train_modulus=10,scope='both',weights=(.5,1.),iterations=500):
    from catboost import CatBoostClassifier, Pool
    if train_modulus<1 or iterations<1 or not weights or any(not 0<=w<=1 for w in weights):raise ValueError('Invalid experiment settings')
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    with (output/'started.json').open('x') as f:json.dump({'test_rows_read':0,'objective':'min(accuracy, balanced_accuracy)'},f)
    started=time.monotonic();incumbent=joblib.load(ROOT/'models/annotation_prediction_v4/predictor.joblib')
    features=tuple(incumbent['features'])
    with closing(sqlite3.connect(Path(database).resolve().as_uri()+'?mode=ro',uri=True)) as con:
        overlap=con.execute("SELECT COUNT(*) FROM (SELECT chrom,pos FROM variants WHERE split IN ('train','validation') GROUP BY chrom,pos HAVING COUNT(DISTINCT split)>1)").fetchone()[0]
        if overlap:raise ValueError('Locus leakage')
    scopes=['baseline','protein_lm'] if scope=='both' else [scope]
    report={'status':'training','test_rows_read':0,'calibration_rows_read':0,'production_model_changed':False,
            'train_modulus':train_modulus,'validation_modulus':5,'position_overlap':overlap,'candidates':[],
            'scope':'Internal validation reused for development; not independent test accuracy or clinical validation.',
            'llm_use':'Precomputed ESM1b protein language-model and MutFormer BERT scores; no chat LLM, synthetic labels, API calls, or local transformer training.',
            'upstream_overlap':'Protein sequence models and supervised source scores may overlap benchmark data. No independent clinical accuracy claim.'}
    def save():
        report['seconds']=time.monotonic()-started
        (output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False),encoding='utf-8')
    best=None
    for current_scope in scopes:
        extra=feature_database if current_scope=='protein_lm' else None
        print(f'Loading {current_scope} development data',flush=True)
        x,y=load_development(database,extra,'train',features,train_modulus)
        xv,yv=load_development(database,extra,'validation',features,5)
        report.update(train_rows=len(y),validation_rows=len(yv))
        gene_codes=xv['gene'].map(incumbent['gene_vocab']).to_numpy(dtype=np.float32)
        oldx=np.column_stack((xv[list(features)].to_numpy(),gene_codes))
        report['incumbent_same_validation']=evaluate(yv,incumbent['classifier'].predict_proba(oldx),incumbent['classifier'].classes_,incumbent['class_priors'],incumbent['decision_exponent'])
        del oldx
        if current_scope=='protein_lm':
            report['protein_coverage']={name:int(x[name].notna().sum()) for name in PROTEIN_FEATURES}
        validation=Pool(xv,label=yv,cat_features=['gene'])
        save()
        for exponent in weights:
            name=f'{current_scope}_weight_{exponent:g}'
            print(f'Training {name} on {len(y):,} rows',flush=True)
            train=Pool(x,label=y,cat_features=['gene'],weight=class_weights(y,exponent))
            model=CatBoostClassifier(iterations=iterations,depth=7,learning_rate=.05,l2_leaf_reg=10,
                loss_function='MultiClass',random_seed=42,thread_count=4,allow_writing_files=False,
                verbose=100,nan_mode='Min',one_hot_max_size=10)
            model.fit(train)
            np.testing.assert_array_equal(model.classes_,np.sort(CLASSES))
            checkpoints=[]
            for trees in sorted(set([min(100,iterations),min(250,iterations),iterations])):
                p=normalized_probabilities(model.predict_proba(validation,ntree_end=trees))
                checkpoints.append({'trees':trees,**evaluate(yv,p,model.classes_)})
            selected=max(checkpoints,key=objective)
            row={'name':name,'feature_scope':current_scope,'weight_exponent':exponent,'checkpoints':checkpoints,**selected}
            report['candidates'].append(row)
            if best is None or objective(row)>objective(best):
                best=row;model.shrink(selected['trees'])
                model.save_model(str(output/'candidate.cbm'))
                restored=CatBoostClassifier();restored.load_model(str(output/'candidate.cbm'))
                np.testing.assert_allclose(restored.predict_proba(validation),model.predict_proba(validation),atol=1e-12,rtol=1e-12)
                report['artifact_roundtrip_verified']=True
                joblib.dump({'model_file':'candidate.cbm','features':list(x.columns),'base_features':features,
                             'feature_scope':current_scope,'classes':list(model.classes_),
                             'annotation_md5':incumbent['metadata']['annotation']['md5'],
                             'probability_status':'uncalibrated weighted research classifier',
                             'selection':selected,'research_only':True},output/'candidate_metadata.joblib')
            report['selected']=best;save()
            print(json.dumps({k:row[k] for k in ('name','accuracy','balanced_accuracy','macro_f1','trees')}),flush=True)
            del model,train
        del x,y,xv,yv,validation
    report['status']='validation_complete'
    report['both_above_90']=min(best['accuracy'],best['balanced_accuracy'])>.9
    save();return report


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--database',type=Path,default=ROOT/'data/processed/annotation_v3/full_cohort.sqlite')
    p.add_argument('--protein-features',type=Path,default=ROOT/'data/processed/protein_lm_v1/features.sqlite')
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--train-modulus',type=int,default=10)
    p.add_argument('--scope',choices=['both','baseline','protein_lm'],default='both')
    p.add_argument('--weights',type=float,nargs='+',default=[.5,1.])
    p.add_argument('--iterations',type=int,default=500)
    a=p.parse_args();run(a.database,a.protein_features,a.output,a.train_modulus,a.scope,a.weights,a.iterations)
