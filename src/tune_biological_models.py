"""Development-only CatBoost search and feature-specific soft-voting experiments."""
import argparse
from pathlib import Path
import joblib
import numpy as np
from catboost import CatBoostClassifier
from sklearn.model_selection import ParameterSampler
from src.biological_features import load_development
from src.robust_variant_training import metrics, save_json

ROOT=Path(__file__).resolve().parents[1]


def predict_bundle(bundle, frame):
    """Reloadable research predictor; preserves class order and missing categories."""
    probabilities=[]
    for member in bundle['models']:
        values=frame[member['columns']].copy()
        for column in member['categorical']:
            values[column]=values[column].fillna('__MISSING__').astype(str)
        classifier=member['classifier']
        order=[list(classifier.classes_).index(c) for c in bundle['classes']]
        probabilities.append(classifier.predict_proba(values)[:,order])
    return np.average(probabilities,axis=0,weights=bundle['weights'])


def run(output,modulus=10,iterations=400,resume=False):
    out=Path(output); out.mkdir(parents=True,exist_ok=resume)
    if resume:
        previous=__import__('json').loads((out/'protocol.json').read_text())
        if previous['modulus']!=modulus or previous['iterations']!=iterations:
            raise ValueError('Resume requires the original modulus and iteration count')
    cohort=ROOT/'data/processed/annotation_v3/full_cohort.sqlite'
    cache=ROOT/'data/processed/biological_features_20261029/features.sqlite'
    print('Loading development features',flush=True)
    x,y,cats=load_development(cohort,cache,'train',modulus)
    xv,yv,_=load_development(cohort,cache,'validation',modulus)
    for frame in (x,xv):
        for c in cats: frame[c]=frame[c].fillna('__MISSING__').astype(str)
    classes=np.unique(y)
    # Models fit training only; the existing locus-disjoint validation selects candidates.
    # Keep argmax fixed so no additional decision-tuning split is needed.
    profiles={'all':list(x),
              'population_sequence':[c for c in x if not c.startswith(('clingen_','inheritance_','brca1_','tp53_','assay_'))],
              'biological_context':[c for c in x if c.startswith(('clingen_','inheritance_','brca1_','tp53_','assay_'))]}
    configs=list(ParameterSampler({'depth':[5,7,9],'weight_power':[.5,1.,1.5],
                                   'l2_leaf_reg':[10,30],'learning_rate':[.03,.06]},n_iter=6,random_state=42))
    configs=[dict(depth=6,weight_power=0.,l2_leaf_reg=10,learning_rate=.05)]+configs
    results=[]; trained=[]; probs=[]
    save_json(out/'protocol.json',{'train_rows':len(y),'validation_rows':len(yv),'modulus':modulus,
        'iterations':iterations,'objective':['balanced_accuracy','macro_f1'],
        'test_rows_read':0,'selection_status':'reused development validation; not independent performance',
        'configurations':configs,'profiles':profiles,'decision_rule':'argmax; no threshold tuning'})
    for i,config in enumerate(configs+[configs[2],configs[2]]):
        profile='all' if i<len(configs) else list(profiles)[i-len(configs)+1]
        cols=profiles[profile]; cat=[c for c in cats if c in cols]
        counts=y.value_counts(); weights=np.array([(len(y)/(len(classes)*counts[v]))**config['weight_power'] for v in y])
        print(f'Training candidate {i}: {profile} {config}',flush=True)
        estimator=CatBoostClassifier(iterations=iterations,depth=config['depth'],learning_rate=config['learning_rate'],
            l2_leaf_reg=config['l2_leaf_reg'],thread_count=4,random_seed=42,verbose=False,allow_writing_files=False)
        saved=out/f'candidate_{i}.joblib'
        if resume and saved.exists():
            member=joblib.load(saved)
            if member['columns']!=cols:
                raise ValueError('Saved candidate profile differs; start a new run')
            estimator=member['classifier']
            print(f'Reusing completed candidate {i}',flush=True)
        else:
            estimator.fit(x[cols],y,cat_features=cat,sample_weight=weights)
        order=[list(estimator.classes_).index(c) for c in classes]
        p=estimator.predict_proba(xv[cols])[:,order]
        report=metrics(yv,classes[p.argmax(axis=1)],classes)
        results.append({'id':i,'profile':profile,'parameters':config,**report})
        trained.append({'classifier':estimator,'columns':cols,'categorical':cat}); probs.append(p)
        joblib.dump(trained[-1],out/f'candidate_{i}.joblib',compress=3)
        save_json(out/'development.json',results)
        print({k:report[k] for k in ('accuracy','balanced_accuracy','macro_f1')},flush=True)
    best_all=max(range(len(configs)),key=lambda i:(results[i]['balanced_accuracy'],results[i]['macro_f1']))
    members=[best_all,len(configs),len(configs)+1]
    mixtures=[]
    for weights in ((1,1,1),(2,1,1),(1,2,1),(3,1,0)):
        p=np.average([probs[i] for i in members],axis=0,weights=weights)
        mixtures.append({'members':members,'weights':list(weights),**metrics(yv,classes[p.argmax(axis=1)],classes)})
    save_json(out/'ensembles.json',mixtures)
    candidates=[('single',i,r) for i,r in enumerate(results)]+[('ensemble',i,r) for i,r in enumerate(mixtures)]
    kind,index,selected=max(candidates,key=lambda t:(t[2]['balanced_accuracy'],t[2]['macro_f1']))
    chosen=[index] if kind=='single' else members
    weights=[1] if kind=='single' else mixtures[index]['weights']
    joblib.dump({'models':[trained[i] for i in chosen],'weights':weights,'classes':classes,
                 'metadata':{'status':'development_candidate','test_evaluated':False,'selected':selected}},out/'selected.joblib',compress=3)
    save_json(out/'selection.json',{'kind':kind,'index':index,'selected':selected,'production_model_changed':False,
                                   'test_evaluated':False,'caution':'Selection-set metrics are optimistic; next step is frozen independent evaluation.'})
    print('Selected',kind,{k:selected[k] for k in ('accuracy','balanced_accuracy','macro_f1')},flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('--output',required=True)
    p.add_argument('--modulus',type=int,default=10); p.add_argument('--iterations',type=int,default=400)
    p.add_argument('--resume',action='store_true')
    a=p.parse_args()
    if a.modulus<1 or a.iterations<1:p.error('positive modulus/iterations required')
    run(a.output,a.modulus,a.iterations,a.resume)

