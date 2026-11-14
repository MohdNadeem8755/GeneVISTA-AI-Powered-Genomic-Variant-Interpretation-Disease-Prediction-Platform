"""Isolated, development-only SMOTENC experiment; never deploys a model."""
import json
from pathlib import Path
import joblib
import numpy as np
import pandas as pd
from catboost import CatBoostClassifier
from imblearn.over_sampling import SMOTENC
from src.biological_features import load_development
from src.retune_biological_models import stopping_mask, passes_gate
from src.robust_variant_training import metrics, save_json, file_hash
from src.tune_biological_models import predict_bundle

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'reports/smotenc_trial_20261115_v3'

def run():
    OUT.mkdir(exist_ok=False)
    baseline_dir=ROOT/'reports/biological_corrected_20261109'
    baseline=json.loads((baseline_dir/'selection.json').read_text())
    digest=file_hash(baseline_dir/'selected.joblib')
    assert digest==baseline['model_sha256']
    save_json(OUT/'protocol.json',dict(seed=20261115,oversampling=1.5,
        targets=['Likely benign','Likely pathogenic'],neighbors=3,
        grouping='gene_identity',iterations=800,patience=100,depth=6,
        baseline_sha256=digest,test_rows_read=0,production_model_changed=False,
        limitation='Synthetic feature vectors are not observed variants or validated biological evidence.',
        acceptance='Existing strict no-decrease gate; both Likely recalls must improve.'))
    print('Loading fixed development partitions',flush=True)
    cohort=ROOT/'data/processed/annotation_v3/full_cohort.sqlite'
    cache=ROOT/'data/processed/biological_features_20261029/features.sqlite'
    x,y,cats,g=load_development(cohort,cache,'train',5,return_groups=True,feature_profile='variant_identity_v2')
    xv,yv,_,vg=load_development(cohort,cache,'validation',5,return_groups=True,feature_profile='variant_identity_v2')
    assert not set(g)&set(vg)
    stop=stopping_mask(g)
    old=joblib.load(baseline_dir/'selected.joblib')
    bp=predict_bundle(old,xv)
    reproduced=metrics(yv,old['classes'][bp.argmax(axis=1)],old['classes'])
    for key in ('accuracy','balanced_accuracy','macro_f1'):
        assert abs(reproduced[key]-baseline['selected'][key])<1e-10
    fit=x.loc[~stop].copy().reset_index(drop=True)
    fy=y[~stop].reset_index(drop=True)
    frames=[fit,x.loc[stop].copy().reset_index(drop=True),xv.copy()]
    # Parameters learned exclusively from fitting records; categorical flags are never interpolated.
    cats=list(cats)+[c for c in fit if c not in cats and fit[c].dropna().isin([0,1]).all()]
    numeric=[c for c in fit if c not in cats]
    med=fit[numeric].median().fillna(0)
    scale=fit[numeric].std().replace(0,1).fillna(1)
    missing=[c for c in numeric if fit[c].isna().any()]
    for frame in frames:
        for c in missing:
            frame['smote_missing__'+c]=frame[c].isna().astype(str)
        frame[numeric]=(frame[numeric].fillna(med)-med)/scale
        for c in cats:
            frame[c]=frame[c].fillna('__MISSING__').astype(str)
    cats += ['smote_missing__'+c for c in missing]
    assert len(cats)==len(set(cats)) and fit.columns.is_unique
    additions=[];labels=[]
    for gene,indices in fit.groupby('gene_identity',sort=True).groups.items():
        if gene=='__MISSING__':continue
        gx=fit.loc[indices]; gy=fy.loc[indices]; counts=gy.value_counts()
        if len(counts)<2:continue
        target={label:int(counts[label]*1.5) for label in ('Likely benign','Likely pathogenic') if counts.get(label,0)>=4}
        if not target:continue
        sampler=SMOTENC(categorical_features=cats,sampling_strategy=target,k_neighbors=3,random_state=20261115)
        sx,sy=sampler.fit_resample(gx,gy)
        additions.append(sx.iloc[len(gx):]);labels.append(sy.iloc[len(gy):])
    if not additions:raise ValueError('No eligible within-gene neighbors')
    synthetic=pd.concat(additions,ignore_index=True)
    synthetic_y=pd.concat(labels,ignore_index=True)
    # Restore missingness so interpolation does not turn unavailable evidence into measured evidence.
    for c in missing:
        synthetic.loc[synthetic['smote_missing__'+c]=='True',c]=np.nan
        for frame in frames:
            frame.loc[frame['smote_missing__'+c]=='True',c]=np.nan
    sx=pd.concat([fit,synthetic],ignore_index=True); sy=pd.concat([fy,synthetic_y],ignore_index=True)
    save_json(OUT/'data_audit.json',dict(fit_real=len(fy),synthetic=len(synthetic_y),
        synthetic_support=synthetic_y.value_counts().to_dict(),validation_rows=len(yv),
        stopping_rows=int(stop.sum()),locus_overlap=0,training_only_resampling=True))
    classes=np.unique(y); results=[]
    for name,tx,ty in [('control',fit,fy),('smotenc',sx,sy)]:
        print('Training',name,len(ty),flush=True)
        model=CatBoostClassifier(iterations=800,depth=6,learning_rate=.05,l2_leaf_reg=20,
            thread_count=16,random_seed=20261115,loss_function='MultiClass',
            eval_metric='TotalF1:average=Macro',allow_writing_files=False)
        model.fit(tx,ty,cat_features=cats,eval_set=(frames[1],y[stop]),early_stopping_rounds=100,verbose=100)
        prediction=model.predict(frames[2]).reshape(-1)
        result=metrics(yv,prediction,classes)
        result.update(name=name,iterations=model.tree_count_,passes_gate=passes_gate(result,baseline['selected']))
        results.append(result)
        joblib.dump(dict(classifier=model,columns=list(fit),categorical=cats,median=med,scale=scale,missing=missing),OUT/(name+'.joblib'),compress=3)
        save_json(OUT/'development.json',results)
    assert file_hash(baseline_dir/'selected.joblib')==digest
    save_json(OUT/'result.json',dict(status='complete',baseline=baseline['selected'],results=results,
        accepted=results[-1]['passes_gate'],production_model_changed=False,test_rows_read=0))
    print('COMPLETED',[(r['name'],r['accuracy'],r['balanced_accuracy'],r['passes_gate']) for r in results],flush=True)

if __name__=='__main__':run()
