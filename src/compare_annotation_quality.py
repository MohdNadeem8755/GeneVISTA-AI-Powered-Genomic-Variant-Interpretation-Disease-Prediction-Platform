"""Review-quality experiment on training/validation only; no production promotion."""
import json
from pathlib import Path
import numpy as np
import joblib
from sklearn.ensemble import HistGradientBoostingClassifier
from src.annotation_pipeline import ROOT, connect
from src.audit_annotation_labels import quality
from src.train_annotation_candidates import matrix,metrics
from src.variant_features import CLASSES
out=ROOT/'reports/annotation_v3/quality_comparison.json'
if out.exists():
    raise ValueError('Experiment already exists')
base=joblib.load(ROOT/'reports/annotation_v2/full_candidates/validation_candidate.joblib')
features=base['features']
splits={}
with connect(ROOT/'data/processed/annotation_v2/download_verified_cohort.sqlite') as con:
    for split in ('train','validation'):
        raw=con.execute('''SELECT v.label,a.scores,v.gene,v.name,v.variant_type,v.review_status
            FROM variants v JOIN annotations a USING(chrom,pos,ref,alt) WHERE v.split=?
            ORDER BY v.chrom,v.pos,v.ref,v.alt''',(split,)).fetchall()
        rows=[r[:5] for r in raw]
        splits[split]=(matrix(rows,features),np.array([r[0] for r in rows]),
            np.array([quality(r[5]) in ('expert','multiple_no_conflicts') for r in raw]))
x,y,mask=splits['train']
xv,yv,mv=splits['validation']
counts={s:{'all':len(y),'high_quality':int(m.sum()),'high_quality_classes':{c:int(np.sum(y[m]==c)) for c in CLASSES}} for s,(_,y,m) in splits.items()}
if any(counts['train']['high_quality_classes'][c]<10 for c in CLASSES):
    raise ValueError('Insufficient high-quality five-class support')
results=[]
def evaluate(name,model):
    entry={'candidate':name,'validation':metrics(yv,model.predict_proba(xv),model.classes_),
           'high_quality_validation':metrics(yv[mv],model.predict_proba(xv[mv]),model.classes_)}
    results.append(entry)
    print(json.dumps({'candidate':name,'accuracy':entry['validation']['accuracy'],'macro_f1':entry['validation']['macro_f1'],'hq_macro_f1':entry['high_quality_validation']['macro_f1']}),flush=True)
evaluate('existing_annotation_candidate',base['classifier'])
for exponent in (0.,.5):
    _,class_counts=np.unique(y[mask],return_counts=True)
    counts_by_label=dict(zip(*np.unique(y[mask],return_counts=True)))
    weights=np.array([(mask.sum()/(5*counts_by_label[c]))**exponent for c in y[mask]])
    model=HistGradientBoostingClassifier(max_iter=250,max_leaf_nodes=15,min_samples_leaf=20,l2_regularization=10,learning_rate=.07,early_stopping=False,random_state=42)
    model.fit(x[mask],y[mask],sample_weight=weights)
    evaluate(f'high_quality_weight_{exponent}',model)
    joblib.dump({'classifier':model,'features':features},ROOT/f'reports/annotation_v3/quality_candidate_{exponent}.joblib',compress=3)
out.write_text(json.dumps({'status':'validation_only','counts':counts,'candidates':results,'test_rows_read':0,'calibration_rows_read':0,'production_model_changed':False,'scope':'Existing validation reused for exploratory comparison. No independent test claim.'},indent=2))
