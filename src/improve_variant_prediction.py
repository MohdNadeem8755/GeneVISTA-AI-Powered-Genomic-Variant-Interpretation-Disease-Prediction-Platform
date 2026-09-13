"""Compare nonlinear candidates on a locus-disjoint inner validation split."""
import os,sys,json,hashlib,gc,argparse
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import numpy as np
import joblib
from sklearn.pipeline import make_pipeline
from sklearn.feature_selection import SelectKBest,chi2
from sklearn.preprocessing import FunctionTransformer
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.calibration import CalibratedClassifierCV
from sklearn.frozen import FrozenEstimator
from sklearn.metrics import accuracy_score,balanced_accuracy_score,classification_report,confusion_matrix,log_loss
from src.train_variant_prediction import prepare,estimator
from src.variant_features import variant_features,CLASSES
from src.nonlinear_prediction import dense_float32
from src.disease_ranking import csr32
OUT=ROOT/'models/variant_prediction_v1'
parser=argparse.ArgumentParser()
parser.add_argument('--iterations',type=int,default=150)
args=parser.parse_args()
if args.iterations<1:parser.error('--iterations must be positive')
frame,_,_=prepare()
model=joblib.load(OUT/'predictor.joblib')
train=frame[frame.split.eq('train')].reset_index(drop=True)
inner=np.array([int.from_bytes(hashlib.sha256(('inner-v2|'+s).encode()).digest()[:4],'big')%5==0 for s in train.site])
x=csr32(model['vectorizer'].transform(variant_features(r) for r in train.to_dict('records')))
y=train.label.to_numpy()
def scores(m,xx,yy):
    p=m.predict(xx)
    r=classification_report(yy,p,output_dict=True,zero_division=0)
    return {'accuracy':float(accuracy_score(yy,p)),'balanced_accuracy':float(balanced_accuracy_score(yy,p)),'macro_f1':r['macro avg']['f1-score']}
base=estimator(x[~inner],y[~inner])
results=[dict(name='linear',**scores(base,x[inner],y[inner]))]
print(results[-1],flush=True)
best=None;best_score=results[0]['macro_f1'];incumbent='linear'
previous=model['metadata'].get('training_experiments',{})
if args.iterations>150 and previous.get('selected'):
    matched=[r for r in previous['candidates'] if r['name']==previous['selected']]
    if matched:
        reference=dict(matched[0]);reference['name']='incumbent_'+reference['name'];reference['source']='previous measurement on identical inner validation split'
        results.append(reference);best_score=max(best_score,reference['macro_f1']);incumbent=reference['name']
        print('Incumbent comparison',reference,flush=True)
for leaves,power in ([(31,0.35)] if args.iterations>150 else [(15,0.0),(31,0.35)]):
    counts={label:int((y[~inner]==label).sum()) for label in np.unique(y)}
    weights={i:(len(y[~inner])/(5*counts[label]))**power for i,label in enumerate(sorted(counts))}
    candidate=make_pipeline(SelectKBest(chi2,k=256),FunctionTransformer(dense_float32),HistGradientBoostingClassifier(max_iter=args.iterations,verbose=1,max_leaf_nodes=leaves,learning_rate=.12,l2_regularization=5,min_samples_leaf=30,early_stopping=False,class_weight=weights,random_state=42))
    candidate.fit(x[~inner],y[~inner])
    result=dict(name=f'histogram_{leaves}_{power}_{args.iterations}_iterations',iterations=int(candidate[-1].n_iter_),**scores(candidate,x[inner],y[inner]))
    results.append(result);print(result,flush=True)
    if result['macro_f1']>best_score+.005:best=candidate;best_score=result['macro_f1']
    gc.collect()
report={'selection':'Candidates compared only on a deterministic locus-disjoint 20% inner validation partition of the original training split. Macro F1 selects a candidate if improvement exceeds 0.005. Original test was previously inspected for baseline, so final scores are reused internal test estimates, not fresh external validation.','inner_train_count':int((~inner).sum()),'inner_validation_count':int(inner.sum()),'candidates':results,'selected':incumbent if best is None else results[[r['macro_f1'] for r in results].index(best_score)]['name'],'target_accuracy':.99,'requested_iterations':args.iterations}
if best is not None:
    print('Refitting selected nonlinear model',flush=True)
    best.fit(x,y)
    del x;gc.collect()
    cal=frame[frame.split.eq('calibration')]
    xc=csr32(model['vectorizer'].transform(variant_features(r) for r in cal.to_dict('records')))
    clf=CalibratedClassifierCV(FrozenEstimator(best),method='sigmoid',ensemble=False).fit(xc,cal.label.to_numpy())
    del xc;gc.collect()
    test=frame[frame.split.eq('test')]
    xt=csr32(model['vectorizer'].transform(variant_features(r) for r in test.to_dict('records')))
    probs=clf.predict_proba(xt);pred=clf.classes_[probs.argmax(axis=1)];yt=test.label.to_numpy()
    old=model['metadata']['metrics']
    report['previous_test_metrics']={k:old[k] for k in ['accuracy','balanced_accuracy']}
    metrics=dict(old)
    metrics.update(accuracy=float(accuracy_score(yt,pred)),balanced_accuracy=float(balanced_accuracy_score(yt,pred)),log_loss=float(log_loss(yt,probs,labels=clf.classes_)),classification_report=classification_report(yt,pred,labels=CLASSES,output_dict=True,zero_division=0),confusion_matrix=confusion_matrix(yt,pred,labels=CLASSES).tolist())
    onehot=(yt[:,None]==clf.classes_[None,:]).astype(float)
    metrics['multiclass_brier']=float(np.mean(np.sum((probs-onehot)**2,axis=1)))
    reliability=[];ece=0.
    for i in range(10):
        mask=(probs.max(axis=1)>=i/10)&((probs.max(axis=1)<(i+1)/10) if i<9 else (probs.max(axis=1)<=1))
        if mask.any():
            confidence=float(probs[mask].max(axis=1).mean());accuracy=float((pred[mask]==yt[mask]).mean())
            ece+=float(mask.mean())*abs(confidence-accuracy)
            reliability.append(dict(lower=i/10,upper=(i+1)/10,count=int(mask.sum()),confidence=confidence,accuracy=accuracy))
    metrics.update(top_label_ece=ece,reliability=reliability)
    report['final_test_metrics']={k:metrics[k] for k in ['accuracy','balanced_accuracy']}
    report['final_test_metrics']['macro_f1']=metrics['classification_report']['macro avg']['f1-score']
    model['classifier']=clf;model['metadata']['metrics']=metrics
    model['metadata']['training_experiments']=report
    reuse_note=' The test has been reused after an inner-validation model comparison; these are internal estimates, not independent external validation.'
    if reuse_note not in model['metadata']['split']:model['metadata']['split']+=reuse_note
    limitation='The requested 99% accuracy has not been established. Per-class recall is substantially lower than aggregate accuracy.'
    if limitation not in model['metadata']['limitations']:model['metadata']['limitations'].append(limitation)
    joblib.dump(model,OUT/'predictor.joblib',compress=3)
    (OUT/'model_card.json').write_text(json.dumps(model['metadata'],indent=2),encoding='utf-8')
if best is None:
    report['final_test_metrics']={k:model['metadata']['metrics'][k] for k in ['accuracy','balanced_accuracy']}
    report['final_test_metrics']['macro_f1']=model['metadata']['metrics']['classification_report']['macro avg']['f1-score']
    report['predictor_updated']=False
else:report['predictor_updated']=True
report['achieved_99_percent']=model['metadata']['metrics']['accuracy']>=.99
(OUT/f'training_experiments_{args.iterations}.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
(OUT/'training_experiments.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps(report,indent=2),flush=True)
