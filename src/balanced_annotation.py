"""Validation-selected, class-balanced annotation experiments with a sealed test extension."""
import argparse
import hashlib
import json
from pathlib import Path
import sqlite3
import time
import joblib
import numpy as np
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.calibration import CalibratedClassifierCV
from sklearn.frozen import FrozenEstimator
from sklearn.metrics import accuracy_score, balanced_accuracy_score, classification_report, log_loss
from src.annotation_pipeline import FEATURES, bucket, meta
from src.train_annotation_candidates import SEQUENCE, INDEPENDENT, matrix
from src.variant_features import CLASSES

ROOT=Path(__file__).resolve().parents[1]
ALL_FEATURES=FEATURES+SEQUENCE


def decision_indices(probabilities, priors, exponent):
    p=np.asarray(probabilities,dtype=float);prior=np.asarray(priors,dtype=float)
    if p.ndim!=2 or prior.shape!=(p.shape[1],) or not np.isfinite(p).all() or not np.isfinite(prior).all() or (prior<=0).any() or not 0<=exponent<=1:
        raise ValueError('Invalid probabilities, class priors or decision exponent')
    if (p<0).any() or (p>1).any() or not np.allclose(p.sum(axis=1),1):raise ValueError('Probabilities must sum to one')
    return (p / prior[None,:]**exponent).argmax(axis=1)


def evaluate(y,p,classes,priors=None,exponent=0):
    classes=np.asarray(classes)
    if set(y)!=set(CLASSES) or set(classes)!=set(CLASSES):raise ValueError('All five classes required for evaluation')
    priors=np.ones(len(classes))/len(classes) if priors is None else priors
    prediction=classes[decision_indices(p,priors,exponent)]
    report=classification_report(y,prediction,labels=CLASSES,output_dict=True,zero_division=0)
    return {'accuracy':float(accuracy_score(y,prediction)),'balanced_accuracy':float(balanced_accuracy_score(y,prediction)),
            'macro_f1':float(report['macro avg']['f1-score']),'log_loss':float(log_loss(y,p,labels=classes)),'per_class':report}


def enhanced_matrix(rows,features,gene_vocab):
    base=matrix(rows,features)
    genes=np.asarray([gene_vocab.get(r[2],np.nan) for r in rows],dtype=np.float32)
    return np.column_stack((base,genes)).astype(np.float32)


def read_split(con,split,features,gene_vocab,modulus=1,test_authorized=False):
    if split not in ('train','validation','calibration','test'):raise ValueError('Invalid split')
    if split=='test' and not test_authorized:raise ValueError('Test is sealed until selection is frozen')
    con.create_function('sample_group',2,lambda chrom,pos:bucket(f'balanced-v4-pilot|{chrom}|{pos}',modulus))
    con.create_function('legacy_sample',2,lambda chrom,pos:bucket(f'genevista-annotation-v2-sample|{chrom}|{pos}',10))
    extra=" AND legacy_sample(v.chrom,v.pos)<>0" if split=='test' else ''
    cursor=con.execute('''SELECT v.label,a.scores,v.gene,v.name,v.variant_type FROM variants v
        JOIN annotations a USING(chrom,pos,ref,alt) WHERE v.split=? AND sample_group(v.chrom,v.pos)=0'''+extra+' ORDER BY v.chrom,v.pos,v.ref,v.alt',(split,))
    arrays=[];labels=[]
    while rows:=cursor.fetchmany(25000):
        arrays.append(enhanced_matrix(rows,features,gene_vocab));labels.extend(r[0] for r in rows)
    if not arrays or set(labels)!=set(CLASSES):raise ValueError(f'Missing class support in {split}')
    return np.concatenate(arrays),np.asarray(labels)


def estimator(leaves,iterations,features):
    return HistGradientBoostingClassifier(max_iter=iterations,max_leaf_nodes=leaves,min_samples_leaf=30,
        l2_regularization=10,learning_rate=.07,early_stopping=False,random_state=42,
        categorical_features=[False]*len(features)+[True])


def fit_weighted(model,x,y,exponent):
    labels,counts=np.unique(y,return_counts=True);lookup=dict(zip(labels,counts))
    weights=np.asarray([(len(y)/(5*lookup[v]))**exponent for v in y])
    return model.fit(x,y,sample_weight=weights)


def run(database,output):
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    started=time.monotonic()
    def save(name,value): (output/name).write_text(json.dumps(value,indent=2,allow_nan=False),encoding='utf-8')
    with (output/'started.json').open('x') as stream:json.dump({'objective':'validation_balanced_accuracy','test_used_for_selection':False},stream)
    con=sqlite3.connect(Path(database).resolve().as_uri()+'?mode=ro',uri=True)
    try:
        if meta(con,'annotation_status')!='ready':raise ValueError('Verified annotations required')
        annotation=meta(con,'annotation')
        overlap=con.execute("SELECT COUNT(*) FROM (SELECT chrom,pos FROM variants WHERE split NOT IN ('excluded','reserved') GROUP BY chrom,pos HAVING COUNT(DISTINCT split)>1)").fetchone()[0]
        if overlap:raise ValueError('Position leakage across partitions')
        genes=con.execute("SELECT v.gene,COUNT(*) n FROM variants v JOIN annotations a USING(chrom,pos,ref,alt) WHERE v.split='train' AND v.gene NOT IN ('','-') GROUP BY v.gene ORDER BY n DESC,v.gene LIMIT 128").fetchall()
        vocab={g:i for i,(g,_) in enumerate(genes)}
        feature_audit=annotation['feature_nonmissing_counts']
        save('audit.json',{'position_overlap':overlap,'features':list(ALL_FEATURES),'feature_nonmissing_counts':feature_audit,'gene_vocabulary_source':'training only; most frequent 128 genes','gene_vocabulary':vocab,
             'excluded_inputs':['label','clinical significance','review status','disease annotations','VariationID'],
             'pretrained_deep_learning_feature':'AlphaMissense_score','upstream_training_overlap':'Possible overlap with upstream ClinVar-trained predictors; this is internal annotation prediction, not independent clinical validation.'})
        print('Preparing locus-sampled development matrices',flush=True)
        x,y=read_split(con,'train',ALL_FEATURES,vocab,5);xv,yv=read_split(con,'validation',ALL_FEATURES,vocab,5)
        candidates=[];best=None
        configs=[('all_unweighted',0.,15,ALL_FEATURES),('all_sqrt_weight',.5,15,ALL_FEATURES),('all_balanced',1.,15,ALL_FEATURES),('all_balanced_deeper',1.,31,ALL_FEATURES),('population_conservation_balanced',1.,15,INDEPENDENT+SEQUENCE)]
        for name,weight,leaves,features in configs:
            indices=[ALL_FEATURES.index(f) for f in features]+[len(ALL_FEATURES)]
            model=fit_weighted(estimator(leaves,160,features),x[:,indices],y,weight)
            result={'name':name,'weight':weight,'leaves':leaves,'features':list(features),**evaluate(yv,model.predict_proba(xv[:,indices]),model.classes_)}
            candidates.append(result)
            print(json.dumps({k:result[k] for k in ('name','accuracy','balanced_accuracy','macro_f1')}),flush=True)
            if best is None or tuple(result[k] for k in ('balanced_accuracy','macro_f1','accuracy'))>tuple(best[k] for k in ('balanced_accuracy','macro_f1','accuracy')):best=result
            save('development.json',{'status':'pilot','candidates':candidates,'selected':best,'train_rows':len(y),'validation_rows':len(yv),'test_rows_read':0})
        del x,y,xv,yv,model
        features=tuple(best['features'])
        print('Training selected configuration on all matched training loci',flush=True)
        x,y=read_split(con,'train',features,vocab)
        model=fit_weighted(estimator(best['leaves'],250,features),x,y,best['weight'])
        classes=model.classes_;prior=np.asarray([(y==c).mean() for c in classes]);training_count=len(y)
        del x,y
        xc,yc=read_split(con,'calibration',features,vocab)
        calibrated=CalibratedClassifierCV(FrozenEstimator(model),method='sigmoid').fit(xc,yc)
        calibration_count=len(yc);del xc,yc
        xv,yv=read_split(con,'validation',features,vocab)
        probabilities=calibrated.predict_proba(xv)
        policies=[{'exponent':e,**evaluate(yv,probabilities,classes,prior,e)} for e in (0.,.25,.5,.75,1.)]
        policy=max(policies,key=lambda r:(r['balanced_accuracy'],r['macro_f1'],r['accuracy']))
        validation_count=len(yv)
        print('Validation decision policies:',json.dumps([{k:r[k] for k in ('exponent','accuracy','balanced_accuracy','macro_f1')} for r in policies]),flush=True)
        candidate={'classifier':calibrated,'features':features,'gene_vocab':vocab,'class_priors':prior,'decision_exponent':policy['exponent']}
        joblib.dump(candidate,output/'candidate.joblib',compress=3)
        restored=joblib.load(output/'candidate.joblib')
        np.testing.assert_allclose(restored['classifier'].predict_proba(xv[:1000]),probabilities[:1000],rtol=1e-12,atol=1e-12)
        save('selection_frozen.json',{'candidate_sha256':hashlib.sha256((output/'candidate.joblib').read_bytes()).hexdigest(),'selected_configuration':best,'selected_policy':policy,'policies':policies,'selection_metric':'validation_balanced_accuracy','test_rows_read':0})
        del xv,yv,probabilities,restored
        print('Selection frozen; evaluating previously unscored test extension',flush=True)
        with (output/'test_started.json').open('x') as stream:json.dump({'selection_frozen':True},stream)
        xt,yt=read_split(con,'test',features,vocab,test_authorized=True)
        p=calibrated.predict_proba(xt);test=evaluate(yt,p,classes,prior,policy['exponent'])
        # Existing model evaluated on the identical extension, never for candidate selection.
        incumbent=joblib.load(ROOT/'models/annotation_prediction_v2/predictor.joblib')
        oldx,oldy=read_split(con,'test',tuple(incumbent['features']),{},test_authorized=True)
        np.testing.assert_array_equal(yt,oldy)
        previous=evaluate(yt,incumbent['classifier'].predict_proba(oldx[:,:-1]),incumbent['classifier'].classes_)
        del oldx,oldy
        report={'model_id':'annotation_prediction_v4','status':'research_evaluated','annotation':annotation,'features':list(features),'test':test,'incumbent_same_test':previous,
            'test_count':len(yt),'training_count':training_count,'validation_count':validation_count,'calibration_count':calibration_count,
            'test_total_eligible':con.execute("SELECT COUNT(*) FROM variants WHERE split='test' AND legacy_sample(chrom,pos)<>0").fetchone()[0],
            'validation_selection':policy,'decision_exponent':policy['exponent'],'decision_policy':'Class-balanced decision using calibrated probability divided by training prior to the selected power. Predicted class may differ from the highest raw probability.',
            'scope':'Internal GRCh38 annotation-matched SNVs. Same genes can occur across splits; not external or clinical validation.',
            'upstream_overlap':'Pretrained annotation predictors may overlap ClinVar labels; no clinical accuracy or patient disease-risk claim.',
            'test_policy':'Original v2 sampled cohort excluded from this test extension; extension was not used for model/policy selection. Future tuning must not reuse these scores as an untouched test.',
            'pretrained_dl_used':'AlphaMissense_score' in features,'split_position_overlap':overlap,'production_model_changed':False,'seconds':round(time.monotonic()-started,1)}
        candidate['metadata']=report
        joblib.dump(candidate,output/'predictor.joblib',compress=3)
        save('model_card.json',report)
        print(json.dumps({'test':{k:test[k] for k in ('accuracy','balanced_accuracy','macro_f1')},'incumbent_same_test':{k:previous[k] for k in ('accuracy','balanced_accuracy','macro_f1')}}),flush=True)
        return report
    finally:con.close()

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--database',type=Path,default=ROOT/'data/processed/annotation_v3/full_cohort.sqlite');parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();run(args.database,args.output)
