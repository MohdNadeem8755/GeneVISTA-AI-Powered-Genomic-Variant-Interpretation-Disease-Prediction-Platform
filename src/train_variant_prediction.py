"""Train exploratory five-class and disease-annotation models using cached ClinVar."""
from pathlib import Path
import csv,gzip,hashlib,json,re,sys
from collections import Counter,defaultdict
import numpy as np
import pandas as pd
import joblib
from sklearn.feature_extraction import DictVectorizer
from sklearn.linear_model import SGDClassifier
from sklearn.calibration import CalibratedClassifierCV
from sklearn.frozen import FrozenEstimator
from sklearn.metrics import accuracy_score,balanced_accuracy_score,classification_report,confusion_matrix,log_loss,average_precision_score,brier_score_loss

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from src.variant_features import variant_features,CLASSES
from src.disease_ranking import csr32
OUT=ROOT/'models/variant_prediction_v1'
CACHE=ROOT/'data/processed/variant_prediction_v1'
SOURCE=ROOT/'data/processed/clinvar_clean_2026-09.csv.gz'
SCHEMA=1

def site(row):return '|'.join(str(row[k]) for k in ('Chromosome','Start','Stop'))

def prepare():
    CACHE.mkdir(parents=True,exist_ok=True)
    cached=CACHE/'cohort.joblib'
    if cached.exists():return joblib.load(cached)
    columns=['VariationID','Type','Name','GeneSymbol','ClinicalSignificance','PhenotypeIDS','PhenotypeList','Chromosome','Start','Stop','OriginSimple']
    reserved=set();found=set()
    for chunk in pd.read_csv(SOURCE,usecols=columns,dtype=str,keep_default_na=False,chunksize=100000):
        for row in chunk[chunk.VariationID.isin(['17660','37565'])].to_dict('records'):
            reserved.add(site(row));found.add(row['VariationID'])
        if len(found)==2:break
    assert len(found)==2,'Reserved evaluation examples not found'
    records=[];names={};seen=set();audits=Counter();split_sites=defaultdict(set)
    for chunk in pd.read_csv(SOURCE,usecols=columns,dtype=str,keep_default_na=False,chunksize=150000):
        audits['source_rows']+=len(chunk)
        for row in chunk.to_dict('records'):
            if row['ClinicalSignificance'] not in CLASSES:audits['excluded_non_exact_label']+=1;continue
            if row['OriginSimple'] not in ('germline','unknown'):audits['excluded_non_germline']+=1;continue
            if not row['Start'].isdigit() or int(row['Start'])<=0:audits['excluded_no_locus']+=1;continue
            key=site(row);digest=hashlib.sha256(('genevista-v1|'+key).encode()).digest()
            if key not in reserved and int.from_bytes(digest[:4],'big')%10!=0:continue
            if row['VariationID'] in seen:continue
            seen.add(row['VariationID'])
            bucket=int.from_bytes(digest[4:8],'big')%100
            split='test' if key in reserved or bucket>=85 else 'calibration' if bucket>=70 else 'train'
            split_sites[split].add(key)
            diseases=sorted(set('MONDO:'+x for x in re.findall(r'MONDO:(?:MONDO:)?(\d+)',row['PhenotypeIDS'])))
            groups=row['PhenotypeIDS'].split('|');labels=row['PhenotypeList'].split('|')
            if len(groups)==len(labels):
                for group,label in zip(groups,labels):
                    for value in re.findall(r'MONDO:(?:MONDO:)?(\d+)',group):names.setdefault('MONDO:'+value,label.replace('_',' '))
            records.append({'variation_id':int(row['VariationID']),'gene_symbol':row['GeneSymbol'],'name':row['Name'],
                'variant_type':row['Type'],'label':row['ClinicalSignificance'],'diseases':diseases,'split':split,'site':key})
        print('Prepared',audits['source_rows'],'source rows;',len(records),'sampled records',flush=True)
    for a,b in [('train','calibration'),('train','test'),('calibration','test')]:assert not split_sites[a]&split_sites[b]
    frame=pd.DataFrame(records)
    with gzip.open(ROOT/'data/processed/clingen_gene_evidence_mapped.csv.gz','rt',encoding='utf-8-sig') as f:
        for row in csv.DictReader(f):names[row['mondo_id']]=row['disease_name']
    data=(frame,names,dict(audits))
    joblib.dump(data,cached,compress=3)
    return data

def estimator(x,y,seed=42):
    model=SGDClassifier(loss='log_loss',alpha=0.00001,max_iter=60,tol=.001,average=True,random_state=seed)
    model.fit(x,y)
    return model

def fit_calibrated(x,y,xc,yc):
    model=estimator(x,y)
    calibrated=CalibratedClassifierCV(FrozenEstimator(model),method='sigmoid',ensemble=False)
    calibrated.fit(xc,yc)
    return calibrated

def train():
    OUT.mkdir(parents=True,exist_ok=True)
    frame,names,audits=prepare()
    if (OUT/'disease_names.json').exists():
        names.update(json.loads((OUT/'disease_names.json').read_text(encoding='utf-8')))
    splits={key:frame[frame.split.eq(key)].reset_index(drop=True) for key in ['train','calibration','test']}
    vectorizer=DictVectorizer(dtype=np.float32)
    matrices={}
    for split,part in splits.items():
        records=(variant_features(row) for row in part.to_dict('records'))
        matrices[split]=csr32(vectorizer.fit_transform(records) if split=='train' else vectorizer.transform(records))
        print('Features',split,matrices[split].shape,part.label.value_counts().to_dict(),flush=True)
    y={split:part.label.to_numpy() for split,part in splits.items()}
    classifier=fit_calibrated(matrices['train'],y['train'],matrices['calibration'],y['calibration'])
    probabilities=classifier.predict_proba(matrices['test']);predictions=classifier.classes_[probabilities.argmax(axis=1)]
    onehot=np.array([[int(value==label) for label in classifier.classes_] for value in y['test']])
    train_counts=Counter(y['train']);total=len(y['train'])
    prior=np.array([train_counts[label]/total for label in classifier.classes_]);baseline=np.tile(prior,(len(y['test']),1))
    ece=0.0;reliability=[]
    for low in np.arange(0,1,.1):
        mask=(probabilities.max(axis=1)>=low)&(probabilities.max(axis=1)<low+.1+1e-12)
        if not mask.any():continue
        confidence=float(probabilities[mask].max(axis=1).mean());accuracy=float((predictions[mask]==y['test'][mask]).mean())
        ece+=float(mask.mean())*abs(confidence-accuracy)
        reliability.append({'lower':float(low),'upper':float(min(1,low+.1)),'count':int(mask.sum()),'confidence':confidence,'accuracy':accuracy})
    metrics={'accuracy':float(accuracy_score(y['test'],predictions)),'balanced_accuracy':float(balanced_accuracy_score(y['test'],predictions)),
        'log_loss':float(log_loss(y['test'],probabilities,labels=classifier.classes_)),
        'prior_baseline_log_loss':float(log_loss(y['test'],baseline,labels=classifier.classes_)),
        'multiclass_brier':float(np.mean(np.sum((probabilities-onehot)**2,axis=1))),
        'top_label_ece':ece,'reliability':reliability,
        'classification_report':classification_report(y['test'],predictions,labels=CLASSES,output_dict=True,zero_division=0),
        'confusion_matrix':confusion_matrix(y['test'],predictions,labels=CLASSES).tolist()}
    print('CLASSIFICATION METRICS',json.dumps({k:v for k,v in metrics.items() if isinstance(v,float)}),flush=True)
    eligible={split:part.diseases.map(bool).to_numpy() for split,part in splits.items()}
    counts=Counter(d for labels in splits['train'].loc[eligible['train'],'diseases'] for d in labels)
    disease_ids=[d for d,count in counts.most_common() if count>=200][:30]
    disease_models=[];disease_metrics=[];genes_by_disease={}
    for disease in disease_ids:
        truth={split:np.array([int(disease in ds) for ds in part.loc[eligible[split],'diseases']],dtype=np.int8) for split,part in splits.items()}
        if min(truth['calibration'].sum(),len(truth['calibration'])-truth['calibration'].sum())<20:continue
        model=fit_calibrated(matrices['train'][eligible['train']],truth['train'],matrices['calibration'][eligible['calibration']],truth['calibration'])
        scores=model.predict_proba(matrices['test'][eligible['test']])[:,1];yt=truth['test']
        gene_counts=Counter(g.strip() for gene in splits['train'].loc[splits['train'].diseases.map(lambda labels:disease in labels),'gene_symbol'] for g in re.split(r'[;,|]',gene) if g.strip() and g.strip()!='-')
        genes_by_disease[disease]=sorted(gene for gene,count in gene_counts.items() if count>=3)
        disease_models.append((disease,model))
        disease_metrics.append({'disease_id':disease,'disease_name':names.get(disease,disease),'train_positive':int(truth['train'].sum()),'test_positive':int(yt.sum()),'test_rows':len(yt),
            'average_precision':float(average_precision_score(yt,scores)) if yt.sum() else None,'prevalence':float(yt.mean()),'brier':float(brier_score_loss(yt,scores))})
        print('Disease calibrated',disease,names.get(disease),flush=True)
    meta={'schema_version':SCHEMA,'model_id':'variant_prediction_v1','source_snapshot':'2026-09',
        'target':'Predict the exact recorded ClinVar class among five eligible labels from gene and parsed sequence/protein-change features. This is an exploratory annotation classifier, not a clinical diagnosis.',
        'classes':CLASSES,'input_fields':['gene_symbol','variant_type','name'],
        'excluded_inputs':['clinical_significance','clinical_state','review_status','number_submitters_raw','phenotype','disease labels','VariationID','patient information'],
        'sampling':'Deterministic approximately 10% sample by GRCh38 chromosome/start/stop locus. Exact five-label germline/unknown-origin records only; combined/conflicting labels excluded.',
        'split':'70% training, 15% sigmoid calibration, 15% internal test, disjoint genomic loci. The loci containing examples 17660 and 37565 are reserved for test. Same genes can occur across splits; this is not external or unseen-gene validation.',
        'split_counts':{key:len(part) for key,part in splits.items()},'label_counts':{key:part.label.value_counts().to_dict() for key,part in splits.items()},'preparation_audit':audits,
        'disease_target':'Independent estimates that a MONDO disease annotation is listed, conditional on a record having at least one MONDO annotation. Missing annotations are not confirmed disease negatives. These are not patient disease probabilities and need not sum to 100%.',
        'disease_candidate_rule':'Up to 30 training-supported MONDO classes (at least 200 records); at least 20 calibration positives/negatives. At inference show only diseases linked to an input gene in at least three training annotations.',
        'limitations':['Limited cached features omit population allele frequency, functional assays, segregation, zygosity, and patient phenotype.',
          'Five-class label probabilities do not implement ACMG/AMP evidence criteria and are not evidence of clinical pathogenicity.',
          'A base substitution alone does not specify a variant; gene, position/transcript, and other evidence matter.',
          'Internal calibration does not establish clinical validity or patient penetrance.',
          'Disease annotation matches can reflect investigation or reporting patterns, including for benign variants; they do not establish causality.'],
        'metrics':metrics,'disease_metrics':disease_metrics}
    saved={'vectorizer':vectorizer,'classifier':classifier,'disease_models':disease_models,'genes_by_disease':genes_by_disease,'disease_names':names,'metadata':meta,
        'cohort_ids':{split:part.variation_id.to_numpy(dtype=np.int64) for split,part in splits.items()}}
    joblib.dump(saved,OUT/'predictor.joblib',compress=3)
    (OUT/'model_card.json').write_text(json.dumps(meta,indent=2),encoding='utf-8')
    loaded=joblib.load(OUT/'predictor.joblib')
    np.testing.assert_allclose(classifier.predict_proba(matrices['test'][:20]),loaded['classifier'].predict_proba(matrices['test'][:20]))
    print('PREDICTOR SAVED AND RELOAD VERIFIED',OUT,flush=True)
    return meta

if __name__=='__main__':train()
