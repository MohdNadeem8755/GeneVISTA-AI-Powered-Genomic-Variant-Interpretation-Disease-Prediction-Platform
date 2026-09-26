"""Exploratory temporal ranking of future ClinVar annotation events."""
from pathlib import Path
import re, json, gzip, csv, hashlib, warnings
from sklearn.exceptions import ConvergenceWarning
from src.feature_validation import load_label_manifest, load_labels
from src.training_validation import validate_manifest
from collections import Counter
import numpy as np
import pandas as pd
import joblib
from sklearn.feature_extraction import DictVectorizer
from sklearn.linear_model import SGDClassifier
from sklearn.metrics import average_precision_score

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"models/disease_ranking_v1"
CACHE=ROOT/"data/processed/disease_ranking_v1"
INTERVALS=[("2022-01","2023-01","train"),("2023-01","2024-01","train"),("2024-01","2025-01","validation"),("2025-01","2026-01","test")]

def input_signature():
    paths = [ROOT/'data/processed/model_labels/label_manifest.csv', ROOT/'data/processed/model_features/feature_manifest.csv', ROOT/'data/processed/clingen_gene_evidence_mapped.csv.gz']
    for start,end,split in INTERVALS:
        paths.extend([ROOT/f'data/processed/model_features/vus_features_{start}.csv.gz', ROOT/f'data/processed/model_labels/vus_labels_{start}_to_{end}.csv.gz'])
    paths.extend(ROOT/f'data/raw/clinvar/variant_summary_{snapshot}.txt.gz' for snapshot in sorted({x for start,end,_ in INTERVALS for x in (start,end)}))
    # Freshness guard: source path, size and nanosecond modification time.
    return {'version':2, 'inputs':{str(path.relative_to(ROOT)):[path.stat().st_size,path.stat().st_mtime_ns] for path in paths}}


def prepare():
    signature=input_signature()
    stamp=CACHE/'prepared_inputs.json'
    valid_cache=stamp.exists() and json.loads(stamp.read_text())==signature
    if valid_cache and (CACHE/'prepared.joblib').exists():
        print('Loading prepared historical cache', flush=True)
        return joblib.load(CACHE/'prepared.joblib')
    OUT.mkdir(parents=True,exist_ok=True);CACHE.mkdir(parents=True,exist_ok=True)
    label_manifest=load_label_manifest(ROOT/'data/processed/model_labels',ROOT/'data/processed')
    manifest=pd.read_csv(ROOT/'data/processed/model_features/feature_manifest.csv',dtype='string',keep_default_na=False)
    validate_manifest(manifest,label_manifest)
    expected=[(r.feature_snapshot,r.outcome_snapshot,r.split) for r in label_manifest.itertuples(index=False)]
    if expected!=INTERVALS:raise ValueError('Disease interval configuration differs from validated label manifest')
    frames={}; positives={}; diseases={}; names={}; genes={}
    for start,end,split in INTERVALS:
        f=pd.read_csv(ROOT/f"data/processed/model_features/vus_features_{start}.csv.gz",dtype='string',keep_default_na=False)
        row=next(r for r in label_manifest.itertuples(index=False) if r.feature_snapshot==start)
        labels=load_labels(ROOT/'data/processed/model_labels',row)
        if len(f)!=int(row.rows) or not f.feature_snapshot.eq(start).all():raise ValueError('Feature count/snapshot mismatch')
        f=f.merge(labels[["VariationID","feature_snapshot","target"]],on=["VariationID","feature_snapshot"],how='outer',validate="one_to_one",indicator=True)
        if not f.pop('_merge').eq('both').all():raise ValueError('Unmatched disease feature/label IDs')
        f["split"]=split;frames[start]=f
        positives[end]=set(f.loc[f.target.eq("became_pathogenic"),"VariationID"])
        print(start,len(f),"future pathogenic records:",len(positives[end]),flush=True)
    for snapshot in sorted(set(frames)|set(positives)):
        cache=CACHE/f"source_{snapshot}.joblib"
        if valid_cache and cache.exists():
            g,d,n=joblib.load(cache)
        else:
            g={};d={};n={};wanted=set(frames[snapshot].VariationID) if snapshot in frames else set()
            future=positives.get(snapshot,set())
            source=ROOT/f"data/raw/clinvar/variant_summary_{snapshot}.txt.gz"
            for chunk in pd.read_csv(source,sep="\t",usecols=["VariationID","GeneSymbol","PhenotypeIDS","PhenotypeList","Assembly"],dtype=str,chunksize=150000,keep_default_na=False):
                chunk=chunk.loc[chunk.Assembly.eq("GRCh38")]
                for row in chunk.loc[chunk.VariationID.isin(wanted|future)].itertuples(index=False):
                    if row.VariationID in wanted:g.setdefault(row.VariationID,row.GeneSymbol)
                    if row.VariationID in future:
                        ids=set(re.findall(r"MONDO:(?:MONDO:)?(\d+)",row.PhenotypeIDS))
                        d.setdefault(row.VariationID,set()).update("MONDO:"+value for value in ids)
                        groups=row.PhenotypeIDS.split("|");labels=row.PhenotypeList.split("|")
                        if len(groups)==len(labels):
                            for group,label in zip(groups,labels):
                                for value in re.findall(r"MONDO:(?:MONDO:)?(\d+)",group):n.setdefault("MONDO:"+value,label.replace("_"," "))
            joblib.dump((g,d,n),cache)
        genes[snapshot]=g;diseases[snapshot]=d;names.update(n)
        print("Read historical annotations:",snapshot,"genes",len(g),"labeled outcomes",len(d),flush=True)
    with gzip.open(ROOT/"data/processed/clingen_gene_evidence_mapped.csv.gz","rt",encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):names[row["mondo_id"]]=row["disease_name"]
    datasets={}
    for start,end,split in INTERVALS:
        f=frames[start];f["gene_symbol"]=f.VariationID.map(genes[start]).fillna("unknown")
        f["diseases"]=[sorted(diseases[end].get(value,set())) for value in f.VariationID]
        datasets.setdefault(split,[]).append(f)
    datasets={key:pd.concat(value,ignore_index=True) for key,value in datasets.items()}
    joblib.dump((datasets,names),CACHE/"prepared.joblib")
    return datasets,names

def feature_dict(row):
    def number(key):
        try:
            value = float(row.get(key,0) or 0)
            return float(np.log1p(max(0,value))) if np.isfinite(value) else 0.0
        except (TypeError,ValueError):return 0.0
    genes=[g.strip() for g in re.split(r"[;,|]",str(row.get("gene_symbol","unknown"))) if g.strip()]
    data={"gene="+gene:1.0 for gene in genes}
    data.update({"type="+str(row.get("variant_type","unknown")):1.0,"review="+str(row.get("review_status","unknown")):1.0,"log_submitters":number("number_submitters")})
    return data

def csr32(matrix):
    matrix = matrix.tocsr()
    if max(matrix.shape + (matrix.nnz,)) > np.iinfo(np.int32).max:
        raise ValueError('Feature matrix exceeds SGD sparse index limit')
    matrix.indices = matrix.indices.astype(np.int32, copy=False)
    matrix.indptr = matrix.indptr.astype(np.int32, copy=False)
    return matrix

def train(datasets,names,output_dir=None):
    output_dir=Path(output_dir) if output_dir is not None else OUT
    output_dir.mkdir(parents=True,exist_ok=True)
    for split in ('train','validation','test'):
        frame=datasets[split]
        if frame.empty or frame.duplicated(['VariationID','feature_snapshot']).any():raise ValueError('Empty or duplicate disease cohort')
        if not frame['split'].eq(split).all():raise ValueError('Incorrect disease split metadata')
        if any(not isinstance(labels,list) or len(labels)!=len(set(labels)) for labels in frame.diseases):raise ValueError('Invalid disease target lists')
        if any(labels and target!='became_pathogenic' for labels,target in zip(frame.diseases,frame.target)):raise ValueError('Disease target attached to nonpathogenic outcome')
    counts=Counter(value for labels in datasets["train"].diseases for value in labels)
    classes=[key for key,count in sorted(counts.items(),key=lambda item:(-item[1],item[0])) if 30<=count<len(datasets["train"])][:20]
    if not classes:raise ValueError("No disease has 30 positive training events; no model was trained.")
    vectorizer=DictVectorizer(dtype=np.float32)
    matrices={}
    for split in ('train', 'validation', 'test'):
        frame = datasets[split]
        records=(feature_dict(row) for row in (dict(zip(["gene_symbol","variant_type","review_status","number_submitters"], values)) for values in frame[["gene_symbol","variant_type","review_status","number_submitters"]].itertuples(index=False,name=None)))
        matrices[split]=vectorizer.fit_transform(records) if split=="train" else vectorizer.transform(records)
        matrices[split] = csr32(matrices[split])
        print("Feature matrix:",split,matrices[split].shape,flush=True)
    classifiers=[];metrics=[]
    for disease in classes:
        y=np.array([disease in row for row in datasets["train"].diseases],dtype=np.int8)
        model=SGDClassifier(loss="log_loss",alpha=0.00001,max_iter=1000,tol=0.001,class_weight="balanced",random_state=42,average=True)
        with warnings.catch_warnings():
            warnings.simplefilter('error',ConvergenceWarning)
            model.fit(matrices["train"],y)
        classifiers.append(model)
        for split in ["validation","test"]:
            truth=np.array([disease in row for row in datasets[split].diseases],dtype=np.int8)
            score=model.decision_function(matrices[split]);k=min(1000,len(truth))
            if not np.isfinite(score).all():raise ValueError('Nonfinite disease scores')
            tie_order=np.random.default_rng(42).permutation(len(truth))
            chosen=tie_order[np.argsort(-score[tie_order],kind="stable")[:k]]
            metrics.append({"split":split,"disease_id":disease,"disease_name":names.get(disease,disease),"train_positive_events":int(y.sum()),"positive_events":int(truth.sum()),"rows":len(truth),"prevalence":float(truth.mean()),"average_precision":float(average_precision_score(truth,score)) if truth.sum() else None,"review_size":k,"precision_at_1000":float(truth[chosen].mean()),"iterations":int(model.n_iter_)})
        print("Trained:",disease,names.get(disease),"positive events",int(y.sum()),flush=True)
    metadata={"target":"A starting VUS becomes pathogenic/likely pathogenic at the next snapshot AND a MONDO disease is listed in that future variant record. Record-level association does not establish disease-specific causality.","classes":classes,"names":{key:names.get(key,key) for key in classes},"training_counts":{key:counts[key] for key in classes},"score_type":"Uncalibrated SGD decision score; higher ranks first; not probability or clinical confidence.","class_selection":"Training-only support >=30 positive events; up to 20 most frequent MONDO IDs.","training_intervals":["2022-01 to 2023-01","2023-01 to 2024-01"],"validation_interval":"2024-01 to 2025-01","test_interval":"2025-01 to 2026-01","test_status":"Exploratory temporal evaluation. This time interval was already examined for earlier experiments; it is not a fresh external test.","limitations":["Unreported annotations are treated as unobserved events, not confirmed disease negatives.","Only modeled MONDO categories can be ranked; this is not all possible diseases.","Gene identity can dominate ranking; performance may not generalize to unseen genes or changing annotation practices.","Overlapping variants across snapshots mean observations are not independent.","The latest snapshot inference is an extrapolation beyond the training period.","No five-class clinical variant classifier or causal protein model is trained here."]}
    bundle={"vectorizer":vectorizer,"classifiers":classifiers,"metadata":metadata}
    joblib.dump(bundle,output_dir/"disease_ranker.joblib")
    (output_dir/"model_card.json").write_text(json.dumps(metadata,indent=2),encoding="utf-8")
    results=pd.DataFrame(metrics);results.to_csv(output_dir/"temporal_evaluation.csv",index=False)
    sample=matrices["validation"][:10];loaded=joblib.load(output_dir/"disease_ranker.joblib")
    assert loaded['metadata']==metadata and len(loaded['classifiers'])==len(classifiers)
    for original,restored in zip(classifiers,loaded['classifiers']):
        np.testing.assert_allclose(original.decision_function(sample),restored.decision_function(sample),rtol=1e-12,atol=1e-12)
    print("MODEL SAVED AND RELOAD VERIFIED",output_dir,flush=True)
    return results
