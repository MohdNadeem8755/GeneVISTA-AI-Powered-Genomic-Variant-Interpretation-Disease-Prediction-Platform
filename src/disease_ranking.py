"""Exploratory temporal ranking of future ClinVar annotation events."""
from pathlib import Path
import re, json, gzip, csv
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

def prepare():
    if (CACHE/'prepared.joblib').exists():
        print('Loading prepared historical cache', flush=True)
        return joblib.load(CACHE/'prepared.joblib')
    OUT.mkdir(parents=True,exist_ok=True);CACHE.mkdir(parents=True,exist_ok=True)
    frames={}; positives={}; diseases={}; names={}; genes={}
    for start,end,split in INTERVALS:
        f=pd.read_csv(ROOT/f"data/processed/model_features/vus_features_{start}.csv.gz",dtype={"VariationID":str})
        labels=pd.read_csv(ROOT/f"data/processed/model_labels/vus_labels_{start}_to_{end}.csv.gz",dtype={"VariationID":str})
        f=f.merge(labels[["VariationID","target"]],on="VariationID",validate="one_to_one")
        f["split"]=split;frames[start]=f
        positives[end]=set(f.loc[f.target.eq("became_pathogenic"),"VariationID"])
        print(start,len(f),"future pathogenic records:",len(positives[end]),flush=True)
    for snapshot in sorted(set(frames)|set(positives)):
        cache=CACHE/f"source_{snapshot}.joblib"
        if cache.exists():
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

def train(datasets,names):
    OUT.mkdir(parents=True,exist_ok=True)
    counts=Counter(value for labels in datasets["train"].diseases for value in labels)
    classes=[key for key,count in sorted(counts.items(),key=lambda item:(-item[1],item[0])) if count>=30][:20]
    if not classes:raise ValueError("No disease has 30 positive training events; no model was trained.")
    vectorizer=DictVectorizer(dtype=np.float32)
    matrices={}
    for split in ('train', 'validation', 'test'):
        frame = datasets[split]
        records=(feature_dict(row) for row in frame[["gene_symbol","variant_type","review_status","number_submitters"]].to_dict("records"))
        matrices[split]=vectorizer.fit_transform(records) if split=="train" else vectorizer.transform(records)
        matrices[split] = csr32(matrices[split])
        print("Feature matrix:",split,matrices[split].shape,flush=True)
    classifiers=[];metrics=[]
    for disease in classes:
        y=np.array([disease in row for row in datasets["train"].diseases],dtype=np.int8)
        model=SGDClassifier(loss="log_loss",alpha=0.00001,max_iter=30,tol=0.001,class_weight="balanced",random_state=42,average=True)
        model.fit(matrices["train"],y);classifiers.append(model)
        for split in ["validation","test"]:
            truth=np.array([disease in row for row in datasets[split].diseases],dtype=np.int8)
            score=model.decision_function(matrices[split]);k=min(1000,len(truth));chosen=np.argsort(-score,kind="stable")[:k]
            metrics.append({"split":split,"disease_id":disease,"disease_name":names.get(disease,disease),"train_positive_events":int(y.sum()),"positive_events":int(truth.sum()),"rows":len(truth),"prevalence":float(truth.mean()),"average_precision":float(average_precision_score(truth,score)) if truth.sum() else None,"review_size":k,"precision_at_1000":float(truth[chosen].mean()),"iterations":int(model.n_iter_)})
        print("Trained:",disease,names.get(disease),"positive events",int(y.sum()),flush=True)
    metadata={"target":"A starting VUS becomes pathogenic/likely pathogenic at the next snapshot AND a MONDO disease is listed in that future variant record. Record-level association does not establish disease-specific causality.","classes":classes,"names":{key:names.get(key,key) for key in classes},"training_counts":{key:counts[key] for key in classes},"score_type":"Uncalibrated SGD decision score; higher ranks first; not probability or clinical confidence.","class_selection":"Training-only support >=30 positive events; up to 20 most frequent MONDO IDs.","training_intervals":["2022-01 to 2023-01","2023-01 to 2024-01"],"validation_interval":"2024-01 to 2025-01","test_interval":"2025-01 to 2026-01","test_status":"Exploratory temporal evaluation. This time interval was already examined for earlier experiments; it is not a fresh external test.","limitations":["Unreported annotations are treated as unobserved events, not confirmed disease negatives.","Only modeled MONDO categories can be ranked; this is not all possible diseases.","Gene identity can dominate ranking; performance may not generalize to unseen genes or changing annotation practices.","Overlapping variants across snapshots mean observations are not independent.","The latest snapshot inference is an extrapolation beyond the training period.","No five-class clinical variant classifier or causal protein model is trained here."]}
    bundle={"vectorizer":vectorizer,"classifiers":classifiers,"metadata":metadata}
    joblib.dump(bundle,OUT/"disease_ranker.joblib")
    (OUT/"model_card.json").write_text(json.dumps(metadata,indent=2),encoding="utf-8")
    results=pd.DataFrame(metrics);results.to_csv(OUT/"temporal_evaluation.csv",index=False)
    sample=matrices["validation"][:10];loaded=joblib.load(OUT/"disease_ranker.joblib")
    np.testing.assert_allclose(classifiers[0].decision_function(sample),loaded["classifiers"][0].decision_function(sample))
    print("MODEL SAVED AND RELOAD VERIFIED",OUT,flush=True)
    return results
