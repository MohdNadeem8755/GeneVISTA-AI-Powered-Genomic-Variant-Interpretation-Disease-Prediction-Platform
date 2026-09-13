"""Exploratory specialists versus gene-aware pooled controls on identical validation rows."""
import json
from pathlib import Path
import joblib
import numpy as np
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.pipeline import Pipeline
from src.annotation_pipeline import ROOT, connect, meta
from src.train_annotation_candidates import matrix, metrics
from src.variant_features import CLASSES


def load_rows(con, genes, split):
    if split not in ('train', 'validation'):
        raise ValueError('Specialist experiments may read training/validation only')
    if not genes:
        raise ValueError('No eligible genes')
    marks=','.join('?' for _ in genes)
    return con.execute(f'''SELECT v.label,a.scores,v.gene,v.name,v.variant_type
        FROM variants v JOIN annotations a USING(chrom,pos,ref,alt)
        WHERE v.split=? AND v.gene IN ({marks})
        ORDER BY v.gene,v.chrom,v.pos,v.ref,v.alt''', (split,*genes)).fetchall()


def weights(labels, exponent):
    classes,counts=np.unique(labels,return_counts=True)
    if set(classes)!=set(CLASSES):
        raise ValueError('All five training classes required')
    lookup=dict(zip(classes,counts))
    return np.asarray([(len(labels)/(len(CLASSES)*lookup[label]))**exponent for label in labels])


class NonconstantFeatures(TransformerMixin, BaseEstimator):
    """Select variable columns from training only, retaining native missing values."""
    def fit(self, X, y=None):
        self.keep_ = np.array([len(np.unique(column[np.isfinite(column)])) > 1 for column in X.T])
        if not self.keep_.any():
            raise ValueError('No variable training features')
        return self
    def transform(self, X):
        return X[:, self.keep_]


def estimator():
    return Pipeline([('features', NonconstantFeatures()), ('model', HistGradientBoostingClassifier(max_iter=150,max_leaf_nodes=7,min_samples_leaf=15,
        l2_regularization=10,learning_rate=.07,early_stopping=False,random_state=42))])


def train(database, audit_path, output):
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    with (output/'started.json').open('x') as stream:
        json.dump({'scope':'validation_only','production_model_changed':False},stream)
    audit=json.loads(Path(audit_path).read_text())
    genes=sorted({row['gene'] for row in audit['eligible_genes'] if row['cohort']=='all'})
    incumbent=joblib.load(ROOT/'reports/annotation_v2/full_candidates/validation_candidate.joblib')
    features=incumbent['features']
    with connect(database) as con:
        if meta(con,'annotation_status')!='ready':raise ValueError('Import not ready')
        if meta(con,'annotation')['md5']!=incumbent['annotation']['md5']:
            raise ValueError('Annotation source differs from incumbent')
        rows={s:load_rows(con,genes,s) for s in ('train','validation')}
    x={s:matrix(r,features) for s,r in rows.items()}
    y={s:np.asarray([r[0] for r in values]) for s,values in rows.items()}
    g={s:np.asarray([r[2] for r in values]) for s,values in rows.items()}
    # Pooled controls know gene identity, matching the information specialists receive.
    pooled_x={s:np.column_stack([x[s],*(g[s]==gene for gene in genes)]) for s in rows}
    reports=[];aggregate=[]
    pooled_models={}
    for exponent in (0.,.5):
        model=estimator().fit(pooled_x['train'],y['train'],model__sample_weight=weights(y['train'],exponent))
        pooled_models[exponent]=model
        scores=model.predict_proba(pooled_x['validation'])
        aggregate.append({'candidate':f'pooled_{exponent}',**metrics(y['validation'],scores,model.classes_)})
    incumbent_scores=incumbent['classifier'].predict_proba(x['validation'])
    aggregate.append({'candidate':'previous_global_candidate',**metrics(y['validation'],incumbent_scores,incumbent['classifier'].classes_)})
    specialist_scores={e:np.empty_like(incumbent_scores) for e in (0.,.5)}
    artifacts={}
    for gene in genes:
        it=g['train']==gene;iv=g['validation']==gene
        entry={'gene':gene,'train_count':int(it.sum()),'validation_count':int(iv.sum()),'candidates':[]}
        entry['candidates'].append({'candidate':'previous_global_candidate',**metrics(y['validation'][iv],incumbent_scores[iv],incumbent['classifier'].classes_)})
        for exponent in (0.,.5):
            pooled=pooled_models[exponent]
            entry['candidates'].append({'candidate':f'pooled_{exponent}',**metrics(y['validation'][iv],pooled.predict_proba(pooled_x['validation'][iv]),pooled.classes_)})
            specialist=estimator().fit(x['train'][it],y['train'][it],model__sample_weight=weights(y['train'][it],exponent))
            if not np.array_equal(specialist.classes_,incumbent['classifier'].classes_):raise ValueError('Class ordering differs')
            scores=specialist.predict_proba(x['validation'][iv])
            specialist_scores[exponent][iv]=scores
            entry['candidates'].append({'candidate':f'specialist_{exponent}',**metrics(y['validation'][iv],scores,specialist.classes_)})
            artifacts[(gene,exponent)]=specialist
        reports.append(entry)
        print(json.dumps({'gene':gene,'results':[{k:r[k] for k in ('candidate','accuracy','macro_f1')} for r in entry['candidates']]}),flush=True)
    for exponent,scores in specialist_scores.items():
        aggregate.append({'candidate':f'specialist_{exponent}',**metrics(y['validation'],scores,incumbent['classifier'].classes_)})
    report={'status':'validation_complete','genes':genes,'rows':{s:len(r) for s,r in rows.items()},
        'aggregate':aggregate,'per_gene':reports,'test_rows_read':0,'calibration_rows_read':0,
        'production_model_changed':False,'scope':'Exploratory validation on eight support-selected genes, not overall catalog accuracy or an independent test.',
        'comparison':'Pooled and specialist models use identical expanded training rows and validation rows. Pooled models also receive gene identity. Previous global candidate used fewer training records.',
        'limitations':'Eligibility inspected validation label counts. Prior test was already evaluated; no new test claim. Source predictor overlap remains unresolved. No probability calibration or deployment promotion performed.'}
    joblib.dump({'specialists':artifacts,'pooled':pooled_models,'genes':genes,'features':features},output/'candidates.joblib',compress=3)
    (output/'comparison.json').write_text(json.dumps(report,indent=2))
    print(json.dumps({'aggregate':[{k:r[k] for k in ('candidate','accuracy','balanced_accuracy','macro_f1')} for r in aggregate]}),flush=True)
    return report


if __name__=='__main__':
    train(ROOT/'data/processed/annotation_v3/full_cohort.sqlite',
        ROOT/'reports/annotation_v3/full_cohort/label_audit.json', ROOT/'reports/annotation_v3/specialists_variable_features')

