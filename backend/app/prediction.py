"""Research predictions from sequence-derived features, with explicit model scope."""
from functools import lru_cache
from pathlib import Path
import re
import joblib
import numpy as np
from src.variant_features import variant_features,parse_change,CLASSES
from src.disease_ranking import csr32

MODEL_DIR=Path(__file__).resolve().parents[2]/'models/variant_prediction_v1'
@lru_cache(maxsize=1)
def bundle():
    result=joblib.load(MODEL_DIR/'predictor.joblib')
    result['cohort_ids']={key:np.sort(values) for key,values in result['cohort_ids'].items()}
    return result

def membership(model,variant_id):
    if variant_id is None:return 'unlisted_input'
    for split,ids in model['cohort_ids'].items():
        index=np.searchsorted(ids,variant_id)
        if index<len(ids) and ids[index]==variant_id:return split
    return 'outside_sampled_cohort'

def predict(details,variant_id=None,include_diseases=True):
    if not details:
        return {'status':'insufficient_evidence','reason':'Sequence and gene details are missing for this variant.'}
    if not (MODEL_DIR/'predictor.joblib').exists():
        return {'status':'model_unavailable','reason':'The five-class research model has not been installed.'}
    model=bundle();features=variant_features(details)
    known_gene=any(key.startswith('gene=') and key in model['vectorizer'].vocabulary_ for key in features)
    if not known_gene:
        return {'status':'insufficient_evidence','reason':'This gene was not observed during model training. More evidence is needed.'}
    x=csr32(model['vectorizer'].transform([features]))
    scores=model['classifier'].predict_proba(x)[0]
    probabilities={str(label):float(value) for label,value in zip(model['classifier'].classes_,scores)}
    predicted=max(probabilities,key=probabilities.get)
    observed=details.get('clinical_significance')
    disagreement=observed in CLASSES and predicted!=observed
    change=parse_change(details.get('name',''))
    warnings=[]
    if probabilities[predicted]<.7:warnings.append('The model does not strongly separate the possible classes. Collect additional evidence.')
    if predicted=='Uncertain significance':warnings.append('The top estimated class is VUS; pathogenicity is unresolved.')
    if disagreement:warnings.append('The model result differs from the recorded ClinVar classification. Do not override the recorded evidence with this estimate.')
    if not change['nucleotide'] and not change['protein']:warnings.append('No single nucleotide or protein substitution could be parsed; this estimate relies on gene, variant type, and event features.')
    result={'status':'available','model_id':'variant_prediction_v1','predicted_class':predicted,
        'probabilities':[{'class_name':label,'probability':probabilities[label]} for label in CLASSES],
        'top_probability':probabilities[predicted],'change':change,'cohort_membership':membership(model,variant_id),
        'disagrees_with_recorded':disagreement,'warnings':warnings,
        'next_step':'Review the original evidence with a qualified genetics professional. These research estimates do not establish a diagnosis.',
        'probability_scope':'Sigmoid-calibrated model estimates for exact ClinVar annotation classes; not patient disease risk.',
        'evaluation_scope':model['metadata']['split']}
    if include_diseases:
        genes={g.strip() for g in re.split(r'[;,|]',str(details.get('gene_symbol',''))) if g.strip()}
        rows=[];metrics={row['disease_id']:row for row in model['metadata']['disease_metrics']}
        for disease,classifier in model['disease_models']:
            if not genes.intersection(model['genes_by_disease'][disease]):continue
            probability=float(classifier.predict_proba(x)[0,1])
            metric=metrics[disease]
            rows.append({'disease_id':disease,'disease_name':model['disease_names'].get(disease,disease),
                'match_probability':probability,'training_records':metric['train_positive'],
                'test_positive':metric['test_positive'],'test_average_precision':metric['average_precision']})
        rows.sort(key=lambda row:(-row['match_probability'],row['disease_id']))
        for i,row in enumerate(rows):row['rank']=i+1
        result['diseases']={'status':'available' if rows else 'insufficient_evidence','items':rows,
            'reason':None if rows else 'No sufficiently supported disease category for this gene is included in the model. Consult the curated associations below.',
            'probability_scope':model['metadata']['disease_target']}
    return result
