"""Diagnose active-model class-boundary errors on validation only; no fitting."""
import argparse
from collections import Counter
import csv
from pathlib import Path
import sqlite3
import joblib
import numpy as np
import pandas as pd
from src.audit_annotation_labels import quality
from src.balanced_annotation import enhanced_matrix, decision_indices
from src.robust_variant_training import metrics, save_json, file_hash

ROOT = Path(__file__).resolve().parents[1]
PAIRS = {'Benign':'Likely benign', 'Likely benign':'Benign',
         'Pathogenic':'Likely pathogenic', 'Likely pathogenic':'Pathogenic'}


def error_kind(actual, predicted):
    if actual == predicted:
        return 'correct'
    if PAIRS.get(actual) == predicted:
        return 'within_pair'
    if predicted == 'Uncertain significance':
        return 'to_vus'
    if actual in PAIRS and predicted in PAIRS:
        return 'opposite_polarity'
    return 'other'


def groups(frame, keys):
    result = frame.groupby(keys, dropna=False).agg(
        support=('actual','size'), correct=('correct','sum'),
        within_pair=('within_pair','sum'), to_vus=('to_vus','sum'),
        opposite_polarity=('opposite_polarity','sum'),
        high_confidence_wrong=('high_confidence_wrong','sum'))
    result['recall'] = result.correct / result.support
    result['within_pair_rate'] = result.within_pair / result.support
    return result.reset_index()


def run(output):
    out = Path(output); out.mkdir(parents=True, exist_ok=False)
    artifact = ROOT/'models/annotation_prediction_v4/predictor.joblib'
    model = joblib.load(artifact)
    classes = model['classifier'].classes_
    queue = ROOT/'reports/label_quality_20261029/review_queue.csv'
    changed = {}
    with queue.open(newline='', encoding='utf-8') as stream:
        for row in csv.DictReader(stream):
            if row['split'] == 'validation':
                key = (row['variation_id'],row['chrom'],int(row['pos']),row['ref'],row['alt'])
                changed[key] = row
    con = sqlite3.connect((ROOT/'data/processed/annotation_v3/full_cohort.sqlite').as_uri()+'?mode=ro',uri=True)
    con.execute('ATTACH DATABASE ? AS bio',((ROOT/'data/processed/biological_features_20261029/features.sqlite').as_uri()+'?mode=ro',))
    training = {(g,c):n for g,c,n in con.execute("SELECT v.gene,v.label,COUNT(*) FROM variants v JOIN annotations a USING(chrom,pos,ref,alt) WHERE v.split='train' GROUP BY v.gene,v.label")}
    cursor = con.execute("SELECT v.label,a.scores,v.gene,v.name,v.variant_type,v.review_status,v.variation_id,v.chrom,v.pos,v.ref,v.alt,b.scores FROM variants v JOIN annotations a USING(chrom,pos,ref,alt) LEFT JOIN bio.features b USING(chrom,pos,ref,alt) WHERE v.split='validation' ORDER BY v.chrom,v.pos,v.ref,v.alt")
    import json
    records = []
    print('Scoring full matched validation set',flush=True)
    while rows := cursor.fetchmany(10000):
        x = enhanced_matrix([r[:5] for r in rows],model['features'],model['gene_vocab'])
        p = model['classifier'].predict_proba(x)
        selected = decision_indices(p,model['class_priors'],model['decision_exponent'])
        raw = p.argmax(axis=1)
        for i,r in enumerate(rows):
            actual,encoded,gene,name,kind,review,vid,chrom,pos,ref,alt,bio = r
            pred = str(classes[selected[i]])
            raw_pred = str(classes[raw[i]])
            flags = changed.get((vid,chrom,pos,ref,alt),{})
            reasons = flags.get('reason','')
            scores = json.loads(encoded); biology = json.loads(bio) if bio else {}
            rec = {'variation_id':vid,'chrom':chrom,'pos':pos,'ref':ref,'alt':alt,'gene':gene,
                   'actual':actual,'predicted':pred,'raw_argmax':raw_pred,'review_tier':quality(review),
                   'error_kind':error_kind(actual,pred),'correct':actual==pred,
                   'within_pair':error_kind(actual,pred)=='within_pair',
                   'to_vus':error_kind(actual,pred)=='to_vus',
                   'opposite_polarity':error_kind(actual,pred)=='opposite_polarity',
                   'raw_max_probability':float(p[i].max()),'selected_probability':float(p[i,selected[i]]),
                   'high_confidence_wrong':actual!=pred and p[i,selected[i]]>=.8,
                   'policy_changed_prediction':pred!=raw_pred,
                   'policy_corrected':pred==actual and raw_pred!=actual,
                   'policy_harmed':pred!=actual and raw_pred==actual,
                   'label_changed_since_snapshot':'changed_label' in reasons,
                   'latest_conflict':'latest_conflicting_classification' in reasons,
                   'old_evaluation':'old_evaluation_date' in reasons,
                   'latest_label':flags.get('latest_label',actual),
                   'review_flags':reasons,'training_gene_class_support':training.get((gene,actual),0),
                   'gene_in_model_vocabulary':gene in model['gene_vocab'],
                   'new_assay_available':biology.get('assay_match_available',0),
                   'new_clingen_available':biology.get('clingen_context_available',0)}
            for field in ('gnomAD4.1_joint_AF','REVEL_score','AlphaMissense_score','CADD_phred'):
                rec[field] = scores.get(field)
            records.append(rec)
        if len(records)%50000==0:
            print(f'Scored {len(records):,}',flush=True)
    con.close()
    frame = pd.DataFrame(records)
    frame.to_csv(out/'validation_predictions.csv.gz',index=False)
    focus = frame[frame.actual.isin(PAIRS)].copy()
    by_class = groups(focus,['actual'])
    by_class.to_csv(out/'class_errors.csv',index=False)
    groups(focus,['actual','review_tier']).to_csv(out/'review_tier_errors.csv',index=False)
    groups(focus,['actual','label_changed_since_snapshot']).to_csv(out/'label_drift_errors.csv',index=False)
    genes = groups(focus,['actual','gene'])
    genes['errors'] = genes.support-genes.correct
    genes.sort_values('errors',ascending=False).to_csv(out/'gene_errors.csv',index=False)
    for field in ('gnomAD4.1_joint_AF','REVEL_score','AlphaMissense_score','CADD_phred'):
        focus[field+'__missing'] = pd.to_numeric(focus[field],errors='coerce').isna()
        groups(focus,['actual',field+'__missing']).to_csv(out/(field+'_coverage_errors.csv'),index=False)
    numeric=['gnomAD4.1_joint_AF','REVEL_score','AlphaMissense_score','CADD_phred','training_gene_class_support','selected_probability']
    focus.groupby(['actual','error_kind'])[numeric].agg(['count','median']).to_csv(out/'feature_error_profiles.csv')
    errors=focus[~focus.correct].copy()
    # Bounded review list balances classes and destinations; not a new training set.
    examples=errors.sort_values(['latest_conflict','label_changed_since_snapshot','selected_probability'],ascending=False).groupby(['actual','error_kind']).head(20)
    examples.to_csv(out/'priority_examples.csv',index=False)
    report={'model_id':model['metadata']['model_id'],'model_sha256':file_hash(artifact),
            'scope':'Full annotation-matched validation, previously reused for model selection; descriptive analysis only.',
            'test_labels_read':0,'model_changed':False,'thresholds_tuned':False,
            'selected_policy':metrics(frame.actual,frame.predicted,classes),
            'raw_argmax_reference':metrics(frame.actual,frame.raw_argmax,classes),
            'class_errors':by_class.to_dict('records'),
            'policy_effects':frame.groupby('actual')[['policy_changed_prediction','policy_corrected','policy_harmed']].sum().reset_index().to_dict('records'),
            'focus_errors':len(errors),'within_pair_errors':int(errors.within_pair.sum()),
            'error_label_drift':{'changed':int(errors.label_changed_since_snapshot.sum()),'conflicting':int(errors.latest_conflict.sum())},
            'new_feature_coverage_among_errors':{'assay':int(errors.new_assay_available.sum()),'clingen':int(errors.new_clingen_available.sum())},
            'cautions':['Association is not causation; genes/review tiers are confounded with class and variant mix.',
                        'Current predictions use old features; new biological features are coverage diagnostics only.',
                        'High-confidence errors use selected-class probability >=0.8; decision-adjusted scores are not probabilities.',
                        'Missing queue entries mean no flagged audit issue, not independently verified label correctness.',
                        'No class merging, label replacement, test selection or model fitting performed.']}
    save_json(out/'report.json',report)
    print(json.dumps({'class_errors':report['class_errors'],'policy_effects':report['policy_effects']},indent=2),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',required=True)
    run(parser.parse_args().output)
