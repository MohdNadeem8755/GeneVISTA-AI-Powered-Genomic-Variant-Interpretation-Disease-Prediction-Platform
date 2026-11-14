"""One-shot evaluation of a frozen candidate on novel variants from a later snapshot.

Workflow: freeze -> prepare -> annotate -> evaluate. Each stage writes an immutable
marker. Labels are first read only after freeze. Evaluation refuses to run twice.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import json
import math
from pathlib import Path
import sqlite3

import joblib
import numpy as np
import pandas as pd

from src.annotation_pipeline import FEATURES, snv_key
from src.biological_features import clingen_index, enrich, gene_context
from src.robust_variant_training import file_hash, metrics, save_json
from src.train_annotation_candidates import INDEPENDENT, SEQUENCE, matrix
from src.tune_biological_models import predict_bundle
from src.variant_features import CLASSES

ROOT=Path(__file__).resolve().parents[1]
LATER_URL='https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz'
MD5_URL=LATER_URL+'.md5'
HIGH_REVIEW={'reviewed by expert panel','practice guideline','criteria provided, multiple submitters, no conflicts'}


def freeze(output,snapshot='2026-10-26'):
    out=Path(output); out.mkdir(parents=True,exist_ok=False)
    candidate=ROOT/'reports/biological_tuning_20261030/selected.joblib'
    selection=ROOT/'reports/biological_tuning_20261030/selection.json'
    if not candidate.is_file() or not selection.is_file():raise FileNotFoundError('Frozen Step 4 candidate missing')
    chosen=json.loads(selection.read_text())
    if chosen.get('test_evaluated') or chosen.get('production_model_changed'):
        raise ValueError('Candidate is not an unevaluated development artifact')
    protocol={'schema':1,'status':'frozen_before_labels','frozen_utc':datetime.now(timezone.utc).isoformat(),
      'snapshot':snapshot,'source_url':LATER_URL,'checksum_url':MD5_URL,
      'candidate':str(candidate),'candidate_sha256':file_hash(candidate),
      'selection_sha256':file_hash(selection),'baseline_snapshot':'2026-10-05',
      'eligibility':['GRCh38','germline or unknown origin','single-nucleotide substitution','exactly one of five target labels',
        'unambiguous label and identity per exact allele','VariationID absent from October 5 snapshot',
        'allele absent from October 5 snapshot','allele absent from original September model cohort','exact dbNSFP allele match'],
      'primary_cohort':'All review tiers meeting eligibility; no confidence threshold or abstention.',
      'secondary_cohort':'Expert/practice-guideline or multiple-submitters-no-conflicts; reported only when all five classes occur.',
      'primary_objective':['balanced_accuracy','macro_f1'],'reported_metrics':['accuracy','balanced_accuracy','macro precision','macro recall','macro F1','weighted metrics','per-class metrics','confusion matrix'],
      'decision_rule':'Frozen Step 4 soft-voting argmax. No threshold, weight, feature, class, or cohort tuning after labels.',
      'minimum_support':'Primary requires every class and at least 30 records per class; otherwise report insufficient evidence.',
      'comparison':'Descriptive later-snapshot concordance only; no model selection or production replacement.',
      'test_read_count':0}
    save_json(out/'frozen_protocol.json',protocol)
    print(json.dumps({'status':'frozen','candidate_sha256':protocol['candidate_sha256']}))


def verified_source(out):
    protocol=json.loads((out/'frozen_protocol.json').read_text())
    source=out/f"variant_summary_{protocol['snapshot']}.txt.gz"
    checksum=out/(source.name+'.md5')
    if not source.is_file() or not checksum.is_file():raise FileNotFoundError('Download source and checksum after freezing')
    token=next((x for x in checksum.read_text().split() if len(x)==32 and all(c in '0123456789abcdefABCDEF' for c in x)),None)
    if not token:raise ValueError('No MD5 in checksum file')
    with source.open('rb') as stream:actual=hashlib.file_digest(stream,'md5').hexdigest()
    if actual!=token.lower():raise ValueError('Later snapshot MD5 mismatch')
    return protocol,source,actual


def prepare(output):
    out=Path(output); protocol,source,digest=verified_source(out)
    final=out/'cohort.sqlite'; staging=out/'cohort.building.sqlite'
    if final.exists():raise FileExistsError('Prepared cohort already exists')
    resuming=staging.exists();db=sqlite3.connect(staging); db.execute('PRAGMA journal_mode=OFF'); db.execute('PRAGMA synchronous=OFF')
    counts=Counter();baseline=ROOT/'data/raw/clinvar/variant_summary_2026-10-05.txt.gz'
    if not resuming:
      db.executescript('CREATE TABLE later(chrom TEXT,pos INTEGER,ref TEXT,alt TEXT,variation_id TEXT,gene TEXT,name TEXT,variant_type TEXT,label TEXT,review TEXT); CREATE TABLE baseline_id(id TEXT PRIMARY KEY); CREATE TABLE baseline_allele(chrom TEXT,pos INTEGER,ref TEXT,alt TEXT,PRIMARY KEY(chrom,pos,ref,alt));')
      exact={x.lower():x for x in CLASSES}; batch=[]
      print('Reading frozen later snapshot',flush=True)
      with gzip.open(source,'rt',encoding='utf-8-sig',newline='') as stream:
        reader=csv.DictReader(stream,delimiter='\t')
        for row in reader:
          counts['later_source_rows']+=1
          if row.get('OriginSimple') not in ('germline','unknown'):continue
          key=snv_key(row); label=exact.get(row.get('ClinicalSignificance','').strip().lower())
          if not key or not label:continue
          batch.append((*key,row['VariationID'],row['GeneSymbol'],row['Name'],row['Type'],label,row['ReviewStatus']))
          if len(batch)>=20000:db.executemany('INSERT INTO later VALUES (?,?,?,?,?,?,?,?,?,?)',batch);db.commit();batch.clear()
      db.executemany('INSERT INTO later VALUES (?,?,?,?,?,?,?,?,?,?)',batch);db.commit()
      batch_id=[];batch_key=[]
      print('Indexing pre-test October 5 snapshot',flush=True)
      with gzip.open(baseline,'rt',encoding='utf-8-sig',newline='') as stream:
        for row in csv.DictReader(stream,delimiter='\t'):
          counts['baseline_source_rows']+=1; batch_id.append((row['VariationID'],))
          key=snv_key(row)
          if key:batch_key.append(key)
          if len(batch_id)>=25000:
            db.executemany('INSERT OR IGNORE INTO baseline_id VALUES (?)',batch_id);batch_id.clear()
            db.executemany('INSERT OR IGNORE INTO baseline_allele VALUES (?,?,?,?)',batch_key);batch_key.clear();db.commit()
      db.executemany('INSERT OR IGNORE INTO baseline_id VALUES (?)',batch_id);db.executemany('INSERT OR IGNORE INTO baseline_allele VALUES (?,?,?,?)',batch_key);db.commit()
    else:
      required={'later','baseline_id','baseline_allele'}
      present={r[0] for r in db.execute("SELECT name FROM sqlite_master WHERE type='table'")}
      if not required.issubset(present):raise ValueError('Incomplete staged cohort cannot be resumed')
      counts['later_eligible_source_rows']=db.execute('SELECT COUNT(*) FROM later').fetchone()[0]
      counts['baseline_distinct_ids']=db.execute('SELECT COUNT(*) FROM baseline_id').fetchone()[0]
      counts['baseline_distinct_snv_alleles']=db.execute('SELECT COUNT(*) FROM baseline_allele').fetchone()[0]
      print('Resuming indexed cohort finalization',flush=True)
    db.execute('ATTACH DATABASE ? AS original',(str((ROOT/'data/processed/annotation_v3/full_cohort.sqlite').resolve()),))
    db.executescript('''CREATE TABLE eligible AS WITH grouped AS (
      SELECT chrom,pos,ref,alt,MIN(variation_id) variation_id,MIN(gene) gene,MIN(name) name,MIN(variant_type) variant_type,
      MIN(label) label,MIN(review) review,COUNT(DISTINCT label) labels,COUNT(DISTINCT variation_id) ids
      FROM later GROUP BY chrom,pos,ref,alt)
      SELECT g.* FROM grouped g WHERE labels=1 AND ids=1
      AND NOT EXISTS(SELECT 1 FROM baseline_id b WHERE b.id=g.variation_id)
      AND NOT EXISTS(SELECT 1 FROM baseline_allele b WHERE (b.chrom,b.pos,b.ref,b.alt)=(g.chrom,g.pos,g.ref,g.alt))
      AND NOT EXISTS(SELECT 1 FROM original.variants v WHERE (v.chrom,v.pos,v.ref,v.alt)=(g.chrom,g.pos,g.ref,g.alt));
      CREATE UNIQUE INDEX eligible_allele ON eligible(chrom,pos,ref,alt);
      DROP TABLE later; DROP TABLE baseline_id; DROP TABLE baseline_allele;''')
    support=dict(db.execute('SELECT label,COUNT(*) FROM eligible GROUP BY label').fetchall())
    db.commit();db.close();staging.replace(final)
    report={'status':'prepared_labels_sealed','snapshot':protocol['snapshot'],'source_md5':digest,'baseline_sha256':file_hash(baseline),'counts':dict(counts),'eligible_support_before_annotation':support,'candidate_unchanged':file_hash(protocol['candidate'])==protocol['candidate_sha256']}
    save_json(out/'prepare_report.json',report);print(json.dumps(report,indent=2))


def annotate(output):
    """Extract the six non-clinical dbNSFP fields using its existing tabix index."""
    import pysam
    out=Path(output); protocol=json.loads((out/'frozen_protocol.json').read_text())
    if (out/'annotations_complete.json').exists():raise FileExistsError('Annotations already complete')
    db=sqlite3.connect(out/'cohort.sqlite'); targets=db.execute('SELECT chrom,pos,ref,alt FROM eligible ORDER BY chrom,pos').fetchall()
    if not targets:raise ValueError('No eligible later variants')
    source=ROOT/'data/incoming/dbNSFP5.4a_grch38.gz'; tabix=pysam.TabixFile(str(source))
    header=next(x for x in tabix.header if x.startswith('#chr\t')).split('\t'); index={x:i for i,x in enumerate(header)}
    fields=list(INDEPENDENT); required=['#chr','pos(1-based)','ref','alt']+fields
    if any(x not in index for x in required):raise ValueError('dbNSFP field missing')
    db.execute('CREATE TABLE annotations(chrom TEXT,pos INTEGER,ref TEXT,alt TEXT,scores TEXT,PRIMARY KEY(chrom,pos,ref,alt))')
    matched=Counter();batch=[]
    print(f'Querying dbNSFP for {len(targets):,} frozen targets',flush=True)
    for n,(chrom,pos,ref,alt) in enumerate(targets,1):
      values=defaultdict(list)
      try: rows=tabix.fetch(chrom,pos-1,pos)
      except ValueError:rows=[]
      for line in rows:
        parts=line.split('\t')
        if parts[index['pos(1-based)']]!=str(pos) or parts[index['ref']]!=ref or parts[index['alt']]!=alt:continue
        matched['source_rows']+=1
        for field in fields:
          for token in parts[index[field]].replace(';',',').split(','):
            try:value=float(token)
            except ValueError:continue
            if math.isfinite(value):values[field].append(value)
      if values:
        # Match original annotation aggregation: maximum finite transcript/source value.
        batch.append((chrom,pos,ref,alt,json.dumps({f:max(values[f]) if values[f] else None for f in fields},allow_nan=False)))
      if len(batch)>=5000:db.executemany('INSERT INTO annotations VALUES (?,?,?,?,?)',batch);db.commit();batch.clear()
      if n%10000==0:print(f'Queried {n:,}; matched {matched["source_rows"]:,} source rows',flush=True)
    db.executemany('INSERT INTO annotations VALUES (?,?,?,?,?)',batch);db.commit();tabix.close()
    observed=db.execute('SELECT COUNT(*) FROM annotations').fetchone()[0];support=dict(db.execute('SELECT e.label,COUNT(*) FROM eligible e JOIN annotations a USING(chrom,pos,ref,alt) GROUP BY e.label').fetchall());db.close()
    report={'status':'ready','target_variants':len(targets),'matched_variants':observed,'matched_source_rows':matched['source_rows'],'support':support,'features':fields,'dbnsfp_sha256':file_hash(source),'candidate_unchanged':file_hash(protocol['candidate'])==protocol['candidate_sha256']}
    save_json(out/'annotations_complete.json',report);print(json.dumps(report,indent=2))


def evaluate(output):
    out=Path(output);protocol=json.loads((out/'frozen_protocol.json').read_text())
    if (out/'final_results.json').exists() or (out/'evaluation_started.json').exists():raise FileExistsError('One-shot evaluation already started')
    if file_hash(protocol['candidate'])!=protocol['candidate_sha256']:raise ValueError('Frozen candidate changed')
    prepared=json.loads((out/'prepare_report.json').read_text())
    support=prepared['eligible_support_before_annotation'];missing=[c for c in CLASSES if int(support.get(c,0))<30]
    if missing:
      result={'status':'insufficient_support','snapshot':protocol['snapshot'],'classes_below_30':missing,'support':support,
        'minimum_required_per_class':30,'candidate_sha256':protocol['candidate_sha256'],'candidate_unchanged':True,
        'predictions_generated':False,'performance_metrics_computed':False,'test_reused_for_selection':False,
        'production_model_changed':False,'reason':'The pre-registered minimum class support was not met. The cohort definition and threshold were not changed after support was observed.',
        'next_valid_action':'Wait for a later pre-registered snapshot to accumulate enough novel examples in every class; do not repeatedly score this sealed cohort.'}
      save_json(out/'insufficient_support.json',result);print(json.dumps(result,indent=2));return
    annotation=json.loads((out/'annotations_complete.json').read_text())
    save_json(out/'evaluation_started.json',{'started_utc':datetime.now(timezone.utc).isoformat(),'candidate_sha256':protocol['candidate_sha256']})
    brca=json.load(gzip.open(ROOT/'data/processed/brca1_mavedb_v1.json.gz','rt'));tp53=json.load(gzip.open(ROOT/'data/processed/tp53_evidence_v1.json.gz','rt'))
    context_rows=clingen_index(ROOT/'data/processed/clingen_gene_evidence_mapped.csv.gz',datetime.strptime('2026-10-26','%Y-%m-%d').date());contexts={g:gene_context(g,context_rows) for g in context_rows}
    feature_meta=json.loads((ROOT/'data/processed/biological_features_20261029/report.json').read_text())
    base_names=tuple(INDEPENDENT)+SEQUENCE
    db=sqlite3.connect(out/'cohort.sqlite');rows=db.execute('SELECT e.label,a.scores,e.gene,e.name,e.variant_type,e.review,e.variation_id,e.chrom,e.pos,e.ref,e.alt FROM eligible e JOIN annotations a USING(chrom,pos,ref,alt) ORDER BY e.chrom,e.pos,e.ref,e.alt').fetchall();db.close()
    base=pd.DataFrame(matrix([r[:5] for r in rows],base_names),columns=base_names)
    added=pd.DataFrame([enrich(r[2],r[3],json.loads(r[1]),contexts,brca,tp53) for r in rows]).reindex(columns=feature_meta['numeric']+feature_meta['categorical'])
    base=base.drop(columns=[n for n in added if n in base],errors='ignore');x=pd.concat([base,added],axis=1)
    bundle=joblib.load(protocol['candidate']);p=predict_bundle(bundle,x);pred=bundle['classes'][p.argmax(axis=1)];y=np.asarray([r[0] for r in rows]);reviews=np.asarray([r[5] for r in rows])
    primary=metrics(y,pred,bundle['classes']);mask=np.array([r in HIGH_REVIEW for r in reviews]);secondary=None
    if set(y[mask])==set(CLASSES):secondary=metrics(y[mask],pred[mask],bundle['classes'])
    result={'status':'evaluated_once','snapshot':protocol['snapshot'],'evaluation_utc':datetime.now(timezone.utc).isoformat(),'candidate_sha256':protocol['candidate_sha256'],'primary':primary,'high_review_secondary':secondary,'high_review_support':dict(Counter(y[mask])),'feature_rows':len(x),'test_reused_for_selection':False,'production_model_changed':False,'target_accuracy_95_met':primary['accuracy']>=.95,'target_balanced_accuracy_90_met':primary['balanced_accuracy']>=.90,'limitations':['Later-snapshot concordance, not clinical validation.','Functional and gene-level evidence may also contribute to ClinVar assertions.','Only novel variants with exact dbNSFP matches are covered; this is a selected subset of the later release.','No confidence filtering, threshold tuning, retraining or post-label model selection was performed.']}
    save_json(out/'final_results.json',result);pd.DataFrame({'variation_id':[r[6] for r in rows],'actual':y,'predicted':pred,'review':[r[5] for r in rows]}).to_csv(out/'predictions.csv.gz',index=False)
    print(json.dumps({'primary':{k:primary[k] for k in ('n','accuracy','balanced_accuracy','macro_f1')},'high_review_n':int(mask.sum())},indent=2))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('stage',choices=['freeze','prepare','annotate','evaluate']);p.add_argument('--output',required=True);p.add_argument('--snapshot',default='2026-10-26');a=p.parse_args()
    {'freeze':lambda:freeze(a.output,a.snapshot),'prepare':lambda:prepare(a.output),'annotate':lambda:annotate(a.output),'evaluate':lambda:evaluate(a.output)}[a.stage]()


if __name__=='__main__':main()
