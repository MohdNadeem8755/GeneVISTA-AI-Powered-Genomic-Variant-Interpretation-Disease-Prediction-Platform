"""Offline GRCh38 SNV cohort, streamed dbNSFP merge, validation-only comparison."""
import argparse
from collections import Counter
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
import re
import sqlite3
import time
from src.variant_features import CLASSES

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DB = ROOT / 'data/processed/annotation_v2/cohort.sqlite'
FEATURES = ('gnomAD4.1_joint_AF', 'gnomAD4.1_joint_POPMAX_AF',
    'gnomAD4.1_joint_SAS_AF', 'gnomAD4.1_joint_AN', 'CADD_phred',
    'REVEL_score', 'SIFT_converted_rankscore', 'Polyphen2_HVAR_rankscore',
    'AlphaMissense_score', 'GERP_92_mammals', 'phyloP100way_vertebrate')

def bucket(text, modulus):
    return int.from_bytes(hashlib.sha256(text.encode()).digest()[:4], 'big') % modulus

class ClosingConnection(sqlite3.Connection):
    def __exit__(self, *args):
        try:
            return super().__exit__(*args)
        finally:
            self.close()

def connect(path):
    con = sqlite3.connect(path, factory=ClosingConnection)
    con.execute('PRAGMA cache_size=-16384')
    con.execute('PRAGMA temp_store=MEMORY')
    return con

def put_meta(con, key, value):
    con.execute('INSERT OR REPLACE INTO metadata VALUES (?,?)', (key, json.dumps(value)))

def meta(con, key):
    row = con.execute('SELECT value FROM metadata WHERE key=?', (key,)).fetchone()
    return json.loads(row[0]) if row else None

def snv_key(row):
    if row.get('Assembly') != 'GRCh38': return None
    chrom, pos = row.get('Chromosome', '').removeprefix('chr'), row.get('PositionVCF', '')
    ref, alt = row.get('ReferenceAlleleVCF', ''), row.get('AlternateAlleleVCF', '')
    if chrom not in [str(i) for i in range(1, 23)] + ['X', 'Y']: return None
    if not pos.isdigit() or int(pos) <= 0 or ref not in ('A','C','G','T') or alt not in ('A','C','G','T') or ref == alt: return None
    return chrom, int(pos), ref, alt

def assign_split(chrom, pos):
    n = bucket(f'genevista-annotation-v2-split|{chrom}|{pos}', 100)
    split = 'train' if n < 65 else 'validation' if n < 80 else 'calibration' if n < 90 else 'test'
    old_sample = bucket(f'genevista-v1|{chrom}|{pos}|{pos}', 10) == 0
    return 'train' if split == 'test' and old_sample else split

def prepare(source, output, sample_modulus=10):
    if sample_modulus < 1: raise ValueError('sample-modulus must be positive')
    output, source = Path(output), Path(source)
    staging = output.with_suffix('.building')
    if output.exists() or staging.exists(): raise ValueError('Output/staging already exists; choose a new output path')
    output.parent.mkdir(parents=True, exist_ok=True)
    counts = Counter()
    con = connect(staging)
    try:
        con.executescript('''CREATE TABLE metadata(key TEXT PRIMARY KEY,value TEXT NOT NULL);
        CREATE TABLE variants(chrom TEXT,pos INTEGER,ref TEXT,alt TEXT,variation_id TEXT,
          gene TEXT,name TEXT,variant_type TEXT,label TEXT,review_status TEXT,split TEXT,
          conflict INTEGER NOT NULL DEFAULT 0,PRIMARY KEY(chrom,pos,ref,alt));
        CREATE TABLE annotations(chrom TEXT,pos INTEGER,ref TEXT,alt TEXT,
          scores TEXT NOT NULL,source_rows INTEGER NOT NULL,PRIMARY KEY(chrom,pos,ref,alt));''')
        statement = '''INSERT INTO variants VALUES (?,?,?,?,?,?,?,?,?,?,?,0)
        ON CONFLICT(chrom,pos,ref,alt) DO UPDATE SET
          conflict=MAX(variants.conflict,variants.label<>excluded.label,variants.variation_id<>excluded.variation_id),
          split=CASE WHEN excluded.split='reserved' THEN 'reserved' ELSE variants.split END'''
        with gzip.open(source, 'rt', encoding='utf-8-sig', newline='') as stream:
            for row in csv.DictReader(stream, delimiter='\t'):
                counts['source_rows'] += 1
                if row.get('Assembly') != 'GRCh38':
                    counts['other_assembly'] += 1; continue
                if row.get('OriginSimple') not in ('germline','unknown'):
                    counts['other_origin'] += 1; continue
                key = snv_key(row)
                if key is None:
                    counts['unsupported_coordinates_or_variant_type'] += 1; continue
                chrom, pos, ref, alt = key
                if bucket(f'genevista-annotation-v2-sample|{chrom}|{pos}', sample_modulus): continue
                split = 'reserved' if row['VariationID'] in ('17660','37565') else assign_split(chrom, pos)
                con.execute(statement, (*key, row['VariationID'], row['GeneSymbol'], row['Name'], row['Type'], row['ClinicalSignificance'], row['ReviewStatus'], split))
                counts['sampled_source_rows'] += 1
                if counts['sampled_source_rows'] % 20000 == 0:
                    con.commit(); print(f"Prepared {counts['sampled_source_rows']:,} sampled rows", flush=True)
        con.execute('''UPDATE variants SET conflict=1 WHERE variation_id IN
          (SELECT variation_id FROM variants GROUP BY variation_id HAVING COUNT(*)>1)''')
        placeholders = ','.join('?' for _ in CLASSES)
        con.execute(f"UPDATE variants SET split='excluded' WHERE conflict=1 OR label NOT IN ({placeholders})", CLASSES)
        con.execute("UPDATE variants SET split='reserved' WHERE split<>'excluded' AND (chrom,pos) IN (SELECT chrom,pos FROM variants WHERE variation_id IN ('17660','37565'))")
        con.executescript('CREATE INDEX variant_split ON variants(split); CREATE INDEX variant_id ON variants(variation_id);')
        overlap = con.execute("SELECT COUNT(*) FROM (SELECT chrom,pos FROM variants WHERE split NOT IN ('excluded','reserved') GROUP BY chrom,pos HAVING COUNT(DISTINCT split)>1)").fetchone()[0]
        if overlap: raise ValueError('Genomic-position split overlap detected')
        audit = {'schema':1, 'assembly':'GRCh38', 'scope':'Autosomal/X/Y SNVs only; not indels, CNVs or all ClinVar variants.',
            'source':{'name':source.name,'size_bytes':source.stat().st_size}, 'sample_modulus':sample_modulus,
            'counts':dict(counts), 'split_label_counts':[list(r) for r in con.execute('SELECT split,label,COUNT(*) FROM variants GROUP BY split,label ORDER BY split,label')],
            'conflicting_keys':con.execute('SELECT COUNT(*) FROM variants WHERE conflict=1').fetchone()[0],
            'duplicate_source_rows_collapsed':counts['sampled_source_rows']-con.execute('SELECT COUNT(*) FROM variants').fetchone()[0],
            'split_position_overlap':overlap,
            'split_policy':'Genomic-position groups; nominal 65/15/10/10 train/validation/calibration/test. v1 sampled SNV loci cannot enter test; published examples reserved. Internal same-snapshot validation only.',
            'test_policy':'compare never reads test or calibration rows. Upstream predictor training overlap still needs audit.'}
        put_meta(con, 'cohort', audit); con.commit()
    finally: con.close()
    staging.replace(output)
    output.with_suffix('.audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
    return audit

def check_md5(source, checksum):
    match = re.search(r'\b[0-9a-fA-F]{32}\b', Path(checksum).read_text(encoding='utf-8-sig'))
    if not match: raise ValueError('Checksum file contains no MD5 digest')
    with Path(source).open('rb') as stream: actual = hashlib.file_digest(stream,'md5').hexdigest()
    if actual != match[0].lower(): raise ValueError('MD5 mismatch: incomplete/corrupt download; no annotations imported')
    return actual

def number(value):
    values = []
    for token in re.split('[;,]',value):
        try:
            n = float(token)
            if math.isfinite(n): values.append(n)
        except ValueError: pass
    return max(values) if values else None

def merge(source, checksum, database):
    source, database = Path(source), Path(database)
    if 'grch38' not in source.name.lower() or '5.4a' not in source.name: raise ValueError('Expected versioned GRCh38 dbNSFP 5.4a filename')
    if not database.is_file(): raise ValueError('Prepare the cohort first')
    verified_md5 = check_md5(source,checksum)
    started = time.monotonic()
    con = connect(database)
    try:
        if meta(con,'annotation_status') == 'ready': raise ValueError('Annotations already ready; use a new cohort for a different source')
        con.execute('DELETE FROM annotations'); put_meta(con,'annotation_status','importing'); con.commit()
        targets = set(con.execute("SELECT chrom,pos,ref,alt FROM variants WHERE split<>'excluded'"))
        seen = matched = 0
        with gzip.open(source,'rt',encoding='utf-8-sig') as stream:
            for line in stream:
                if line.startswith('#chr\t'):
                    header = line.rstrip('\r\n').split('\t'); break
            else: raise ValueError('dbNSFP #chr header missing')
            positions = {name:i for i,name in enumerate(header)}
            required = ('#chr','pos(1-based)','ref','alt')
            if not all(name in positions for name in required): raise ValueError('GRCh38 coordinate columns missing')
            available = [name for name in FEATURES if name in positions]
            if len(available)<5 or 'gnomAD4.1_joint_AF' not in available: raise ValueError('Expected score/frequency schema missing')
            for line in stream:
                if not line.strip() or line.startswith('#'): continue
                parts = line.rstrip('\r\n').split('\t')
                if len(parts)!=len(header): raise ValueError('Truncated or malformed annotation row')
                seen += 1
                chrom,pos,ref,alt = [parts[positions[name]] for name in required]
                if not pos.isdigit(): raise ValueError('Invalid dbNSFP position')
                key = (chrom.removeprefix('chr'),int(pos),ref,alt)
                if key not in targets: continue
                scores = {name:number(parts[positions[name]]) for name in available}
                prior = con.execute('SELECT scores,source_rows FROM annotations WHERE chrom=? AND pos=? AND ref=? AND alt=?',key).fetchone()
                rows = 1
                if prior:
                    rows += prior[1]
                    for name,value in json.loads(prior[0]).items():
                        if value is not None: scores[name]=max(value,scores[name]) if scores.get(name) is not None else value
                con.execute('INSERT OR REPLACE INTO annotations VALUES (?,?,?,?,?,?)',(*key,json.dumps(scores,allow_nan=False),rows))
                matched += 1
                if matched%10000==0:
                    con.commit(); print(f'Imported {matched:,} matching source rows',flush=True)
        coverage = Counter()
        for (scores,) in con.execute('SELECT scores FROM annotations'): coverage.update(k for k,v in json.loads(scores).items() if v is not None)
        result = {'source':source.name,'version':'5.4a','assembly':'GRCh38','md5':verified_md5,
            'rows_scanned':seen,'matched_source_rows':matched,'matched_variants':con.execute('SELECT COUNT(*) FROM annotations').fetchone()[0],
            'target_variants':len(targets),'feature_nonmissing_counts':dict(coverage),'features':available,'seconds':round(time.monotonic()-started,2),
            'aggregation':'Maximum finite value over transcript entries and duplicate genomic-key rows; variant-level exploratory aggregation, not transcript-specific inference.',
            'missing_policy':'Null remains missing; absent allele frequency is never assumed zero.',
            'evaluation_caution':'Upstream predictors may have trained on overlapping ClinVar records. Audit provenance before claiming independent accuracy.'}
        if not result['matched_variants']: raise ValueError('No exact allele matches; annotations not marked ready')
        put_meta(con,'annotation',result); put_meta(con,'annotation_status','ready'); con.commit()
        return result
    finally: con.close()

def compare(database, output, max_per_split=50000):
    import numpy as np
    from sklearn.ensemble import HistGradientBoostingClassifier
    from sklearn.impute import SimpleImputer
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import accuracy_score,balanced_accuracy_score,classification_report
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler
    if max_per_split<5: raise ValueError('max-per-split must be at least five')
    con = connect(database)
    try:
        if meta(con,'annotation_status')!='ready': raise ValueError('Complete and verify the dbNSFP merge first; no model comparison run')
        info=meta(con,'annotation'); splits={}
        for split in ('train','validation'):
            records=[]
            for label in CLASSES:
                records.extend(con.execute('''SELECT v.label,a.scores FROM variants v JOIN annotations a
                  USING(chrom,pos,ref,alt) WHERE v.split=? AND v.label=?
                  ORDER BY v.chrom,v.pos,v.ref,v.alt LIMIT ?''',(split,label,max_per_split//5)).fetchall())
            if set(label for label,_ in records)!=set(CLASSES): raise ValueError(f'All five labels need matched records in {split}')
            x=np.array([[json.loads(scores).get(name,np.nan) for name in FEATURES] for _,scores in records],dtype=np.float32)
            splits[split]=x,np.array([label for label,_ in records])
    finally: con.close()
    models={'linear':make_pipeline(SimpleImputer(add_indicator=True,keep_empty_features=True),StandardScaler(),LogisticRegression(max_iter=500,random_state=42)),
        'histogram':HistGradientBoostingClassifier(max_iter=150,max_leaf_nodes=15,min_samples_leaf=30,l2_regularization=5,random_state=42)}
    rows=[]
    for name,model in models.items():
        model.fit(*splits['train']); predicted=model.predict(splits['validation'][0]); y=splits['validation'][1]
        report=classification_report(y,predicted,labels=CLASSES,output_dict=True,zero_division=0)
        rows.append({'candidate':name,'accuracy':float(accuracy_score(y,predicted)),'balanced_accuracy':float(balanced_accuracy_score(y,predicted)),'macro_f1':report['macro avg']['f1-score'],'per_class':report})
    result={'status':'validation_only','test_rows_read':0,'calibration_rows_read':0,'annotation_md5':info['md5'],'features':list(FEATURES),
        'sampling':'Per-class capped pilot ordered by coordinates; not population-representative accuracy. Further training needs separate calibration and upstream overlap audit.',
        'split_label_counts':{s:dict(Counter(y.tolist())) for s,(_,y) in splits.items()},'candidates':rows,
        'selected_for_further_study':max(rows,key=lambda r:r['macro_f1'])['candidate'],'production_model_changed':False}
    output=Path(output); output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,indent=2),encoding='utf-8')
    return result

def main():
    parser=argparse.ArgumentParser(description=__doc__); sub=parser.add_subparsers(dest='command',required=True)
    p=sub.add_parser('prepare'); p.add_argument('--source',type=Path,default=ROOT/'data/raw/clinvar/variant_summary_2026-09.txt.gz'); p.add_argument('--database',type=Path,default=DEFAULT_DB); p.add_argument('--sample-modulus',type=int,default=10)
    p=sub.add_parser('merge'); p.add_argument('--source',type=Path,required=True); p.add_argument('--checksum',type=Path,required=True); p.add_argument('--database',type=Path,default=DEFAULT_DB)
    p=sub.add_parser('compare'); p.add_argument('--database',type=Path,default=DEFAULT_DB); p.add_argument('--output',type=Path,default=ROOT/'reports/annotation_v2/validation_comparison.json'); p.add_argument('--max-per-split',type=int,default=50000)
    args=parser.parse_args()
    if args.command=='prepare': result=prepare(args.source,args.database,args.sample_modulus)
    elif args.command=='merge': result=merge(args.source,args.checksum,args.database)
    else: result=compare(args.database,args.output,args.max_per_split)
    print(json.dumps(result,indent=2))

if __name__=='__main__': main()
