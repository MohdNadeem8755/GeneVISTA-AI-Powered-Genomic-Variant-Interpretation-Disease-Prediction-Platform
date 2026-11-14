"""Deterministic biological feature enrichment; never reads clinical targets.

Builds a separate cache for train/validation only and audits score availability
across evaluation partitions without reading their labels. No model is fitted.
"""
import argparse
from collections import Counter, defaultdict
from contextlib import closing
import csv
from datetime import date
import gzip
import json
import math
from pathlib import Path
import sqlite3

from src.brca1_mavedb import assess_variant as brca_assess
from src.tp53_evidence import assess_variant as tp53_assess
from src.annotation_pipeline import FEATURES
from src.robust_variant_training import file_hash, save_json

ROOT = Path(__file__).resolve().parents[1]
POPULATION = tuple(n for n in FEATURES if n.startswith('gnomAD'))
PREDICTORS = tuple(n for n in FEATURES if n not in POPULATION and not n.startswith(('GERP','phyloP')))


def finite(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


def clingen_index(path, cutoff):
    """Gene-context only. Preserve disease relations separately; no patient inference."""
    index = defaultdict(list)
    with gzip.open(path, 'rt', encoding='utf-8-sig', newline='') as stream:
        for row in csv.DictReader(stream):
            curated = row.get('classification_date_iso', '')
            released = row.get('source_file_created', '')
            if not curated or not released:
                continue
            if date.fromisoformat(curated) <= cutoff and date.fromisoformat(released) <= cutoff:
                index[row['gene_symbol']].append(row)
    return dict(index)


def gene_context(gene, index):
    # Multi-gene strings are not silently assigned to one gene.
    rows = index.get(gene, [])
    if not rows:
        return {'clingen_context_available': 0}
    relations = {(r['mondo_id'], r['mode_of_inheritance'], r['classification']) for r in rows}
    result = {'clingen_context_available': 1,
              'clingen_disease_count': len({r[0] for r in relations})}
    for level in ('Definitive','Strong','Moderate','Limited','Disputed','Refuted'):
        result['clingen_' + level.lower() + '_count'] = sum(r[2] == level for r in relations)
    for moi in ('AD','AR','XL','XLD','XLR','MT'):
        result['inheritance_' + moi] = int(any(r[1] == moi for r in relations))
    result['inheritance_other_count'] = sum(r[1] not in ('AD','AR','XL','XLD','XLR','MT') for r in relations)
    return result


def enrich(gene, name, scores, contexts, brca, tp53):
    result = dict(contexts.get(gene, {'clingen_context_available': 0}))
    for field in POPULATION:
        value = scores.get(field)
        valid = finite(value) and value >= 0 and (field.endswith('_AN') or value <= 1)
        result[field] = float(value) if valid else None
        result[field + '__missing'] = int(not valid)
        if not field.endswith('_AN'):
            # Missing and observed zero have different meanings.
            result[field + '__neglog10'] = -math.log10(max(value, 1e-12)) if valid else None
    result.update(brca1_function_score=None, brca1_rna_score=None,
                  tp53_spliceai=None, assay_match_available=0)
    record = {'gene_symbol': gene, 'name': name}
    if gene == 'BRCA1':
        evidence = brca_assess(record, brca)
        if evidence['status'] == 'available':
            result['brca1_function_score'] = evidence['score']
            result['brca1_rna_score'] = evidence.get('rna_score')
            result['assay_match_available'] = 1
    elif gene == 'TP53':
        evidence = tp53_assess(record, tp53)
        if evidence['status'] == 'available':
            value = evidence['computational'].get('spliceai')
            result['tp53_spliceai'] = value if finite(value) and 0 <= value <= 1 else None
            functional = evidence['functional']
            if functional['status'] == 'matched':
                for assay, value in functional['assays'].items():
                    result['tp53_assay_' + assay] = None if value in ('NA','No data','') else value
                result['assay_match_available'] = 1
    return result


def build(output, cutoff):
    out = Path(output)
    out.mkdir(parents=True, exist_ok=False)
    sources = {'brca1': ROOT/'data/processed/brca1_mavedb_v1.json.gz',
               'tp53': ROOT/'data/processed/tp53_evidence_v1.json.gz',
               'clingen': ROOT/'data/processed/clingen_gene_evidence_mapped.csv.gz'}
    brca = json.load(gzip.open(sources['brca1'], 'rt'))
    tp53 = json.load(gzip.open(sources['tp53'], 'rt'))
    index = clingen_index(sources['clingen'], cutoff)
    contexts = {gene: gene_context(gene, index) for gene in index}
    database = ROOT/'data/processed/annotation_v3/full_cohort.sqlite'
    con = sqlite3.connect(database.as_uri() + '?mode=ro', uri=True)
    target = sqlite3.connect(out/'features.sqlite')
    target.executescript('CREATE TABLE features(chrom TEXT,pos INTEGER,ref TEXT,alt TEXT,split TEXT,scores TEXT,PRIMARY KEY(chrom,pos,ref,alt)); CREATE TABLE metadata(key TEXT PRIMARY KEY,value TEXT); CREATE TABLE gene_disease(gene TEXT,mondo TEXT,inheritance TEXT,validity TEXT,curated TEXT,source_url TEXT);')
    target.executemany('INSERT INTO gene_disease VALUES (?,?,?,?,?,?)',
        [(g,r['mondo_id'],r['mode_of_inheritance'],r['classification'],r['classification_date_iso'],r['report_url']) for g,rows in index.items() for r in rows])
    counts = Counter(); coverage = defaultdict(Counter); predictor_coverage = defaultdict(Counter)
    feature_names = set(); categorical = set(); batch = []
    # Label is deliberately absent from the SELECT and from the feature cache.
    rows = con.execute("SELECT v.chrom,v.pos,v.ref,v.alt,v.split,v.gene,v.name,a.scores FROM variants v JOIN annotations a USING(chrom,pos,ref,alt) WHERE v.split IN ('train','validation','calibration','test')")
    for chrom,pos,ref,alt,split,gene,name,encoded in rows:
        scores = json.loads(encoded)
        counts[split] += 1
        for field in PREDICTORS:
            if finite(scores.get(field)):
                predictor_coverage[split][field] += 1
        values = enrich(gene,name,scores,contexts,brca,tp53)
        for key,value in values.items():
            if value is not None:
                coverage[split][key] += int(value != 0) if key.endswith('_available') else 1
        if split in ('train','validation'):
            feature_names.update(values)
            categorical.update(k for k,v in values.items() if isinstance(v,str))
            batch.append((chrom,pos,ref,alt,split,json.dumps(values,allow_nan=False)))
        if len(batch) >= 10000:
            target.executemany('INSERT INTO features VALUES (?,?,?,?,?,?)',batch)
            target.commit(); batch.clear()
        if sum(counts.values()) % 250000 == 0:
            print(f'Processed {sum(counts.values()):,} matched variants',flush=True)
    target.executemany('INSERT INTO features VALUES (?,?,?,?,?,?)',batch)
    report = {'status':'ready','cutoff':cutoff.isoformat(),'counts':dict(counts),
              'stored_partitions':['train','validation'],'labels_read':0,
              'coverage':{s:dict(c) for s,c in coverage.items()},
              'pretrained_score_coverage':{s:dict(c) for s,c in predictor_coverage.items()},
              'pretrained_training_overlap':'UNRESOLVED: exact versioned upstream training/evaluation variant manifests are not locally available. Coverage is not proof of overlap or independence.',
              'sources':{k:{'path':str(p),'sha256':file_hash(p)} for k,p in sources.items()},
              'annotation_provenance':json.loads(con.execute("SELECT value FROM metadata WHERE key='annotation'").fetchone()[0]),
              'numeric':sorted(feature_names-categorical),'categorical':sorted(categorical),
              'excluded':['clinical labels','review status','PS3/BS3/PP3/BP4 criterion codes','disease match derived from target'],
              'limitations':['Static research snapshot, not historical as-of training. Assay release availability has not been verified for earlier prediction times.',
                             'SpliceAI coverage is TP53 only. Missing scores are not zero.',
                             'Assays are independently measured, but may be evidence used in ClinVar classifications.',
                             'ClinGen features are gene-level context, not variant pathogenicity or patient inheritance.',
                             'Population frequencies are reused from the existing verified cache; no new population dataset was downloaded.']}
    save_json(out/'report.json',report)
    target.execute('INSERT INTO metadata VALUES (?,?)',('biological_features',json.dumps(report)))
    target.commit(); target.close(); con.close()
    save_json(out/'profiles.json',{'biological_context':sorted(feature_names),
                                 'primary_baseline_additions':'Population/conservation/sequence plus biological_context; exclude pretrained pathogenicity scores.',
                                 'sensitivity_only_pretrained':list(PREDICTORS),
                                 'assay_ablation':'Repeat without brca1_*, tp53_assay_* and assay_match_available to assess shared-evidence dependence.'})
    print(json.dumps({'counts':dict(counts),'feature_count':len(feature_names)},indent=2))


def load_development(cohort, cache, split, modulus=1, include_pretrained=False, return_groups=False, feature_profile='legacy', minority_all=False):
    """Training adapter: preserve existing partitions; deterministic locus sampling.

    Returns (DataFrame, labels, categorical_columns). Imputation/encoding must be
    fitted inside training folds by the caller. Test loading is intentionally denied.
    """
    if split not in ('train','validation') or modulus < 1:
        raise ValueError('Only development splits and a positive modulus are allowed')
    if feature_profile not in ('legacy','variant_identity_v2'):
        raise ValueError('Unknown feature profile')
    if feature_profile!='legacy' and include_pretrained:
        raise ValueError('Revised profile excludes unaudited pretrained scores')
    import pandas as pd
    import numpy as np
    from src.annotation_pipeline import bucket
    from src.train_annotation_candidates import matrix, INDEPENDENT, SEQUENCE
    with closing(sqlite3.connect(Path(cohort).resolve().as_uri()+'?mode=ro',uri=True)) as con:
        con.execute('ATTACH DATABASE ? AS bio',(Path(cache).resolve().as_uri()+'?mode=ro',))
        metadata=json.loads(con.execute("SELECT value FROM bio.metadata WHERE key='biological_features'").fetchone()[0])
        provenance=json.loads(con.execute("SELECT value FROM metadata WHERE key='annotation'").fetchone()[0])
        if metadata['annotation_provenance']['md5'] != provenance['md5']:
            raise ValueError('Annotation provenance mismatch')
        con.create_function('sample_locus',2,lambda chrom,pos:bucket(f'biology-v1|{chrom}|{pos}',modulus))
        con.create_function('stop_locus',2,lambda chrom,pos:bucket(f'step4-redo-stop|{chrom}:{pos}',10))
        base_names=tuple(FEATURES if include_pretrained else INDEPENDENT)+SEQUENCE
        condition='sample_locus(v.chrom,v.pos)=0'
        if minority_all:
            if split!='train':
                raise ValueError('Class sampling allowed only for training')
            condition="(sample_locus(v.chrom,v.pos)=0 OR (v.label<>'Uncertain significance' AND stop_locus(v.chrom,v.pos)<>0))"
        rows=con.execute('SELECT v.label,a.scores,v.gene,v.name,v.variant_type,b.scores,v.chrom,v.pos,v.ref,v.alt FROM variants v JOIN annotations a USING(chrom,pos,ref,alt) LEFT JOIN bio.features b USING(chrom,pos,ref,alt) WHERE v.split=? AND '+condition,(split,)).fetchall()
        if not rows or any(r[5] is None for r in rows):
            raise ValueError('Missing feature rows; build the complete matching cache first')
        base=pd.DataFrame(matrix([r[:5] for r in rows],base_names),columns=base_names)
        categorical=list(metadata['categorical'])
        if feature_profile=='variant_identity_v2':
            from src.variant_feature_revision import revised_base_matrix, CATEGORICAL
            base=revised_base_matrix([(*r[:5],r[8],r[9]) for r in rows])
            categorical.extend(CATEGORICAL)
        added=pd.DataFrame([json.loads(r[5]) for r in rows]).reindex(columns=metadata['numeric']+metadata['categorical'])
        # The enriched population fields replace base values with validated ranges.
        base=base.drop(columns=[n for n in added if n in base],errors='ignore')
        result=(pd.concat([base,added],axis=1),pd.Series([r[0] for r in rows]),categorical)
        if return_groups:
            return (*result, np.asarray([f'{r[6]}:{r[7]}' for r in rows]))
        return result


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',required=True)
    parser.add_argument('--cutoff',default='2026-10-29')
    args=parser.parse_args()
    build(args.output,date.fromisoformat(args.cutoff))
