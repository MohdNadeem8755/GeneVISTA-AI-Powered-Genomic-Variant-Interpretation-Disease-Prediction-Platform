"""Audit training/validation label support without reading held-out partitions."""
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path
from src.annotation_pipeline import connect, meta, ROOT
from src.variant_features import CLASSES

def quality(review):
    if review == 'practice guideline':
        return 'expert'
    if review == 'reviewed by expert panel':
        return 'expert'
    if review == 'criteria provided, multiple submitters, no conflicts':
        return 'multiple_no_conflicts'
    if review == 'criteria provided, single submitter':
        return 'single_submitter'
    return 'other'

def audit(database, output):
    output=Path(output); output.mkdir(parents=True, exist_ok=True)
    groups=Counter(); support=defaultdict(Counter)
    with connect(database) as con:
        if meta(con,'annotation_status')!='ready':
            raise ValueError('Verified annotations required')
        rows=con.execute("""SELECT v.split,v.label,v.review_status,v.gene,
          CASE WHEN a.chrom IS NULL THEN 0 ELSE 1 END
          FROM variants v LEFT JOIN annotations a USING(chrom,pos,ref,alt)
          WHERE v.split IN ('train','validation')""")
        for split,label,review,gene,matched in rows:
            tier=quality(review)
            groups[(split,tier,label,matched)]+=1
            if not matched or not gene or gene=='-' or any(c in gene for c in ';,|'):
                continue
            for cohort in ('all', 'high_quality') if tier in ('expert','multiple_no_conflicts') else ('all',):
                support[(cohort,gene)][(split,label)]+=1
    gene_rows=[]
    for (cohort,gene),counts in sorted(support.items()):
        train=[counts[('train',label)] for label in CLASSES]
        validation=[counts[('validation',label)] for label in CLASSES]
        gene_rows.append({'cohort':cohort,'gene':gene,'train_total':sum(train),
            'validation_total':sum(validation),'min_train_class':min(train),
            'min_validation_class':min(validation),'eligible':min(train)>=30 and min(validation)>=10,
            **{f'train_{label}':count for label,count in zip(CLASSES,train)},
            **{f'validation_{label}':count for label,count in zip(CLASSES,validation)}})
    with (output/'gene_support.csv').open('w',newline='',encoding='utf-8') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(gene_rows[0]) if gene_rows else ['gene'])
        writer.writeheader(); writer.writerows(gene_rows)
    report={'partitions_read':['train','validation'],'test_rows_read':0,'calibration_rows_read':0,
        'eligibility':'At least 30 training and 10 validation examples in EACH of five classes; feasibility screen, not proof of reliability.',
        'quality_definition':'Expert panel/practice guideline, or multiple submitters with no conflicts; exact source labels preserved.',
        'counts':[{'split':s,'quality':q,'label':l,'annotation_matched':bool(m),'count':n}
                  for (s,q,l,m),n in sorted(groups.items())],
        'eligible_genes':[r for r in gene_rows if r['eligible']],
        'top_support':sorted(gene_rows,key=lambda r:(r['min_train_class'],r['min_validation_class']),reverse=True)[:20]}
    (output/'label_audit.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    return report

if __name__=='__main__':
    report=audit(ROOT/'data/processed/annotation_v2/download_verified_cohort.sqlite',ROOT/'reports/annotation_v3')
    print(json.dumps({'eligible_genes':report['eligible_genes'],'top_support':report['top_support'][:3]},indent=2))
