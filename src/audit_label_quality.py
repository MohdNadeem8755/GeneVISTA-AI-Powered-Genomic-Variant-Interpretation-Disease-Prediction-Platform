"""Read-only development-label audit against the latest locally cached snapshot."""
import argparse
from collections import Counter
import csv
import gzip
import json
from pathlib import Path
import sqlite3
import joblib
from src.annotation_pipeline import snv_key
from src.audit_annotation_labels import quality
from src.balanced_annotation import enhanced_matrix, decision_indices
from src.robust_variant_training import metrics, save_json, file_hash

ROOT = Path(__file__).resolve().parents[1]


def run(output):
    out = Path(output)
    out.mkdir(parents=True, exist_ok=False)
    db = ROOT / 'data/processed/annotation_v3/full_cohort.sqlite'
    source = ROOT / 'data/raw/clinvar/variant_summary_2026-10-05.txt.gz'
    con = sqlite3.connect(db.as_uri() + '?mode=ro', uri=True)
    report = {'scope': 'training and validation only; no held-out labels read',
              'as_of': '2026-10-29', 'latest_local_snapshot': '2026-10-05',
              'test_labels_read': 0, 'labels_modified': 0, 'database': str(db), 'source': str(source)}
    report['cohort_provenance'] = json.loads(con.execute("SELECT value FROM metadata WHERE key='cohort'").fetchone()[0])
    report['support'] = [dict(zip(['split','label','review_status','matched','count'], row)) for row in con.execute(
        "SELECT v.split,v.label,v.review_status,a.chrom IS NOT NULL,COUNT(*) FROM variants v LEFT JOIN annotations a USING(chrom,pos,ref,alt) WHERE v.split IN ('train','validation') GROUP BY v.split,v.label,v.review_status,a.chrom IS NOT NULL")]
    for row in report['support']:
        row['quality_tier'] = quality(row['review_status'])
    report['duplicate_allele_groups_in_development'] = con.execute("SELECT COUNT(*) FROM (SELECT chrom,pos,ref,alt FROM variants WHERE split IN ('train','validation') GROUP BY chrom,pos,ref,alt HAVING COUNT(*)>1)").fetchone()[0]
    report['duplicate_variant_ids_in_development'] = con.execute("SELECT COUNT(*) FROM (SELECT variation_id FROM variants WHERE split IN ('train','validation') GROUP BY variation_id HAVING COUNT(*)>1)").fetchone()[0]
    report['conflict_flagged_development_rows'] = con.execute("SELECT COUNT(*) FROM variants WHERE split IN ('train','validation') AND conflict=1").fetchone()[0]
    report['excluded_conflict_flagged_rows'] = con.execute("SELECT COUNT(*) FROM variants WHERE split='excluded' AND conflict=1").fetchone()[0]
    report['duplicate_caution'] = 'Cohort primary key already collapses alleles. Conflict flag combines differing labels and differing IDs; it is not exclusively clinical disagreement.'
    ids = {row[0] for row in con.execute("SELECT variation_id FROM variants WHERE split IN ('train','validation')")}
    stage = sqlite3.connect(out / 'latest_development.sqlite')
    stage.execute('CREATE TABLE latest(id TEXT,chrom TEXT,pos INTEGER,ref TEXT,alt TEXT,label TEXT,review TEXT,evaluated TEXT)')
    batch, scanned = [], 0
    print('Scanning latest local ClinVar snapshot for development IDs', flush=True)
    with gzip.open(source, 'rt', encoding='utf-8-sig', newline='') as stream:
        for row in csv.DictReader(stream, delimiter='\t'):
            scanned += 1
            if row['VariationID'] not in ids or row.get('OriginSimple') not in ('germline','unknown'):
                continue
            key = snv_key(row)
            if key is None:
                continue
            batch.append((row['VariationID'], *key, row['ClinicalSignificance'], row['ReviewStatus'], row['LastEvaluated']))
            if len(batch) >= 20000:
                stage.executemany('INSERT INTO latest VALUES (?,?,?,?,?,?,?,?)', batch)
                stage.commit(); batch.clear()
    stage.executemany('INSERT INTO latest VALUES (?,?,?,?,?,?,?,?)', batch)
    stage.execute('CREATE INDEX latest_key ON latest(id,chrom,pos,ref,alt)')
    stage.commit(); stage.close()
    del ids
    con.execute('ATTACH DATABASE ? AS latestdb', ((out / 'latest_development.sqlite').resolve().as_uri() + '?mode=ro',))
    con.execute('CREATE TEMP VIEW observations AS SELECT id,chrom,pos,ref,alt,COUNT(*) n,COUNT(DISTINCT label) labels,MIN(label) label,MIN(review) review,MAX(evaluated) evaluated FROM latestdb.latest GROUP BY id,chrom,pos,ref,alt')
    stats, transitions = Counter(), Counter()
    from datetime import datetime, date
    print('Comparing exact alleles and label recency', flush=True)
    with (out / 'review_queue.csv').open('w', newline='', encoding='utf-8') as stream:
        writer = csv.writer(stream)
        writer.writerow(['variation_id','split','chrom','pos','ref','alt','old_label','latest_label','latest_review','last_evaluated','reason'])
        rows = con.execute("SELECT v.variation_id,v.split,v.chrom,v.pos,v.ref,v.alt,v.label,o.label,o.review,o.evaluated,o.n,o.labels FROM variants v LEFT JOIN observations o ON v.variation_id=o.id AND v.chrom=o.chrom AND v.pos=o.pos AND v.ref=o.ref AND v.alt=o.alt WHERE v.split IN ('train','validation')")
        for row in rows:
            vid, split, chrom, pos, ref, alt, old, new, review, evaluated, n, labels = row
            stats['development_rows'] += 1
            reasons = []
            if n is None:
                stats['missing_exact_latest_match'] += 1
                reasons.append('missing_exact_latest_match_not_proof_of_withdrawal')
            else:
                stats['matched_rows'] += 1
                stats['latest_source_duplicate_rows_beyond_first'] += n - 1
                if labels > 1:
                    stats['ambiguous_latest_labels'] += 1
                    reasons.append('ambiguous_latest_labels')
                elif new != old:
                    stats['changed_label'] += 1
                    transitions[(old, new)] += 1
                    reasons.append('changed_label')
                if 'conflict' in new.lower() or 'conflict' in review.lower() and 'no conflicts' not in review.lower():
                    stats['latest_conflicting_classification'] += 1
                    reasons.append('latest_conflicting_classification')
                if n == 1:
                    try:
                        dt = datetime.strptime(evaluated, '%b %d, %Y').date()
                        age = (date(2026,10,29) - dt).days
                        if age < 0:
                            stats['future_evaluation_date'] += 1
                        elif age > 5 * 365.25:
                            stats['last_evaluated_over_five_years'] += 1
                            reasons.append('old_evaluation_date_not_proof_label_wrong')
                        else:
                            stats['last_evaluated_within_five_years'] += 1
                    except (ValueError, TypeError):
                        stats['missing_or_unparseable_evaluation_date'] += 1
                else:
                    stats['recency_not_assessed_duplicate_observations'] += 1
            if reasons:
                writer.writerow([vid,split,chrom,pos,ref,alt,old,new,review,evaluated,';'.join(reasons)])
    report['latest_comparison'] = dict(stats)
    report['label_transitions'] = [{'old':a,'latest':b,'count':n} for (a,b),n in transitions.most_common()]
    print('Evaluating expert-reviewed validation subset', flush=True)
    expert_rows = con.execute("SELECT v.label,a.scores,v.gene,v.name,v.variant_type FROM variants v JOIN annotations a USING(chrom,pos,ref,alt) WHERE v.split='validation' AND v.review_status IN ('reviewed by expert panel','practice guideline')").fetchall()
    bundle = joblib.load(ROOT / 'models/annotation_prediction_v4/predictor.joblib')
    if expert_rows:
        x = enhanced_matrix(expert_rows, bundle['features'], bundle['gene_vocab'])
        model = bundle['classifier']
        p = model.predict_proba(x)
        pred = model.classes_[decision_indices(p, bundle['class_priors'], bundle['decision_exponent'])]
        report['expert_validation'] = metrics([r[0] for r in expert_rows], pred, model.classes_)
        total = sum(r['count'] for r in report['support'] if r['split']=='validation' and r['matched'])
        report['expert_validation']['coverage_of_matched_validation'] = len(expert_rows) / total
        report['expert_validation']['caution'] = 'Previously reused validation; source-snapshot expert labels, not an independent external test. Balanced accuracy averages recall over supported classes.'
    report['source_sha256'] = file_hash(source)
    report['source_rows_scanned'] = scanned
    save_json(out / 'audit.json', report)
    con.close()
    print(json.dumps({'comparison':report['latest_comparison'], 'expert_count':len(expert_rows)}, indent=2), flush=True)
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True)
    run(parser.parse_args().output)
