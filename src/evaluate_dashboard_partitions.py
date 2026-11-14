"""Score the frozen fresh ensemble; never fit or select on evaluation results."""
import gzip
import json
import sqlite3
from datetime import date, datetime, timezone
from pathlib import Path
import joblib
import numpy as np
import pandas as pd
from src.annotation_pipeline import bucket
from src.biological_features import clingen_index, gene_context, enrich
from src.retune_biological_models import stopping_mask
from src.robust_variant_training import file_hash, metrics, save_json
from src.train_annotation_candidates import INDEPENDENT, SEQUENCE, matrix
from src.tune_biological_models import predict_bundle

ROOT=Path(__file__).resolve().parents[1]
RUN=ROOT/'reports/biological_retuning_fresh_20261031'
OUTPUT=RUN/'partition_evaluation.json'


def run(directory=RUN):
    RUN=Path(directory).resolve()
    OUTPUT=RUN/'partition_evaluation.json'
    if OUTPUT.exists():
        raise ValueError('Evaluation already exists; preserve the saved result')
    selection=json.loads((RUN/'selection.json').read_text())
    model_path=RUN/'selected.joblib'
    digest=file_hash(model_path)
    assert digest == selection['model_sha256'], 'Frozen model changed'
    feature_path=ROOT/'data/processed/biological_features_20261029/report.json'
    protocol=json.loads((RUN/'protocol.json').read_text())
    assert file_hash(feature_path)==protocol['feature_report_sha256']
    meta=json.loads(feature_path.read_text())
    for source in meta['sources'].values():
        assert file_hash(Path(source['path']))==source['sha256'], 'Feature source changed'
    brca=json.load(gzip.open(meta['sources']['brca1']['path'],'rt'))
    tp53=json.load(gzip.open(meta['sources']['tp53']['path'],'rt'))
    index=clingen_index(meta['sources']['clingen']['path'],date.fromisoformat(meta['cutoff']))
    contexts={g:gene_context(g,index) for g in index}
    bundle=joblib.load(model_path)
    revised=bundle.get('metadata',{}).get('feature_profile')=='variant_identity_v2'
    base_names=tuple(INDEPENDENT)+SEQUENCE
    database=ROOT/'data/processed/annotation_v3/full_cohort.sqlite'
    cache=ROOT/'data/processed/biological_features_20261029/features.sqlite'
    con=sqlite3.connect(database.as_uri()+'?mode=ro',uri=True)
    con.execute('ATTACH DATABASE ? AS bio',(cache.as_uri()+'?mode=ro',))
    provenance=json.loads(con.execute("SELECT value FROM metadata WHERE key='annotation'").fetchone()[0])
    assert provenance['md5']==meta['annotation_provenance']['md5']
    overlap=con.execute("SELECT COUNT(*) FROM (SELECT chrom,pos FROM variants WHERE split IN ('train','validation','test') GROUP BY chrom,pos HAVING COUNT(DISTINCT split)>1)").fetchone()[0]
    assert overlap==0, 'Partition locus overlap'
    con.create_function('sample_locus',2,lambda c,p:bucket(f'biology-v1|{c}|{p}',5))
    con.create_function('legacy_sample',2,lambda c,p:bucket(f'genevista-annotation-v2-sample|{c}|{p}',10))
    save_json(RUN/'partition_evaluation_protocol.json',dict(model_sha256=digest,
        train='Original sampled fitting rows only; stopping rows reported separately',
        validation='Original sampled selection rows; reused development evaluation',
        internal_test='Existing annotation-matched test extension excluding legacy v2 sample, as for active v4; historical test previously inspected in this project, not untouched external validation',
        model_refit=False,decision_rule='Frozen saved-weight soft voting argmax',feature_cutoff=meta['cutoff']))
    actual={k:[] for k in ('training','early_stopping','validation','internal_test')}
    predicted={k:[] for k in actual}
    for split in ('train','validation','test'):
        condition='legacy_sample(v.chrom,v.pos)<>0' if split=='test' else 'sample_locus(v.chrom,v.pos)=0'
        cursor=con.execute(f'''SELECT v.label,a.scores,v.gene,v.name,v.variant_type,
            b.scores,v.chrom,v.pos,v.ref,v.alt FROM variants v JOIN annotations a USING(chrom,pos,ref,alt)
            LEFT JOIN bio.features b USING(chrom,pos,ref,alt)
            WHERE v.split=? AND {condition}''',(split,))
        count=0
        while rows:=cursor.fetchmany(10000):
            base=pd.DataFrame(matrix([r[:5] for r in rows],base_names),columns=base_names)
            if revised:
                from src.variant_feature_revision import revised_base_matrix
                base=revised_base_matrix([(*r[:5],r[8],r[9]) for r in rows])
            if split!='test' and any(r[5] is None for r in rows):
                raise ValueError('Missing cached development features')
            additions=[json.loads(r[5]) if r[5] is not None else enrich(r[2],r[3],json.loads(r[1]),contexts,brca,tp53) for r in rows]
            extra=pd.DataFrame(additions).reindex(columns=meta['numeric']+meta['categorical'])
            frame=pd.concat([base.drop(columns=[c for c in extra if c in base]),extra],axis=1)
            p=predict_bundle(bundle,frame)
            assert np.isfinite(p).all()
            np.testing.assert_allclose(p.sum(axis=1),1,atol=1e-6)
            pred=bundle['classes'][p.argmax(axis=1)]
            y=np.asarray([r[0] for r in rows])
            if split=='train':
                stop=stopping_mask([f'{r[6]}:{r[7]}' for r in rows])
                masks={'training':~stop,'early_stopping':stop}
            else:
                masks={'internal_test' if split=='test' else 'validation':np.ones(len(rows),dtype=bool)}
            for key,mask in masks.items():
                actual[key].extend(y[mask]); predicted[key].extend(pred[mask])
            count+=len(rows)
            print(split,count,flush=True)
    con.close()
    scores={key:metrics(actual[key],predicted[key],bundle['classes']) for key in actual}
    audit=json.loads((RUN/'data_audit.json').read_text())
    assert scores['training']['n']==audit['fit_rows']
    assert scores['early_stopping']['n']==audit['stopping_rows']
    assert scores['validation']['n']==audit['validation_rows']
    for key in ('accuracy','balanced_accuracy','macro_f1'):
        assert abs(scores['validation'][key]-selection['selected'][key])<1e-10, 'Validation reproduction failed'
    assert scores['internal_test']['n']==200134, 'Historical test extension differs'
    assert file_hash(model_path)==digest
    result=dict(status='complete',model_sha256=digest,evaluated_at=datetime.now(timezone.utc).isoformat(),
        scores=scores,model_refit=False,production_model_changed=False,locus_overlap=overlap,
        test_scope='Previously inspected internal test extension; not untouched later/external validation')
    temporary=OUTPUT.with_suffix('.tmp')
    save_json(temporary,result)
    temporary.replace(OUTPUT)
    print({k:{m:v[m] for m in ('n','accuracy','balanced_accuracy','macro_f1')} for k,v in scores.items()},flush=True)


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument('--run',default=str(RUN))
    run(parser.parse_args().run)
