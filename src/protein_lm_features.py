"""Import cached protein language-model features without clinical-label columns."""
import argparse
from collections import Counter
from contextlib import closing
import hashlib
import json
import math
from pathlib import Path
import re
import sqlite3
import statistics
import time
try:
    from isal import igzip
except ImportError:
    import gzip as igzip  # Small tests and portable fallback; install isal for large files.
from src.annotation_pipeline import meta

ROOT = Path(__file__).resolve().parents[1]
POOLED = ('ESM1b_score', 'popEVE_score', 'MutFormer_score')
ALIGNED = ('MPC_score', 'VEST4_score')
FEATURES = tuple(f'{name}_{stat}' for name in POOLED for stat in ('min','median','max')) + tuple(f'{name}_mane' for name in ALIGNED)
REQUIRED = ('#chr','pos(1-based)','ref','alt','Ensembl_transcriptid','MANE') + POOLED + ALIGNED


def finite_values(text):
    values=[]
    for token in re.split('[;,]',text):
        if token in ('','.','NA','nan'):continue
        try:value=float(token)
        except ValueError:raise ValueError('Malformed numeric protein score')
        if not math.isfinite(value):raise ValueError('Non-finite protein score')
        values.append(value)
    return values


def parse_features(row):
    scores={}
    # Pooling is explicit: dbNSFP does not document transcript alignment for all LMs.
    for name in POOLED:
        values=finite_values(row[name])
        for stat,func in [('min',min),('median',statistics.median),('max',max)]:
            scores[f'{name}_{stat}']=func(values) if values else None
    transcripts=row['Ensembl_transcriptid'].split(';')
    mane=row['MANE'].split(';')
    for name in ALIGNED:
        values=row[name].split(';')
        chosen=[]
        if len(mane)==len(transcripts)==len(values):
            for flag,value in zip(mane,values):
                if flag=='Select':chosen.extend(finite_values(value))
        scores[f'{name}_mane']=max(chosen) if chosen else None
    return scores


def source_signature(path):
    stat=path.stat()
    return {'bytes':stat.st_size,'mtime_ns':stat.st_mtime_ns}


def extract(source,checksum,cohort,output):
    source,output=Path(source),Path(output)
    staging=output.with_suffix('.building.sqlite')
    if output.exists() or staging.exists():raise ValueError('Choose a new output; existing data is preserved')
    output.parent.mkdir(parents=True,exist_ok=True)
    started=time.monotonic();signature=source_signature(source)
    with closing(sqlite3.connect(Path(cohort).resolve().as_uri()+'?mode=ro',uri=True)) as con:
        if meta(con,'annotation_status')!='ready':raise ValueError('Source cohort is not ready')
        provenance=meta(con,'annotation')
        # IDs and labels are not loaded by the feature importer.
        targets=set(con.execute("SELECT v.chrom,v.pos,v.ref,v.alt FROM variants v JOIN annotations a USING(chrom,pos,ref,alt) WHERE v.split IN ('train','validation')"))
    expected=re.search(r'\b[0-9a-fA-F]{32}\b',Path(checksum).read_text(encoding='utf-8-sig'))
    if expected is None or expected[0].lower()!=provenance['md5']:raise ValueError('Checksum provenance mismatch')
    print('Verifying compressed-source MD5',flush=True)
    with source.open('rb') as f:digest=hashlib.file_digest(f,'md5').hexdigest()
    if digest!=expected[0].lower():raise ValueError('Corrupt annotation download')
    counts=Counter();coverage=Counter();pending=[]
    report={'status':'extracting','source':source.name,'source_signature':signature,'md5':digest,
            'target_rows':len(targets),'features':list(FEATURES),'test_rows_read':0,
            'aggregation':'LM min/median/max across available isoforms and genomic-key rows; not transcript-specific. MPC/VEST4 only with exact MANE Select vector alignment.',
            'excluded_inputs':['clinical classifications','review status','ClinVar identifiers','disease labels'],
            'source_columns':'https://dist.genos.us/release/dbNSFP5.4a_variant.columns.txt'}
    con=sqlite3.connect(staging)
    try:
        con.executescript('CREATE TABLE features(chrom TEXT,pos INT,ref TEXT,alt TEXT,scores TEXT,PRIMARY KEY(chrom,pos,ref,alt)); CREATE TABLE metadata(key TEXT PRIMARY KEY,value TEXT); CREATE TABLE observations(chrom TEXT,pos INT,ref TEXT,alt TEXT,scores TEXT);')
        with igzip.open(source,'rb') as f:
            header=f.readline().decode().rstrip('\r\n').split('\t')
            indices={name:header.index(name) for name in REQUIRED}
            if tuple(header[:4])!=REQUIRED[:4]:raise ValueError('Unexpected coordinate columns')
            for line in f:
                prefix=line.split(b'\t',4)
                if len(prefix)<5:raise ValueError('Truncated annotation row')
                key=(prefix[0].decode(),int(prefix[1]),prefix[2].decode(),prefix[3].decode())
                counts['source_rows']+=1
                if key in targets:
                    parts=line.decode().rstrip('\r\n').split('\t')
                    if len(parts)!=len(header):raise ValueError('Invalid annotation row width')
                    row={name:parts[i] for name,i in indices.items()}
                    scores=parse_features(row)
                    # Genome rows can repeat for different protein changes; collect each row.
                    pending.append((*key,json.dumps(scores,allow_nan=False)))
                    counts['matched_source_rows']+=1
                if len(pending)>=5000:
                    # Raw observations are combined deterministically after import.
                    con.executemany('INSERT INTO observations VALUES (?,?,?,?,?)',pending);pending=[];con.commit()
                if counts['source_rows']%2000000==0:
                    report.update(counts);report['seconds']=time.monotonic()-started
                    output.with_suffix('.progress.json').write_text(json.dumps(report),encoding='utf-8')
                    print(json.dumps({'rows':counts['source_rows'],'matched':counts['matched_source_rows'],'seconds':round(report['seconds'])}),flush=True)
        if pending:con.executemany('INSERT INTO observations VALUES (?,?,?,?,?)',pending)
        if source_signature(source)!=signature:raise ValueError('Source changed during extraction')
        con.commit()
        con.execute('CREATE INDEX observations_key ON observations(chrom,pos,ref,alt)')
        current=None;group=[];batch=[]
        def aggregate(key,rows):
            values={}
            for name in FEATURES:
                observed=[r[name] for r in rows if r[name] is not None]
                # A median of row medians is explicitly a pooled summary, not an isoform score.
                func=min if name.endswith('_min') else statistics.median if name.endswith('_median') else max
                values[name]=func(observed) if observed else None
            coverage.update(name for name,value in values.items() if value is not None)
            return (*key,json.dumps(values,allow_nan=False))
        for chrom,pos,ref,alt,scores in con.execute('SELECT * FROM observations ORDER BY chrom,pos,ref,alt'):
            key=(chrom,pos,ref,alt)
            if current is not None and current!=key:
                batch.append(aggregate(current,group));group=[]
                if len(batch)>=5000:con.executemany('INSERT INTO features VALUES (?,?,?,?,?)',batch);batch=[]
            current=key;group.append(json.loads(scores))
        if current is not None:batch.append(aggregate(current,group))
        con.executemany('INSERT INTO features VALUES (?,?,?,?,?)',batch)
        matched=con.execute('SELECT COUNT(*) FROM features').fetchone()[0]
        if matched!=len(targets):raise ValueError(f'Unexpected coverage: {matched}/{len(targets)} cached alleles')
        con.execute('DROP TABLE observations')
        report.update(counts);report.update(status='ready',matched_variants=matched,coverage=dict(coverage),seconds=time.monotonic()-started)
        con.execute('INSERT INTO metadata VALUES (?,?)',('protein_features',json.dumps(report,allow_nan=False)));con.commit()
    finally:con.close()
    staging.replace(output)
    output.with_suffix('.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    output.with_suffix('.progress.json').write_text(json.dumps(report),encoding='utf-8')
    print(json.dumps({'status':'ready','matched':matched,'coverage':dict(coverage)}),flush=True)
    return report


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,default=ROOT/'data/incoming/dbNSFP5.4a_grch38.gz')
    p.add_argument('--checksum',type=Path,default=ROOT/'data/incoming/dbNSFP5.4a_grch38.gz.md5')
    p.add_argument('--cohort',type=Path,default=ROOT/'data/processed/annotation_v3/full_cohort.sqlite')
    p.add_argument('--output',type=Path,default=ROOT/'data/processed/protein_lm_v1/features.sqlite')
    a=p.parse_args();extract(a.source,a.checksum,a.cohort,a.output)
