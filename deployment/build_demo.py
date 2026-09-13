"""Build a bounded, reproducible demo from the local read-only scientific database."""
from pathlib import Path
import json
import sqlite3
from contextlib import closing

ROOT = Path(__file__).resolve().parents[1]
GENES = ['BRCA1','BRCA2','TP53','NF1','FBN1','MYH7','MYBPC3','SCN1A','SCN2A',
         'CFTR','LDLR','MLH1','MSH2','PMS2','RET','PTEN','PMP22','MFN2','TGFBR1','TGFBR2']
STATES = ['vus','pathogenic','benign','conflicting','other']

def build():
    source = ROOT/'data/app/genevista.sqlite'
    target = ROOT/'data/demo/genevista.sqlite'
    if not source.is_file(): raise FileNotFoundError(source)
    target.parent.mkdir(parents=True,exist_ok=True)
    temporary = target.with_suffix('.building.sqlite')
    if temporary.exists(): raise FileExistsError('Remove the previous demo build temporary file before retrying.')
    try:
        with closing(sqlite3.connect(temporary,uri=True)) as db:
            db.execute('ATTACH DATABASE ? AS full_data',(source.resolve().as_uri()+'?mode=ro',))
            for name in ['variants','timelines','variant_gene_symbols','metadata']:
                sql=db.execute("SELECT sql FROM full_data.sqlite_master WHERE type='table' AND name=?",(name,)).fetchone()[0]
                db.execute(sql)
            db.execute('CREATE TEMP TABLE selected (variation_id INTEGER PRIMARY KEY)')
            for gene in GENES:
                for state in STATES:
                    db.execute("""INSERT OR IGNORE INTO selected
                        SELECT v.variation_id FROM full_data.variant_gene_symbols g
                        JOIN full_data.variants v ON v.variation_id=g.variation_id
                        WHERE g.gene_symbol=? COLLATE NOCASE AND v.clinical_state=?
                        ORDER BY v.variation_id LIMIT 50""",(gene,state))
            db.execute('INSERT OR IGNORE INTO selected SELECT variation_id FROM full_data.variants WHERE variation_id IN (17660,37565)')
            for name in ['variants','timelines','variant_gene_symbols']:
                db.execute(f'INSERT INTO {name} SELECT original.* FROM full_data.{name} original JOIN selected USING (variation_id)')
            db.execute('CREATE INDEX gene_symbol_lookup ON variant_gene_symbols (gene_symbol COLLATE NOCASE,variation_id)')
            metadata=dict(db.execute('SELECT key,value FROM full_data.metadata'))
            count=lambda table: db.execute(f'SELECT COUNT(*) FROM {table}').fetchone()[0]
            metadata.update({
                'variant_count':str(count('variants')),'timeline_count':str(count('timelines')),
                'gene_symbol_link_count':str(count('variant_gene_symbols')),
                'dataset_mode':'demo','demo_genes':json.dumps(GENES),
                'variant_source':'Subset of cached ClinVar latest retained records',
                'timeline_source':'Unmodified six-snapshot histories for selected variants',
                'scope':'Demo sample: up to 50 variants per selected gene and grouped classification, ordered by VariationID, plus example variants. This is not the full dataset or a representative prevalence sample.'})
            db.executemany('INSERT INTO metadata VALUES (?,?)',metadata.items())
            assert 0<count('variants')<=len(GENES)*len(STATES)*50+2
            assert count('variants')==count('timelines')
            assert db.execute('PRAGMA foreign_key_check').fetchone() is None
            assert db.execute('PRAGMA quick_check').fetchone()[0]=='ok'
            db.commit()
        temporary.replace(target)
        print(f'Demo ready: {metadata["variant_count"]} variants, {target.stat().st_size:,} bytes')
    except Exception:
        temporary.unlink(missing_ok=True)
        raise

if __name__=='__main__': build()
