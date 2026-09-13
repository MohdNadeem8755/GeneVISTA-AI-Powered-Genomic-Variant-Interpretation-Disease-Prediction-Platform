import json
import tempfile
import unittest
from pathlib import Path
from src.audit_annotation_labels import audit, quality
from src.annotation_pipeline import connect, put_meta

class LabelAuditTests(unittest.TestCase):
    def test_review_tiers_require_exact_supported_labels(self):
        self.assertEqual(quality('reviewed by expert panel'),'expert')
        self.assertEqual(quality('criteria provided, multiple submitters, no conflicts'),'multiple_no_conflicts')
        self.assertEqual(quality('criteria provided, conflicting classifications'),'other')
    def test_excludes_test_and_calibration_from_support(self):
        with tempfile.TemporaryDirectory(dir=Path(__file__).resolve().parents[1]/'.runtime') as tmp:
            root=Path(tmp); db=root/'audit.sqlite'
            with connect(db) as con:
                con.executescript('CREATE TABLE metadata(key TEXT PRIMARY KEY,value TEXT); CREATE TABLE variants(chrom TEXT,pos INTEGER,ref TEXT,alt TEXT,split TEXT,label TEXT,review_status TEXT,gene TEXT); CREATE TABLE annotations(chrom TEXT,pos INTEGER,ref TEXT,alt TEXT);')
                put_meta(con,'annotation_status','ready')
                for pos,split in enumerate(('train','validation','test','calibration')):
                    con.execute('INSERT INTO variants VALUES (?,?,?,?,?,?,?,?)',('1',pos,'A','G',split,'Benign','reviewed by expert panel','BRCA1'))
                    con.execute('INSERT INTO annotations VALUES (?,?,?,?)',('1',pos,'A','G'))
            result=audit(db,root/'report')
            self.assertEqual(sum(row['count'] for row in result['counts']),2)
            self.assertEqual({r['split'] for r in result['counts']},{'train','validation'})
            self.assertEqual(result['eligible_genes'],[])
