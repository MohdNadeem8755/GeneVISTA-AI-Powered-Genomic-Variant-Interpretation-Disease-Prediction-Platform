from contextlib import closing
import gzip
import hashlib
import json
from pathlib import Path
import sqlite3
import tempfile
import unittest
from src.protein_lm_features import parse_features, extract, REQUIRED


def example():
    return dict(ESM1b_score='-2;-9',popEVE_score='-3',MutFormer_score='0.2;0.9',
                Ensembl_transcriptid='ENST1;ENST2',MANE='.;Select',MPC_score='1;4',VEST4_score='.1;.8')


class ProteinLMTests(unittest.TestCase):
    def test_direction_and_transcript_alignment(self):
        r=parse_features(example())
        self.assertEqual(r['ESM1b_score_min'],-9)
        self.assertEqual(r['ESM1b_score_max'],-2)
        self.assertEqual(r['ESM1b_score_median'],-5.5)
        self.assertEqual(r['MPC_score_mane'],4)
        row=example();row['MPC_score']='2'
        self.assertIsNone(parse_features(row)['MPC_score_mane'])
    def test_missing_scores_and_labels(self):
        row=example();row['popEVE_score']='.';row['clinvar_clnsig']='Pathogenic'
        r=parse_features(row);self.assertIsNone(r['popEVE_score_min'])
        row['clinvar_clnsig']='Benign';self.assertEqual(r,parse_features(row))
        row['ESM1b_score']='inf'
        with self.assertRaises(ValueError):parse_features(row)
    def test_import_exact_alleles_excludes_test_and_collapses_duplicates(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);source=root/'db.gz';cohort=root/'cohort.sqlite';output=root/'features.sqlite'
            row=example();rows=[]
            for pos,alt in [(1,'G'),(1,'G'),(2,'G'),(1,'T')]:
                item={**row,'#chr':'1','pos(1-based)':str(pos),'ref':'A','alt':alt}
                rows.append('\t'.join(item[k] for k in REQUIRED))
            with gzip.open(source,'wt') as f:f.write('\t'.join(REQUIRED)+'\n'+'\n'.join(rows)+'\n')
            md5=hashlib.md5(source.read_bytes()).hexdigest();checksum=root/'db.md5';checksum.write_text(md5)
            con=sqlite3.connect(cohort)
            con.executescript('CREATE TABLE variants(chrom TEXT,pos INT,ref TEXT,alt TEXT,split TEXT); CREATE TABLE annotations(chrom TEXT,pos INT,ref TEXT,alt TEXT); CREATE TABLE metadata(key TEXT,value TEXT);')
            con.executemany('INSERT INTO variants VALUES (?,?,?,?,?)',[('1',1,'A','G','train'),('1',2,'A','G','test')])
            con.executemany('INSERT INTO annotations VALUES (?,?,?,?)',[('1',1,'A','G'),('1',2,'A','G')])
            con.executemany('INSERT INTO metadata VALUES (?,?)',[('annotation_status',json.dumps('ready')),('annotation',json.dumps({'md5':md5}))]);con.commit();con.close()
            report=extract(source,checksum,cohort,output)
            self.assertEqual(report['matched_variants'],1);self.assertEqual(report['matched_source_rows'],2)
            self.assertEqual(report['test_rows_read'],0)
            with closing(sqlite3.connect(output)) as con:
                self.assertEqual(con.execute('SELECT chrom,pos,ref,alt FROM features').fetchall(),[('1',1,'A','G')])
            with self.assertRaises(ValueError):extract(source,checksum,cohort,output)
