from contextlib import closing
import json
from pathlib import Path
import sqlite3
import tempfile
import unittest
from src.train_protein_catboost import load_development
from src.protein_lm_features import FEATURES
from src.variant_features import CLASSES


class ProteinCatboostTests(unittest.TestCase):
    def test_test_labels_are_inaccessible(self):
        with self.assertRaisesRegex(ValueError,'development'):
            load_development(None,None,'test',(),1)
    def test_cached_feature_order_and_provenance(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);database=root/'base.sqlite';protein=root/'extra.sqlite'
            with closing(sqlite3.connect(database)) as con:
                con.executescript('CREATE TABLE variants(chrom TEXT,pos INT,ref TEXT,alt TEXT,label TEXT,gene TEXT,name TEXT,variant_type TEXT,split TEXT); CREATE TABLE annotations(chrom TEXT,pos INT,ref TEXT,alt TEXT,scores TEXT); CREATE TABLE metadata(key TEXT,value TEXT);')
                con.executemany('INSERT INTO metadata VALUES (?,?)',[('annotation_status',json.dumps('ready')),('annotation',json.dumps({'md5':'verified'}))])
                for i,label in enumerate(CLASSES):
                    con.execute('INSERT INTO variants VALUES (?,?,?,?,?,?,?,?,?)',('1',i+1,'A','G',label,'GENE','c.1A>G','single nucleotide variant','train'))
                    con.execute('INSERT INTO annotations VALUES (?,?,?,?,?)',('1',i+1,'A','G','{}'))
                con.commit()
            with closing(sqlite3.connect(protein)) as con:
                con.executescript('CREATE TABLE features(chrom TEXT,pos INT,ref TEXT,alt TEXT,scores TEXT); CREATE TABLE metadata(key TEXT,value TEXT);')
                con.execute('INSERT INTO metadata VALUES (?,?)',('protein_features',json.dumps({'status':'ready','md5':'verified'})))
                for i in range(5):con.execute('INSERT INTO features VALUES (?,?,?,?,?)',('1',i+1,'A','G',json.dumps({name:-i for name in FEATURES})))
                con.commit()
            x,y=load_development(database,protein,'train',('transition',),1)
            self.assertEqual(list(x.columns),['transition']+list(FEATURES)+['gene'])
            self.assertEqual(x['ESM1b_score_min'].tolist(),[0,-1,-2,-3,-4])
            self.assertEqual(y.tolist(),list(CLASSES))
            self.assertFalse(set(x.columns)&{'label','name','clinical_significance','variation_id'})
            with closing(sqlite3.connect(protein)) as con:
                con.execute('UPDATE metadata SET value=?',(json.dumps({'status':'ready','md5':'wrong'}),));con.commit()
            with self.assertRaisesRegex(ValueError,'provenance'):
                load_development(database,protein,'train',('transition',),1)
