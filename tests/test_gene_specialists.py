import unittest
import sqlite3
import numpy as np
from src.train_gene_specialists import load_rows,weights
from src.variant_features import CLASSES

class SpecialistTests(unittest.TestCase):
    def test_held_out_partitions_rejected_before_query(self):
        for split in ('test','calibration','reserved'):
            with self.assertRaises(ValueError):load_rows(None,['BRCA1'],split)
    def test_exact_gene_and_split_filter(self):
        with sqlite3.connect(':memory:') as con:
            con.executescript('CREATE TABLE variants(chrom TEXT,pos INTEGER,ref TEXT,alt TEXT,label TEXT,gene TEXT,name TEXT,variant_type TEXT,split TEXT); CREATE TABLE annotations(chrom TEXT,pos INTEGER,ref TEXT,alt TEXT,scores TEXT);')
            for i,(gene,split) in enumerate((('BRCA1','train'),('BRCA1','test'),('BRCA2','train'))):
                con.execute('INSERT INTO variants VALUES (?,?,?,?,?,?,?,?,?)',('1',i,'A','G','Benign',gene,'name','SNV',split))
                con.execute('INSERT INTO annotations VALUES (?,?,?,?,?)',('1',i,'A','G','{}'))
            rows=load_rows(con,['BRCA1'],'train')
            self.assertEqual(len(rows),1)
            self.assertEqual(rows[0][2],'BRCA1')
    def test_weights_only_from_training_class_counts(self):
        labels=np.array(CLASSES+[CLASSES[0]]*9)
        np.testing.assert_equal(weights(labels,0),np.ones(len(labels)))
        self.assertLess(weights(labels,.5)[0],weights(labels,.5)[1])
        with self.assertRaises(ValueError):weights(np.array(['Benign']),.5)

class ConstantFeatureTests(unittest.TestCase):
    def test_mask_fitted_on_training_only_preserves_missing(self):
        from src.train_gene_specialists import NonconstantFeatures
        train=np.array([[np.nan,1.,2.],[np.nan,1.,3.],[np.nan,1.,np.nan]])
        selector=NonconstantFeatures().fit(train)
        np.testing.assert_array_equal(selector.keep_,[False,False,True])
        result=selector.transform(np.array([[7.,9.,np.nan]]))
        self.assertEqual(result.shape,(1,1))
        self.assertTrue(np.isnan(result[0,0]))
