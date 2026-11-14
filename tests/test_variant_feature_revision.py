import unittest
import json
import sqlite3
import tempfile
from pathlib import Path
import numpy as np
from src.variant_feature_revision import identity_features, revised_base_matrix
from src.train_annotation_candidates import matrix, INDEPENDENT, SEQUENCE
from src.variant_features import variant_features
from src.biological_features import load_development
from src.variant_feature_revision import CATEGORICAL
from src.tune_biological_models import predict_bundle
from catboost import CatBoostClassifier


class RevisedFeatureTests(unittest.TestCase):
    def test_revised_inputs_fit_and_predict_all_five_classes(self):
        labels=['Benign','Likely benign','Uncertain significance','Likely pathogenic','Pathogenic']
        rows=[(labels[i%5],'{}','GENE'+str(i%3),f'c.{100+i}G>T (p.Ala{34+i}Ser)','SNV','G','T') for i in range(30)]
        frame=revised_base_matrix(rows)
        model=CatBoostClassifier(iterations=3,depth=2,thread_count=1,verbose=False,allow_writing_files=False)
        model.fit(frame,[r[0] for r in rows],cat_features=list(CATEGORICAL))
        bundle={'models':[{'classifier':model,'columns':list(frame),'categorical':list(CATEGORICAL)}],
                'classes':labels,'weights':[1]}
        probabilities=predict_bundle(bundle,frame)
        self.assertEqual(probabilities.shape,(30,5))
        self.assertTrue(np.isfinite(probabilities).all())
        np.testing.assert_allclose(probabilities.sum(axis=1),1)

    def test_gene_names_do_not_create_sequence_events(self):
        for gene,event in [('INS','ins'),('WFS1','fs')]:
            name=f'NM_000000.1({gene}):c.100G>T (p.Ala34Ser)'
            # Preserve historical behavior for existing bundles; v2 corrects it.
            self.assertEqual(variant_features({'name':name,'gene_symbol':gene})['event='+event],1)
            self.assertEqual(identity_features(gene,name,'G','T')['event='+event],0)

    def test_actual_hgvs_events_are_preserved(self):
        r=identity_features('GENE','NM_1:c.100_101delinsAA (p.Arg34SerfsTer2)','G','T')
        for token in ('del','ins','fs'):
            self.assertEqual(r['event='+token],1)

    def test_different_substitutions_no_longer_collapse(self):
        rows=[('Benign','{}','GENE','c.100G>T (p.Ala34Ser)','SNV','G','T'),
              ('Benign','{}','GENE','c.100G>C (p.Ala34Pro)','SNV','G','C')]
        legacy=matrix([r[:5] for r in rows],tuple(INDEPENDENT)+SEQUENCE)
        self.assertTrue(np.array_equal(legacy[0],legacy[1],equal_nan=True))
        revised=revised_base_matrix(rows)
        self.assertNotEqual(revised.iloc[0]['genomic_base_change'],revised.iloc[1]['genomic_base_change'])
        self.assertNotEqual(revised.iloc[0]['amino_change'],revised.iloc[1]['amino_change'])

    def test_labels_cannot_change_features(self):
        row=('Benign','{}','GENE','c.100G>T (p.Ala34Ser)','SNV','G','T')
        changed=('Pathogenic',*row[1:])
        self.assertTrue(revised_base_matrix([row]).equals(revised_base_matrix([changed])))

    def test_missing_or_ambiguous_hgvs_is_explicit(self):
        r=identity_features('GENE','c.100G>T; c.101C>A','G','T')
        self.assertEqual(r['genomic_base_change'],'G>T')
        self.assertEqual(r['transcript_base_change'],'__MISSING__')
        self.assertEqual(r['amino_change'],'__MISSING__')

    def test_new_loader_profile_preserves_legacy_default_and_split_guard(self):
        with tempfile.TemporaryDirectory() as folder:
            cohort=Path(folder)/'cohort.sqlite'; cache=Path(folder)/'bio.sqlite'
            with sqlite3.connect(cohort) as con:
                con.executescript('CREATE TABLE metadata(key TEXT,value TEXT); CREATE TABLE variants(chrom TEXT,pos INTEGER,ref TEXT,alt TEXT,label TEXT,gene TEXT,name TEXT,variant_type TEXT,split TEXT); CREATE TABLE annotations(chrom TEXT,pos INTEGER,ref TEXT,alt TEXT,scores TEXT);')
                con.execute('INSERT INTO metadata VALUES (?,?)',('annotation',json.dumps({'md5':'fixture'})))
                con.execute('INSERT INTO variants VALUES (?,?,?,?,?,?,?,?,?)',('1',100,'G','T','Benign','INS','NM_1(INS):c.100G>T (p.Ala34Ser)','SNV','train'))
                con.execute('INSERT INTO annotations VALUES (?,?,?,?,?)',('1',100,'G','T','{}'))
            con.close()
            with sqlite3.connect(cache) as con:
                con.executescript('CREATE TABLE metadata(key TEXT,value TEXT); CREATE TABLE features(chrom TEXT,pos INTEGER,ref TEXT,alt TEXT,scores TEXT);')
                con.execute('INSERT INTO metadata VALUES (?,?)',('biological_features',json.dumps({'annotation_provenance':{'md5':'fixture'},'numeric':[],'categorical':[]})))
                con.execute('INSERT INTO features VALUES (?,?,?,?,?)',('1',100,'G','T','{}'))
            con.close()
            legacy,y,_=load_development(cohort,cache,'train')
            revised,y2,cats=load_development(cohort,cache,'train',feature_profile='variant_identity_v2')
            self.assertTrue(y.equals(y2))
            self.assertNotIn('genomic_base_change',legacy)
            self.assertEqual(legacy.iloc[0]['event=ins'],1)
            self.assertEqual(revised.iloc[0]['event=ins'],0)
            self.assertEqual(revised.iloc[0]['genomic_base_change'],'G>T')
            self.assertIn('gene_identity',cats)
            with self.assertRaises(ValueError):
                load_development(cohort,cache,'test',feature_profile='variant_identity_v2')
