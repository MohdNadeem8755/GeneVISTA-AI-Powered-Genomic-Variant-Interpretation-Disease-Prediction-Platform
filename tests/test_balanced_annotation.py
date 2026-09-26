import unittest
import numpy as np
from src.balanced_annotation import decision_indices,enhanced_matrix,read_split

class BalancedAnnotationTests(unittest.TestCase):
    def test_evaluation_uses_active_artifact_not_stale_sidecar(self):
        from pathlib import Path
        import tempfile
        from unittest.mock import patch
        from backend.app import annotation_prediction as api
        with tempfile.TemporaryDirectory() as temp:
            artifact=Path(temp)/'predictor.joblib';artifact.write_bytes(b'fixture')
            artifact.with_name('model_card.json').write_text('{"model_id":"old"}')
            metadata={'model_id':'active','test':{'accuracy':.76,'balanced_accuracy':.71}}
            with patch.object(api,'model_path',return_value=artifact),patch.object(api,'bundle',return_value={'metadata':metadata}):
                card=api.model_card()
            self.assertEqual(card['model_id'],'active')
            self.assertEqual(card['test'],metadata['test'])
            self.assertTrue(card['model_revision'])
            self.assertTrue(card['model_saved_at'])

    def test_balanced_decision_does_not_relabel_probabilities(self):
        p=np.array([[.8,.2],[.3,.7]]);original=p.copy()
        np.testing.assert_array_equal(decision_indices(p,[.9,.1],0),[0,1])
        np.testing.assert_array_equal(decision_indices(p,[.9,.1],1),[1,1])
        np.testing.assert_array_equal(p,original)
        with self.assertRaises(ValueError):decision_indices(p,[0,1],1)
        with self.assertRaises(ValueError):decision_indices([[np.nan,1]],[.5,.5],1)
    def test_test_sealed_and_gene_vocab_training_only(self):
        with self.assertRaisesRegex(ValueError,'sealed'):read_split(None,'test',[],{})
        rows=[('Pathogenic','{}','KNOWN','c.1A>G','SNV'),('Benign','{}','NEW','c.1A>G','SNV')]
        x=enhanced_matrix(rows,('transition',),{'KNOWN':0})
        self.assertEqual(x[0,-1],0);self.assertTrue(np.isnan(x[1,-1]))
        self.assertEqual(x[0,0],x[1,0])

    def test_fresh_test_excludes_legacy_sample_loci(self):
        import sqlite3,json
        from src.annotation_pipeline import bucket
        from src.variant_features import CLASSES
        con=sqlite3.connect(':memory:')
        con.executescript('CREATE TABLE variants(chrom TEXT,pos INT,ref TEXT,alt TEXT,label TEXT,gene TEXT,name TEXT,variant_type TEXT,split TEXT); CREATE TABLE annotations(chrom TEXT,pos INT,ref TEXT,alt TEXT,scores TEXT);')
        for pos in range(1,201):
            con.execute('INSERT INTO variants VALUES (?,?,?,?,?,?,?,?,?)',('1',pos,'A','G',CLASSES[pos%5],'GENE','c.1A>G','SNV','test'))
            con.execute('INSERT INTO annotations VALUES (?,?,?,?,?)',('1',pos,'A','G','{}'))
        x,y=read_split(con,'test',('transition',),{'GENE':0},test_authorized=True)
        expected=sum(bucket(f'genevista-annotation-v2-sample|1|{pos}',10)!=0 for pos in range(1,201))
        self.assertEqual(len(y),expected)
        con.close()
    def test_api_preserves_probabilities_and_explains_balanced_decision(self):
        import tempfile,sqlite3,json,os
        from pathlib import Path
        from types import SimpleNamespace
        from unittest.mock import patch
        from backend.app import annotation_prediction as api
        from src.variant_features import CLASSES
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);database=root/'cohort.sqlite';model_path=root/'predictor.joblib';model_path.write_bytes(b'test')
            con=sqlite3.connect(database)
            con.executescript('CREATE TABLE metadata(key TEXT,value TEXT); CREATE TABLE variants(chrom TEXT,pos INT,ref TEXT,alt TEXT,variation_id TEXT,gene TEXT,name TEXT,variant_type TEXT,split TEXT); CREATE TABLE annotations(chrom TEXT,pos INT,ref TEXT,alt TEXT,scores TEXT);')
            con.executemany('INSERT INTO metadata VALUES (?,?)',[('annotation_status',json.dumps('ready')),('annotation',json.dumps({'md5':'test'}))])
            con.execute('INSERT INTO variants VALUES (?,?,?,?,?,?,?,?,?)',('1',1,'A','G','1','GENE','c.1A>G','SNV','validation'))
            con.execute('INSERT INTO annotations VALUES (?,?,?,?,?)',('1',1,'A','G','{}'));con.commit();con.close()
            probabilities=np.array([[.03,.04,.80,.10,.03]])
            classifier=SimpleNamespace(classes_=np.array(CLASSES),predict_proba=lambda x:probabilities)
            model={'classifier':classifier,'features':('transition',),'gene_vocab':{'GENE':0},'class_priors':np.array([.02,.02,.9,.05,.01]),'decision_exponent':1.,
                   'metadata':{'model_id':'annotation_prediction_v4','annotation':{'md5':'test'},'test':{'accuracy':.7,'balanced_accuracy':.7,'macro_f1':.5},'test_count':10,'scope':'research','upstream_overlap':'unknown'}}
            with patch.object(api,'model_path',return_value=model_path),patch.object(api,'bundle',return_value=model),patch.dict(os.environ,{'GENEVISTA_ANNOTATION_DB':str(database)}):
                result=api.predict({'name':'c.1A>G','gene_symbol':'GENE','clinical_significance':'Benign'},1)
            self.assertEqual(result['predicted_class'],'Benign')
            self.assertEqual(result['most_probable_class'],'Uncertain significance')
            self.assertAlmostEqual(result['top_probability'],.03)
            self.assertEqual([r['probability'] for r in result['probabilities']],probabilities[0].tolist())
            self.assertTrue(result['decision_policy'])
