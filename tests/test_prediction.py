import os
import unittest
from unittest.mock import patch
import numpy as np
from fastapi.testclient import TestClient
from backend.app.main import app
from backend.app.database import Database
from backend.app.prediction import predict,bundle
from src.variant_features import variant_features,parse_change,CLASSES

class PredictionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.row=Database('data/demo/genevista.sqlite').detail(17660)['latest_details']
    def test_exact_five_probabilities_and_one_winner(self):
        result=predict(self.row,17660)
        self.assertEqual(result['status'],'available')
        self.assertEqual([x['class_name'] for x in result['probabilities']],CLASSES)
        self.assertAlmostEqual(sum(x['probability'] for x in result['probabilities']),1,places=10)
        self.assertTrue(all(0<=x['probability']<=1 for x in result['probabilities']))
        self.assertEqual(result['predicted_class'],max(result['probabilities'],key=lambda x:x['probability'])['class_name'])
        self.assertEqual(result['cohort_membership'],'test')
        self.assertEqual(result['disagrees_with_recorded'],result['predicted_class']!=self.row['clinical_significance'])
    def test_labels_and_review_metadata_do_not_change_scores(self):
        changed=dict(self.row,clinical_significance='Benign',clinical_state='benign',review_status='fake',number_submitters_raw='999999',variation_id=99999)
        self.assertEqual(variant_features(self.row),variant_features(changed))
        self.assertEqual(predict(self.row)['probabilities'],predict(changed)['probabilities'])
        self.assertEqual(predict(self.row)['diseases']['items'],predict(changed)['diseases']['items'])
    def test_sequence_change_affects_prediction(self):
        changed=dict(self.row,name=self.row['name'].replace('190T>G','190T>A').replace('Cys64Gly','Cys64Ser'))
        self.assertNotEqual(variant_features(self.row),variant_features(changed))
        self.assertNotEqual(predict(self.row)['probabilities'],predict(changed)['probabilities'])
        change=parse_change(self.row['name'])
        self.assertEqual(change['nucleotide']['reference'],'T')
        self.assertEqual(change['nucleotide']['alternate'],'G')
    def test_ranked_conditional_disease_estimates(self):
        result=predict(self.row)['diseases']
        self.assertEqual(result['status'],'available')
        values=[row['match_probability'] for row in result['items']]
        self.assertEqual(values,sorted(values,reverse=True))
        self.assertEqual([row['rank'] for row in result['items']],list(range(1,len(values)+1)))
        self.assertTrue(all(0<=value<=1 for value in values))
        self.assertIn('not patient disease probabilities',result['probability_scope'])
    def test_unseen_gene_needs_evidence(self):
        self.assertEqual(predict(dict(self.row,gene_symbol='UNSEEN_GENE_000'))['status'],'insufficient_evidence')
        self.assertEqual(predict(None)['status'],'insufficient_evidence')
    def test_disjoint_training_calibration_and_test_ids(self):
        ids=bundle()['cohort_ids']
        for a,b in [('train','calibration'),('train','test'),('calibration','test')]:
            self.assertEqual(len(np.intersect1d(ids[a],ids[b])),0)
    def test_api_exact_recorded_and_research_prediction(self):
        with patch.dict(os.environ,{'GENEVISTA_AUTH_ENABLED':'false'}),TestClient(app) as client:
            result=client.get('/api/variants/search?q=17660').json()['items'][0]
            self.assertEqual(result['latest_details']['clinical_significance'],'Pathogenic')
            self.assertIn(result['prediction']['predicted_class'],CLASSES)
            detail=client.get('/api/variants/17660').json()
            self.assertTrue(detail['prediction_support']['clinical_class_probability'])
            self.assertFalse(detail['prediction_support']['disease_probability'])
            self.assertTrue(detail['prediction_support']['disease_annotation_match_probability'])
            self.assertEqual(client.get('/api/model/prediction-evaluation').status_code,200)
        with patch.dict(os.environ,{'GENEVISTA_AUTH_ENABLED':'true'}),TestClient(app) as client:
            self.assertEqual(client.get('/api/model/prediction-evaluation').status_code,401)
