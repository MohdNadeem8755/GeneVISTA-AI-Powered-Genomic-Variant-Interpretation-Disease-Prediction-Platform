import hashlib
import os
import unittest
from unittest.mock import patch
import numpy as np
from scipy.sparse import csr_matrix
from fastapi.testclient import TestClient
from backend.app.main import app, database
from backend.app import auth
from backend.app.ranking import rank, bundle
from src.disease_ranking import csr32, feature_dict

class IntegrationTests(unittest.TestCase):
    def setUp(self):
        self.env=patch.dict(os.environ,{'GENEVISTA_AUTH_ENABLED':'false'})
        self.env.start()
        self.client=TestClient(app)
    def tearDown(self):
        self.client.close();self.env.stop()
    def test_sparse_64_bit_regression(self):
        x=csr_matrix([[0.,2.],[1.,0.]])
        x.indices=x.indices.astype(np.int64);x.indptr=x.indptr.astype(np.int64)
        converted=csr32(x)
        self.assertEqual(converted.indices.dtype,np.int32)
        self.assertEqual(converted.indptr.dtype,np.int32)
        np.testing.assert_array_equal(converted.toarray(),[[0,2],[1,0]])
    def test_model_feature_parity_and_eligibility(self):
        with database.connect() as c:
            row=dict(c.execute("SELECT * FROM variants WHERE gene_symbol='BRCA1' AND clinical_state='vus' LIMIT 1").fetchone())
        result=rank(row)
        self.assertEqual(result['status'],'available')
        model=bundle()
        x=csr32(model['vectorizer'].transform([feature_dict(dict(row,number_submitters=row['number_submitters_raw']))]))
        scores={d:float(clf.decision_function(x)[0]) for d,clf in zip(model['metadata']['classes'],model['classifiers'])}
        self.assertEqual({r['disease_id']:r['score'] for r in result['items']},scores)
        self.assertEqual([r['score'] for r in result['items']],sorted(scores.values(),reverse=True))
        self.assertEqual(rank(dict(row,clinical_state='pathogenic'))['status'],'ineligible')
        self.assertEqual(rank(dict(row,gene_symbol='GENE_NOT_IN_TRAINING'))['status'],'unsupported_gene')
    def test_api_data_and_metrics(self):
        self.assertEqual(self.client.get('/health').status_code,200)
        self.assertEqual(self.client.get('/api/variants/17660/ranking').json()['status'],'ineligible')
        response=self.client.get('/api/model/evaluation')
        self.assertEqual(response.status_code,200)
        rows=response.json()['items']
        self.assertEqual(len(rows),2*len(bundle()['classifiers']))
        for row in rows:
            if row['positive_events']==0:self.assertIsNone(row['average_precision'])
        self.assertEqual(self.client.get('/api/variants/0/ranking').status_code,422)
        self.assertEqual(self.client.get('/api/variants/9223372036854775807/ranking').status_code,404)
    def test_login_protects_api_logout_revokes_and_expiry(self):
        salt=bytes(range(16));password='test-only-long-password'
        digest=hashlib.scrypt(password.encode(),salt=salt,n=16384,r=8,p=1).hex()
        headers={'X-GeneVISTA-Request':'1'}
        with patch.dict(os.environ,{'GENEVISTA_AUTH_ENABLED':'true','GENEVISTA_USERNAME':'researcher','GENEVISTA_PASSWORD_HASH':salt.hex()+':'+digest,'GENEVISTA_COOKIE_SECURE':'false'}):
            self.assertEqual(self.client.get('/api/stats').status_code,401)
            self.assertEqual(self.client.get('/api/model/evaluation').status_code,401)
            self.assertEqual(self.client.get('/health').status_code,200)
            self.assertEqual(self.client.post('/api/auth/login',json={'username':'researcher','password':password}).status_code,403)
            self.assertEqual(self.client.post('/api/auth/login',headers=headers,json={'username':'researcher','password':'incorrect'}).status_code,401)
            response=self.client.post('/api/auth/login',headers=headers,json={'username':'researcher','password':password})
            self.assertEqual(response.status_code,200)
            self.assertIn('HttpOnly',response.headers['set-cookie'])
            token=self.client.cookies.get(auth.COOKIE)
            self.assertEqual(self.client.get('/api/stats').status_code,200)
            self.assertEqual(self.client.post('/api/auth/logout',headers=headers,json={}).status_code,200)
            self.client.cookies.set(auth.COOKIE,token)
            self.assertEqual(self.client.get('/api/stats').status_code,401)
            self.client.cookies.clear()
            self.client.post('/api/auth/login',headers=headers,json={'username':'researcher','password':password})
            with patch('backend.app.auth.time.time',return_value=10**12):
                self.assertEqual(self.client.get('/api/stats').status_code,401)
    def test_login_rate_limit_and_secure_cookie(self):
        auth._attempts.clear()
        with patch.dict(os.environ,{'GENEVISTA_AUTH_ENABLED':'true','GENEVISTA_USERNAME':'r','GENEVISTA_PASSWORD_HASH':'00:00'}):
            for _ in range(10):
                self.client.post('/api/auth/login',headers={'X-GeneVISTA-Request':'1'},json={'username':'r','password':'bad'})
            self.assertEqual(self.client.post('/api/auth/login',headers={'X-GeneVISTA-Request':'1'},json={'username':'r','password':'bad'}).status_code,429)
        auth._attempts.clear()
    def test_built_frontend(self):
        response=self.client.get('/')
        self.assertEqual(response.status_code,200)
        self.assertIn('GeneVISTA',response.text)

if __name__=='__main__':unittest.main()
