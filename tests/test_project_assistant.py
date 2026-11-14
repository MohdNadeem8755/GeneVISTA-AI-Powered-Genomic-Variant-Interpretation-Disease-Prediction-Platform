import unittest
from unittest.mock import patch
from fastapi import FastAPI
from fastapi.testclient import TestClient
from backend.app.assistant import router


class AssistantTests(unittest.TestCase):
    def setUp(self):
        app = FastAPI()
        app.include_router(router)
        self.client = TestClient(app)

    def ask(self, question, mode='guide'):
        return self.client.post('/api/assistant/chat', json={'question': question, 'mode': mode},
                                headers={'X-GeneVISTA-Request': '1'})

    def test_saved_scores_labelled_validation(self):
        result = self.ask('latest scores').json()
        self.assertIn('VALIDATION', result['answer'])
        self.assertIn('Not deployed', result['answer'])
        self.assertTrue(result['sources'])

    def test_unknown_not_fabricated(self):
        self.assertIn('cannot answer', self.ask('xyzxyz').json()['answer'])

    def test_missing_provider_and_input_validation(self):
        with patch.dict('os.environ', {}, clear=True):
            self.assertEqual(self.ask('scores', 'online').status_code, 503)
        self.assertEqual(self.ask('   ').status_code, 422)
        self.assertEqual(self.ask('x'*2001).status_code, 422)
        self.assertEqual(self.client.post('/api/assistant/chat', json={'question':'scores'}).status_code, 403)
