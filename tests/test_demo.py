import os
import sqlite3
import unittest
from unittest.mock import patch
from pathlib import Path
from fastapi.testclient import TestClient
from backend.app.database import Database
from backend.app.main import app
from deployment.start import provision

ROOT=Path(__file__).resolve().parents[1]
DEMO=ROOT/'data/demo/genevista.sqlite'

class DemoTests(unittest.TestCase):
    def test_size_integrity_and_scope(self):
        self.assertLess(DEMO.stat().st_size,10*1024*1024)
        database=Database(DEMO)
        stats=database.stats()
        self.assertEqual(stats['dataset_mode'],'demo')
        self.assertEqual(len(stats['demo_genes']),20)
        with database.connect() as c:
            self.assertEqual(c.execute('PRAGMA quick_check').fetchone()[0],'ok')
            self.assertIsNone(c.execute('PRAGMA foreign_key_check').fetchone())
            self.assertEqual(stats['variant_count'],c.execute('SELECT count(*) FROM variants').fetchone()[0])
            self.assertEqual(stats['timeline_count'],stats['variant_count'])
    def test_records_preserve_source(self):
        source=ROOT/'data/app/genevista.sqlite'
        if not source.exists():self.skipTest('Full local source not installed')
        with Database(source).connect() as original, Database(DEMO).connect() as demo:
            for table in ('variants','timelines'):
                for row in demo.execute(f'SELECT * FROM {table}'):
                    expected=original.execute(f'SELECT * FROM {table} WHERE variation_id=?',(row['variation_id'],)).fetchone()
                    self.assertEqual(tuple(row),tuple(expected))
    def test_demo_api_rankings_and_unknown_record(self):
        with patch('backend.app.main.database',Database(DEMO)),patch.dict(os.environ,{'GENEVISTA_AUTH_ENABLED':'false'}),TestClient(app) as client:
            self.assertEqual(client.get('/health').status_code,200)
            self.assertEqual(client.get('/api/stats').json()['dataset_mode'],'demo')
            self.assertTrue(client.get('/api/variants/search?q=TP53').json()['items'])
            self.assertEqual(client.get('/api/variants/37565/ranking').json()['status'],'available')
            self.assertEqual(client.get('/api/variants/17660/ranking').json()['status'],'ineligible')
            self.assertEqual(client.get('/api/variants/9223372036854775807').status_code,404)
    def test_startup_needs_no_download_or_disk(self):
        with patch.dict(os.environ,{'GENEVISTA_DB':str(DEMO),'GENEVISTA_DB_URL':'','GENEVISTA_DB_SHA256':'','GENEVISTA_AUTH_ENABLED':'false'}),patch('urllib.request.urlopen',side_effect=AssertionError('Demo must not download data')):
            provision()
