import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from fastapi import HTTPException
from backend.app import main


class PartitionEndpointTests(unittest.TestCase):
    def test_latest_candidate_requires_matching_frozen_model(self):
        with tempfile.TemporaryDirectory() as root:
            self.fixture(root)
            directory=Path(root)/'reports/biological_corrected_20261109'
            directory.mkdir()
            (directory/'selection.json').write_text(json.dumps({'selected':{'accuracy':.91},'model_sha256':'latest'}))
            report={'status':'complete','model_sha256':'latest','scores':{'internal_test':{'accuracy':.9}}}
            path=directory/'partition_evaluation.json'
            path.write_text(json.dumps(report))
            with patch.object(main,'project_root',Path(root)):
                candidate=main.latest_training_evaluation()['corrected_candidate']
                self.assertTrue(candidate['test_evaluated'])
                self.assertFalse(candidate['deployed'])
                self.assertEqual(candidate['partition_evaluation'],report)
                report['model_sha256']='old'
                path.write_text(json.dumps(report))
                with self.assertRaises(HTTPException):
                    main.latest_training_evaluation()

    def fixture(self, root, evaluation=None):
        directory=Path(root)/'reports/biological_retuning_fresh_20261031'
        directory.mkdir(parents=True)
        (directory/'selection.json').write_text(json.dumps({
            'selected':{'accuracy':.81},'model_sha256':'correct'}))
        (directory/'data_audit.json').write_text('{}')
        if evaluation is not None:
            (directory/'partition_evaluation.json').write_text(json.dumps(evaluation))

    def test_missing_evaluation_is_not_reported_as_tested(self):
        with tempfile.TemporaryDirectory() as root:
            self.fixture(root)
            with patch.object(main,'project_root',Path(root)):
                result=main.latest_training_evaluation()
            self.assertFalse(result['test_evaluated'])
            self.assertIsNone(result['partition_evaluation'])

    def test_report_for_different_model_is_rejected(self):
        with tempfile.TemporaryDirectory() as root:
            self.fixture(root,{'status':'complete','model_sha256':'wrong'})
            with patch.object(main,'project_root',Path(root)):
                with self.assertRaises(HTTPException) as caught:
                    main.latest_training_evaluation()
            self.assertEqual(caught.exception.status_code,503)

    def test_complete_matching_report_is_served_without_deploying(self):
        report={'status':'complete','model_sha256':'correct','scores':{'training':{'accuracy':.82}}}
        with tempfile.TemporaryDirectory() as root:
            self.fixture(root,report)
            with patch.object(main,'project_root',Path(root)):
                result=main.latest_training_evaluation()
            self.assertTrue(result['test_evaluated'])
            self.assertFalse(result['deployed'])
            self.assertEqual(result['partition_evaluation'],report)
