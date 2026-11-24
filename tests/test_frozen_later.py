import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
import json
from src import evaluate_frozen_later


class FrozenLaterTests(unittest.TestCase):
    def test_freeze_occurs_before_labels_and_is_non_overwritable(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            report=root/'reports/biological_tuning_20261030'
            report.mkdir(parents=True)
            (report/'selected.joblib').write_bytes(b'test candidate')
            (report/'selection.json').write_text(json.dumps({
                'test_evaluated': False,
                'production_model_changed': False
            }), encoding='utf-8')
            target=root/'run'
            with patch.object(evaluate_frozen_later, 'ROOT', root):
                evaluate_frozen_later.freeze(target,'2026-10-26')
            text=(target/'frozen_protocol.json').read_text()
            self.assertIn('test_read_count',text)
            self.assertIn('candidate_sha256',text)
            with patch.object(evaluate_frozen_later, 'ROOT', root):
                with self.assertRaises(FileExistsError):
                    evaluate_frozen_later.freeze(target,'2026-10-26')


if __name__=='__main__':unittest.main()
