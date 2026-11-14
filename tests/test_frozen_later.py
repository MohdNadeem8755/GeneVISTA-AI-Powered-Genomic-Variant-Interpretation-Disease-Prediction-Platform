import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
from src.evaluate_frozen_later import freeze


class FrozenLaterTests(unittest.TestCase):
    def test_freeze_occurs_before_labels_and_is_non_overwritable(self):
        with tempfile.TemporaryDirectory(dir='D:/GeneVISTA/.runtime') as tmp:
            target=Path(tmp)/'run'
            freeze(target,'2026-10-26')
            text=(target/'frozen_protocol.json').read_text()
            self.assertIn('test_read_count',text)
            self.assertIn('candidate_sha256',text)
            with self.assertRaises(FileExistsError):freeze(target,'2026-10-26')


if __name__=='__main__':unittest.main()
