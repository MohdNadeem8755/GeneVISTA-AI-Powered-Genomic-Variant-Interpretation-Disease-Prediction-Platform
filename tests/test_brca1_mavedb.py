import gzip
import json
from pathlib import Path
import unittest

from src.brca1_mavedb import assess_variant, parse_variant, score_classification

ROOT = Path(__file__).resolve().parents[1]


class Brca1MaveDbTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with gzip.open(ROOT / "data/processed/brca1_mavedb_v1.json.gz", "rt", encoding="utf-8") as handle:
            cls.bundle = json.load(handle)

    def test_transcript_and_score_mapping(self):
        record = {"gene_symbol": "BRCA1", "name": "NM_007294.4(BRCA1):c.190T>G (p.Cys64Gly)"}
        result = assess_variant(record, self.bundle)
        self.assertEqual(result["status"], "available")
        self.assertEqual(result["functional_class"], "non-functional")
        self.assertEqual(result["criterion_candidate"], "PS3_Strong")
        self.assertIsNone(result["classification"])
        self.assertIsNone(result["patient_disease_probability"])

    def test_scope_and_unmatched_rows(self):
        self.assertEqual(assess_variant({"gene_symbol": "TP53", "name": "x"}, self.bundle)["status"], "out_of_scope")
        self.assertEqual(assess_variant({"gene_symbol": "BRCA1", "name": "NM_007294.4(BRCA1):c.68_69del (p.Glu23fs)"}, self.bundle)["status"], "unsupported_variant")
        self.assertIsNone(parse_variant("NM_007294.4(BRCA1):c.68_69del (p.Glu23fs)"))

    def test_threshold_boundaries(self):
        thresholds = {"non_functional_upper": -1.328, "functional_upper": -0.748}
        self.assertEqual(score_classification(-1.329, thresholds)["functional_class"], "non-functional")
        self.assertEqual(score_classification(-1.0, thresholds)["functional_class"], "intermediate")
        self.assertEqual(score_classification(-0.748, thresholds)["functional_class"], "functional")


if __name__ == "__main__":
    unittest.main()
