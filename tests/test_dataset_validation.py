import gzip
from pathlib import Path
import tempfile
import unittest
from src.dataset_validation import table_check, ontology_check, validate_datasets

class DatasetValidationTests(unittest.TestCase):
    def test_missing_sources_cannot_pass(self):
        with tempfile.TemporaryDirectory() as folder:
            report = validate_datasets(folder)
        self.assertEqual(report['passed'], 0)
        self.assertEqual(report['expected'], 11)

    def test_shifted_clingen_header_and_separator(self):
        with tempfile.TemporaryDirectory() as folder:
            p = Path(folder)/'table.csv'
            p.write_text('metadata\nmore metadata\nGENE SYMBOL,MOI\n+++,+++\nBRCA1,AD\n')
            result = table_check(p, ('GENE SYMBOL','MOI'), delimiter=',', limit=None, header_first='GENE SYMBOL')
            self.assertEqual(result['rows_checked'], 1)

    def test_bad_schema_and_row(self):
        with tempfile.TemporaryDirectory() as folder:
            p = Path(folder)/'table.tsv'
            for content in ('a\ta\n1\t2\n','a\tb\n1\n','a\tb\n'):
                p.write_text(content)
                with self.assertRaises(ValueError): table_check(p, ('a','b'))

    def test_full_scan_catches_truncated_gzip(self):
        with tempfile.TemporaryDirectory() as folder:
            p = Path(folder)/'table.gz'
            p.write_bytes(gzip.compress(b'a\tb\n1\t2\n')[:-5])
            with self.assertRaises(EOFError): table_check(p, ('a','b'), limit=None)

    def test_sampling_and_missing_are_explicit(self):
        with tempfile.TemporaryDirectory() as folder:
            p = Path(folder)/'table.tsv'; p.write_text('a\tb\n1\t-\n2\t3\n')
            result = table_check(p, ('a','b'), limit=1)
            self.assertEqual(result['rows_checked'], 1)
            self.assertIn('only', result['scope'])
            self.assertEqual(result['missing_values_observed']['b'], 1)

    def test_ontology_invalid_term(self):
        with tempfile.TemporaryDirectory() as folder:
            p = Path(folder)/'hp.obo'; p.write_text('format-version: 1.2\n[Term]\nid: HP:0000001\n')
            with self.assertRaises(ValueError): ontology_check(p)
