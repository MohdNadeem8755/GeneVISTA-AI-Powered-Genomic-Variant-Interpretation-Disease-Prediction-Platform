import unittest
from src.investigate_class_errors import error_kind


class ClassErrorTests(unittest.TestCase):
    def test_boundaries_do_not_merge_distinct_errors(self):
        self.assertEqual(error_kind('Benign','Likely benign'),'within_pair')
        self.assertEqual(error_kind('Likely pathogenic','Pathogenic'),'within_pair')
        self.assertEqual(error_kind('Benign','Pathogenic'),'opposite_polarity')
        self.assertEqual(error_kind('Likely pathogenic','Uncertain significance'),'to_vus')
        self.assertEqual(error_kind('Benign','Benign'),'correct')
        self.assertEqual(error_kind('Uncertain significance','Benign'),'other')


if __name__=='__main__':
    unittest.main()
