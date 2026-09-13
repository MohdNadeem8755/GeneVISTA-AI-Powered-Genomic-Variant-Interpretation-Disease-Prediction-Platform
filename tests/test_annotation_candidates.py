import json
import unittest
import numpy as np
from src.train_annotation_candidates import matrix, INDEPENDENT, SEQUENCE


class AnnotationCandidateTests(unittest.TestCase):
    def test_feature_matrix_preserves_missing_and_ignores_label_and_identity(self):
        scores = json.dumps({'gnomAD4.1_joint_AF': None, 'CADD_phred': 15,
                             'clinvar_clnsig': 'Pathogenic'})
        rows = [('Pathogenic', scores, 'BRCA1', 'NM_1:c.1A>G', 'SNV'),
                ('Benign', scores, 'TP53', 'NM_1:c.1A>G', 'SNV')]
        features = ('gnomAD4.1_joint_AF', 'CADD_phred') + SEQUENCE
        values = matrix(rows, features)
        np.testing.assert_allclose(values[0], values[1], equal_nan=True)
        self.assertTrue(np.isnan(values[0, 0]))
        self.assertEqual(values[0, 1], 15)

    def test_ablation_excludes_supervised_pathogenicity_scores(self):
        self.assertFalse(set(INDEPENDENT) & {'REVEL_score', 'CADD_phred',
                         'AlphaMissense_score', 'Polyphen2_HVAR_rankscore',
                         'SIFT_converted_rankscore'})


if __name__ == '__main__':
    unittest.main()
