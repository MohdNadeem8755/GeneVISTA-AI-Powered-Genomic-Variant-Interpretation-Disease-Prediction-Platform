import unittest
import numpy as np
from src.ann_annotation import (transform_inputs, preprocessor, class_weights,
                                development_split, objective, normalized_probabilities)
from src.variant_features import CLASSES


class ANNAnnotationTests(unittest.TestCase):
    def test_probability_roundoff_only(self):
        p = np.array([[.1, .2, .7]], dtype=np.float32)
        corrected = normalized_probabilities(p)
        self.assertAlmostEqual(float(corrected.sum()), 1, places=15)
        self.assertEqual(corrected.argmax(), p.argmax())
        with self.assertRaises(ValueError):
            normalized_probabilities([[.1, .2, .6]])

    def test_missing_population_is_not_zero(self):
        x = np.array([[np.nan, 0], [0, 1], [.1, 2]], dtype=np.float32)
        result = transform_inputs(x, ('gnomAD4.1_joint_AF',))
        self.assertTrue(np.isnan(result[0, 0]))
        self.assertEqual(result[1, 0], -9)
        self.assertAlmostEqual(result[2, 0], -1)
        self.assertEqual(x[1, 0], 0)
        with self.assertRaises(ValueError):
            transform_inputs([[1.1, 0]], ('gnomAD_AF',))

    def test_training_imputation_and_unknown_gene(self):
        train = np.array([[1., np.nan, 0], [3., np.nan, 1], [np.nan, np.nan, 0]])
        validation = np.array([[100., 2., 99], [np.nan, np.nan, np.nan]])
        prep = preprocessor(2)
        a = prep.fit_transform(train)
        b = prep.transform(validation)
        self.assertEqual(a.shape[1], b.shape[1])
        self.assertTrue(np.isfinite(b).all())
        self.assertEqual(prep.named_transformers_['numeric'][0].statistics_[0], 2)
        np.testing.assert_array_equal(b[:, -2:], 0)

    def test_balanced_class_mass(self):
        y = np.array([c for i, c in enumerate(CLASSES) for _ in range(i + 1)])
        weights = class_weights(y, 1)
        mass = [weights[y == c].sum() for c in CLASSES]
        np.testing.assert_allclose(mass, mass[0])
        self.assertAlmostEqual(float(weights.mean()), 1, places=6)

    def test_no_test_or_calibration_access(self):
        for split in ('test', 'calibration', 'excluded'):
            with self.assertRaisesRegex(ValueError, 'only read'):
                development_split(None, split, (), {}, 1)

    def test_objective_requires_both_metrics(self):
        a = dict(accuracy=.95, balanced_accuracy=.6, macro_f1=.6)
        b = dict(accuracy=.8, balanced_accuracy=.75, macro_f1=.7)
        self.assertGreater(objective(b), objective(a))
