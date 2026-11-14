import unittest
import numpy as np
from src.retrain_annotation_ensemble import AnnotationVote


class FixedModel:
    def __init__(self, classes, probabilities):
        self.classes_ = np.array(classes)
        self.probabilities = probabilities

    def predict_proba(self, X):
        return np.tile(self.probabilities, (len(X), 1))


class AnnotationVoteTests(unittest.TestCase):
    def test_aligns_class_orders_before_averaging(self):
        first = FixedModel(['a', 'b'], [.8, .2])
        second = FixedModel(['b', 'a'], [.6, .4])
        vote = AnnotationVote([first, second], [1, 3])
        np.testing.assert_allclose(vote.predict_proba([[0]]), [[.5, .5]])
        self.assertEqual(vote.predict([[0]])[0], 'a')


if __name__ == '__main__':
    unittest.main()
