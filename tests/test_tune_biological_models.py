import unittest
import numpy as np
import pandas as pd
from src.tune_biological_models import predict_bundle


class FixedClassifier:
    def __init__(self, classes, probabilities):
        self.classes_=np.asarray(classes)
        self.probabilities=np.asarray(probabilities)
    def predict_proba(self, frame):
        return np.tile(self.probabilities,(len(frame),1))


class TuneBiologicalTests(unittest.TestCase):
    def test_bundle_aligns_classes_and_handles_missing_category(self):
        bundle={'classes':np.asarray(['a','b']),'weights':[1,3],'models':[
            {'classifier':FixedClassifier(['a','b'],[.8,.2]),'columns':['cat'],'categorical':['cat']},
            {'classifier':FixedClassifier(['b','a'],[.6,.4]),'columns':['cat'],'categorical':['cat']}]}
        result=predict_bundle(bundle,pd.DataFrame({'cat':[None,'x']}))
        np.testing.assert_allclose(result,[[.5,.5],[.5,.5]])


if __name__=='__main__':unittest.main()
