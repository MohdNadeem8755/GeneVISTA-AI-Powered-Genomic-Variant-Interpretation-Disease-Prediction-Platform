import contextlib,io,json,tempfile,unittest
from pathlib import Path
import numpy as np
import pandas as pd
import joblib,sklearn
from sklearn.pipeline import Pipeline
from sklearn.compose import ColumnTransformer
from sklearn.impute import SimpleImputer
from sklearn.dummy import DummyClassifier
from tests import test_feature_validation as fixtures
from src.evaluation_validation import resolve_bundle,load_evaluation_data,verify_predictions,sha256
ROOT=Path(__file__).resolve().parents[1]

class EvaluationValidationTests(unittest.TestCase):
    def test_prediction_corruption_detected(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'pred.csv.gz'
            expected=pd.DataFrame({'VariationID':['1','2'],'score_vus':[.1,.9]})
            expected.to_csv(path,index=False);verify_predictions(path,expected,chunksize=1)
            changed=expected.copy();changed.loc[1,'score_vus']=.8;changed.to_csv(path,index=False)
            with self.assertRaises(AssertionError):verify_predictions(path,expected,chunksize=1)

    def test_full_evaluation_without_positive_targets(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);folder,labels=fixtures.FeatureValidationTests().make_fixture(root)
            def execute(name):
                ns={'display':lambda x:None,'__name__':'__main__'}
                n=json.loads((ROOT/'notebooks'/name).read_text(encoding='utf-8'))
                with contextlib.redirect_stdout(io.StringIO()):
                    for c in n['cells']:
                        if c['cell_type']=='code':exec(''.join(c['source']).replace('PROJECT_ROOT = Path(r"D:\\GeneVISTA")',f'PROJECT_ROOT = Path({str(root)!r})'),ns)
                return ns
            ns=execute('03_feature_engineering.ipynb');config=ns['temporal_config']
            bundle=root/'models/temporal_v1/test_run';bundle.mkdir(parents=True)
            plan={'sklearn_version':sklearn.__version__,'feature_config':config,'test_interval':'2025-01_to_2026-01',
                  'test_files':{'baseline_features':'vus_features_2025-01.csv.gz','temporal_features':'vus_temporal_2025-01.csv.gz','labels':'labels_2025-01.csv.gz'},'expected_test_rows':2,
                  'review_sizes':[100,1000],'tie_break_seed':42,'primary_model':'primary','comparator_model':'comparator'}
            X=pd.DataFrame({c:[1.,2.,3.] for c in config['numeric_features']})
            for c in config['categorical_features']:X[c]='x'
            X=X[config['categorical_features']+config['numeric_features']]
            model=Pipeline([('features',ColumnTransformer([('numeric',SimpleImputer(),config['numeric_features'])])),('classifier',DummyClassifier())]).fit(X,config['target_classes'])
            for name in ('primary','comparator'):joblib.dump(model,bundle/f'{name}.joblib')
            (bundle/'evaluation_plan.json').write_text(json.dumps(plan))
            source=folder/'model_features/vus_features_2022-01.csv.gz'
            (bundle/'development_inputs.json').write_text(json.dumps([{'file':str(source.relative_to(root)),'sha256':sha256(source)}]))
            (root/'reports').mkdir();(root/'reports/training04_validation.json').write_text(json.dumps({'status':'passed','evaluation_bundle':'models/temporal_v1/test_run'}))
            self.assertEqual(resolve_bundle(root),bundle)
            result=execute('05_model_evaluation.ipynb')
            self.assertTrue(result['primary_metric']['mean_reclassification_average_precision'].isna().all())
            self.assertEqual(len(result['test_review_lists']),4)
            self.assertTrue(result['evaluation_record']['prediction_artifacts_verified'])
            self.assertFalse(result['evaluation_record']['fitting_performed'])
            plan['test_interval']='2024-01_to_2025-01'
            with self.assertRaisesRegex(ValueError,'interval'):load_evaluation_data(folder,plan)
