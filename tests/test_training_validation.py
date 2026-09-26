import contextlib
import io
import json
from pathlib import Path
from types import SimpleNamespace
import tempfile
import unittest
from unittest.mock import Mock
import numpy as np
import pandas as pd
from sklearn.exceptions import ConvergenceWarning
from src.training_validation import load_split, ranking_metrics, safe_ratio, select_primary, fit_converged, validate_config
from tests import test_feature_validation as feature_tests

ROOT=Path(__file__).resolve().parents[1]

class TrainingValidationTests(unittest.TestCase):
    def test_zero_positive_ranking_and_selection(self):
        prevalence,ap,lift=ranking_metrics([False,False],[.1,.2])
        self.assertEqual(prevalence,0)
        self.assertTrue(np.isnan(ap) and np.isnan(lift))
        self.assertTrue(np.isnan(safe_ratio(0,0)))
        with self.assertRaises(ValueError):ranking_metrics([True],[np.nan])
        ranking=pd.DataFrame({'model':['boost','boost','logistic','logistic'],'target':['became_benign','became_pathogenic']*2,'average_precision':[.1,.2,.5,.6]})
        name,table=select_primary(ranking,{'boost':None,'logistic':None})
        self.assertEqual(name,'logistic')
        with self.assertRaises(ValueError):select_primary(ranking.iloc[:3],{'logistic':None})

    def test_reject_unconverged_and_test_access(self):
        model=Mock();model.fit.side_effect=ConvergenceWarning('did not converge')
        with self.assertRaisesRegex(RuntimeError,'did not converge'):fit_converged(model,[],[])
        with self.assertRaisesRegex(ValueError,'held-out'):load_split('.',None,None,None,'test')

    def test_full_notebook_and_invalid_inputs(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);folder,label_folder=feature_tests.FeatureValidationTests().make_fixture(root)
            manifest=pd.read_csv(label_folder/'label_manifest.csv',dtype='string')
            for row in manifest.itertuples(index=False):
                path=label_folder/row.file; labels=pd.read_csv(path,dtype='string')
                extra=labels.iloc[[0]].copy();extra['VariationID']='3';labels=pd.concat([labels,extra],ignore_index=True)
                labels['target']=['became_benign','became_pathogenic','stayed_vus'];labels.to_csv(path,index=False)
                source=folder/f'clinvar_variants_{row.feature_snapshot}.csv.gz';data=pd.read_csv(source,dtype='string',keep_default_na=False)
                extra=data.iloc[[0]].copy();extra['VariationID']='3';data=pd.concat([data,extra],ignore_index=True);data.to_csv(source,index=False)
            manifest['rows']='3'
            for target in ('stayed_vus','became_benign','became_pathogenic'):manifest[target]='1'
            manifest.to_csv(label_folder/'label_manifest.csv',index=False)
            path=folder/'clinvar_state_timelines.csv.gz';timeline=pd.read_csv(path,dtype='string');extra=timeline.iloc[[0]].copy();extra['VariationID']='3';pd.concat([timeline,extra]).to_csv(path,index=False)
            for notebook in ('03_feature_engineering.ipynb','04_model_training.ipynb'):
                n=json.loads((ROOT/'notebooks'/notebook).read_text(encoding='utf-8'))
                ns={'display':lambda x:None,'__name__':'__main__'}
                with contextlib.redirect_stdout(io.StringIO()):
                    for i,c in enumerate(n['cells']):
                        if c['cell_type']=='code':
                            source=''.join(c['source']).replace('PROJECT_ROOT = Path(r"D:\\GeneVISTA")',f'PROJECT_ROOT = Path({str(root)!r})')
                            exec(compile(source,f'{notebook}:{i}','exec'),ns)
            plan=ns['evaluation_plan'];self.assertFalse(plan['refitted_on_validation'])
            self.assertTrue((ns['evaluation_bundle']/(plan['primary_model']+'.joblib')).exists())
            self.assertEqual(plan['expected_test_rows'],3)
            self.assertEqual(plan['validation_variants_seen_in_training'],3)
            bad=dict(ns['feature_config']);bad['numeric_features']=bad['numeric_features']+['outcome_snapshot']
            with self.assertRaises(ValueError):validate_config(bad)
            source=folder/'model_features'/ns['feature_manifest'].iloc[0].feature_file
            data=pd.read_csv(source,dtype='string',keep_default_na=False);data.loc[0,'VariationID']='999';data.to_csv(source,index=False)
            with self.assertRaisesRegex(ValueError,'IDs do not match'):
                load_split(folder,ns['feature_manifest'],ns['label_manifest'],ns['feature_config'],'train')
