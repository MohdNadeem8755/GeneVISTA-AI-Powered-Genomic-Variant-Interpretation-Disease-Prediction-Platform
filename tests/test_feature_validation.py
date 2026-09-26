import contextlib
import io
import json
from pathlib import Path
import tempfile
import unittest
import pandas as pd
from src.feature_validation import build_baseline_features, load_label_manifest, load_labels, previous_snapshots

ROOT = Path(__file__).resolve().parents[1]

class FeatureValidationTests(unittest.TestCase):
    def test_missing_invalid_and_future_values(self):
        rows = pd.DataFrame({
            'VariationID': ['1','2','3','4','5'], 'Type': ['single nucleotide variant','.','-','NA','deletion'],
            'ReviewStatus': ['reviewed','none','',pd.NA,'reviewed'],
            'NumberSubmitters': ['2','-1','1.5','inf',''],
            'LastEvaluated': ['2021-12-31','2022-01-02','bad','-','2021-12-31T12:00:00Z'],
        }, dtype='string', index=[4,6,8,10,12])
        result = build_baseline_features(rows, '2022-01')
        self.assertEqual(result['number_submitters'].iloc[0], 2)
        self.assertEqual(result['submitter_count_missing'].tolist(), [0,1,1,1,1])
        self.assertTrue(result['variant_type'].iloc[1:4].isna().all())
        self.assertEqual(result['days_since_evaluation'].iloc[0],1)
        self.assertEqual(result['days_since_evaluation'].iloc[4],.5)
        self.assertEqual(result['evaluation_date_after_reference'].tolist(),[0,1,0,0,0])
        self.assertEqual(result['evaluation_date_missing'].tolist(),[0,0,1,1,0])
        self.assertTrue(pd.isna(result['days_since_evaluation'].iloc[1]))
        with self.assertRaises(ValueError):
            build_baseline_features(pd.concat([rows,rows]),'2022-01')

    def make_fixture(self, root):
        folder = root/'data/processed'; labels = folder/'model_labels'; labels.mkdir(parents=True)
        snapshots = ['2022-01','2023-01','2024-01','2025-01','2026-01']
        manifest=[]
        for i, snapshot in enumerate(snapshots[:-1]):
            split=['train','train','validation','test'][i]
            frame=pd.DataFrame({'VariationID':['1','2'], 'clinical_state':['vus','vus'],
                'Type':['single nucleotide variant','.'], 'ReviewStatus':['reviewed','-'],
                'NumberSubmitters':[str(8-i),''], 'LastEvaluated':['2021-12-31','-']})
            frame.to_csv(folder/f'clinvar_variants_{snapshot}.csv.gz',index=False)
            label_file=f'labels_{snapshot}.csv.gz'
            pd.DataFrame({'VariationID':['1','2'],'feature_snapshot':[snapshot]*2,
                'outcome_snapshot':[snapshots[i+1]]*2,'split':[split]*2,'target':['stayed_vus']*2}).to_csv(labels/label_file,index=False)
            manifest.append(dict(file=label_file,split=split,feature_snapshot=snapshot,outcome_snapshot=snapshots[i+1],rows=2,stayed_vus=2,became_benign=0,became_pathogenic=0))
        pd.DataFrame(manifest).to_csv(labels/'label_manifest.csv',index=False)
        pd.DataFrame({'VariationID':['1','2'],**{s:['vus','vus'] for s in snapshots}}).to_csv(folder/'clinvar_state_timelines.csv.gz',index=False)
        return folder, labels

    def test_metadata_validation_and_history(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder, labels = self.make_fixture(Path(tmp))
            manifest = load_label_manifest(labels, folder)
            row = next(manifest.itertuples(index=False))
            self.assertEqual(len(load_labels(labels,row)),2)
            mapping=previous_snapshots(folder/'clinvar_state_timelines.csv.gz',['2022-01','2025-01'])
            self.assertEqual(mapping,{'2022-01':None,'2025-01':'2024-01'})
            manifest.loc[0,'outcome_snapshot']='2026-01'
            manifest.to_csv(labels/'label_manifest.csv',index=False)
            with self.assertRaisesRegex(ValueError,'extend beyond'):
                load_label_manifest(labels,folder)
            manifest.loc[0,'outcome_snapshot']='2023-01'
            manifest.to_csv(labels/'label_manifest.csv',index=False)
            label=pd.read_csv(labels/row.file,dtype='string');label.loc[0,'split']='test';label.to_csv(labels/row.file,index=False)
            with self.assertRaisesRegex(ValueError,'split differs'):
                load_labels(labels,row)

    def test_notebook_clean_run_and_stale_outputs(self):
        notebook=json.loads((ROOT/'notebooks/03_feature_engineering.ipynb').read_text(encoding='utf-8'))
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);folder, labels=self.make_fixture(root)
            for attempt in range(2):
                if attempt:
                    for subdir, pattern in [('model_features','*.csv.gz'),('temporal_features_v1','*.csv.gz')]:
                        for p in (folder/subdir).glob(pattern):p.write_bytes(b'stale')
                    for p in (folder/'model_features').glob('*.json'):p.write_text('{}')
                ns={'display':lambda x: None,'__name__':'__main__'}
                with contextlib.redirect_stdout(io.StringIO()):
                    for c in notebook['cells']:
                        if c['cell_type']=='code':
                            source=''.join(c['source']).replace('PROJECT_ROOT = Path(r"D:\\GeneVISTA")',f'PROJECT_ROOT = Path({str(root)!r})')
                            exec(compile(source,'fixture_notebook03','exec'),ns)
                temporal=pd.read_csv(folder/'temporal_features_v1/vus_temporal_2023-01.csv.gz')
                self.assertEqual(temporal['submitter_count_change'].iloc[0],-1)
                self.assertEqual(temporal['submitter_change_missing'].iloc[1],1)
                baseline=pd.read_csv(folder/'model_features/vus_features_2022-01.csv.gz')
                self.assertTrue(pd.isna(baseline['variant_type'].iloc[1]))
                self.assertFalse(list(folder.rglob('*.partial*')))
