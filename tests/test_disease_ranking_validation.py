import contextlib,io,json,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
import joblib
import pandas as pd
from src import disease_ranking as dr
from tests import test_feature_validation as fixtures
ROOT=Path(__file__).resolve().parents[1]

class DiseaseRankingValidationTests(unittest.TestCase):
    def test_training_only_vocabulary_and_reload(self):
        frames={}
        for split,start in [('train','2022-01'),('validation','2024-01'),('test','2025-01')]:
            frames[split]=pd.DataFrame([{'VariationID':str(i),'feature_snapshot':start,'split':split,
                'target':'became_pathogenic' if i<35 else 'stayed_vus',
                'diseases':['MONDO:0000001']+(['MONDO:0000002'] if split!='train' else []) if i<35 else [],
                'gene_symbol':'A' if i<35 else 'B','variant_type':'SNV','review_status':'reviewed','number_submitters':'2'} for i in range(80)])
        with tempfile.TemporaryDirectory() as tmp,contextlib.redirect_stdout(io.StringIO()):
            results=dr.train(frames,{},output_dir=tmp)
            model=joblib.load(Path(tmp)/'disease_ranker.joblib')
            self.assertEqual(model['metadata']['classes'],['MONDO:0000001'])
            self.assertEqual(len(results),2)
            self.assertTrue(results['average_precision'].notna().all())
            frames['train'].loc[0,'target']='stayed_vus'
            with self.assertRaisesRegex(ValueError,'nonpathogenic'):dr.train(frames,{},output_dir=tmp)

    def test_cache_invalidation_and_unmatched_rows(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);folder,labeldir=fixtures.FeatureValidationTests().make_fixture(root)
            manifest=pd.read_csv(labeldir/'label_manifest.csv',dtype='string')
            for i,row in enumerate(manifest.itertuples(index=False)):
                filename=f'vus_labels_{row.feature_snapshot}_to_{row.outcome_snapshot}.csv.gz'
                (labeldir/filename).write_bytes((labeldir/row.file).read_bytes());manifest.loc[i,'file']=filename
            manifest.to_csv(labeldir/'label_manifest.csv',index=False)
            n=json.loads((ROOT/'notebooks/03_feature_engineering.ipynb').read_text(encoding='utf-8'));ns={'display':lambda x:None}
            with contextlib.redirect_stdout(io.StringIO()):
                for c in n['cells']:
                    if c['cell_type']=='code':exec(''.join(c['source']).replace('PROJECT_ROOT = Path(r"D:\\GeneVISTA")',f'PROJECT_ROOT = Path({str(root)!r})'),ns)
            raw=root/'data/raw/clinvar';raw.mkdir(parents=True)
            for snapshot in ['2022-01','2023-01','2024-01','2025-01','2026-01']:
                pd.DataFrame({'VariationID':['1','2'],'GeneSymbol':['A','B'],'PhenotypeIDS':['MONDO:0000001']*2,'PhenotypeList':['Disease']*2,'Assembly':['GRCh38']*2}).to_csv(raw/f'variant_summary_{snapshot}.txt.gz',sep='\t',index=False)
            pd.DataFrame({'mondo_id':['MONDO:0000001'],'disease_name':['Disease']}).to_csv(folder/'clingen_gene_evidence_mapped.csv.gz',index=False)
            with patch.object(dr,'ROOT',root),patch.object(dr,'CACHE',folder/'disease_ranking_v1'),patch.object(dr,'OUT',root/'models/test'),contextlib.redirect_stdout(io.StringIO()):
                first,_=dr.prepare();second,_=dr.prepare()
                pd.testing.assert_frame_equal(first['train'],second['train'])
                f=folder/'model_features/vus_features_2022-01.csv.gz';data=pd.read_csv(f,dtype='string',keep_default_na=False);data.loc[0,'VariationID']='999';data.to_csv(f,index=False)
                with self.assertRaisesRegex(ValueError,'Unmatched'):dr.prepare()
