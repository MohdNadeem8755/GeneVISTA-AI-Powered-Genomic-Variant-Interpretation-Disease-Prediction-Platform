import gzip
import json
from pathlib import Path
import tempfile
import unittest
from collections import Counter
import pandas as pd
from src.cleaning_validation import exact_class, export_exact_classes, parsed_dates, source_date
ROOT=Path(__file__).resolve().parents[1]

class CleaningValidationTests(unittest.TestCase):
    def test_exact_targets_keep_likely_and_mixed_separate(self):
        for label in ('Pathogenic','Likely pathogenic','Benign','Likely benign','Uncertain significance'):
            self.assertEqual(exact_class(' '+label+' '),label)
        for label in ('Pathogenic/Likely pathogenic','Conflicting classifications of pathogenicity','-','',pd.NA):
            self.assertTrue(pd.isna(exact_class(label)))

    def test_mixed_dates_and_invalid_dates(self):
        result=parsed_dates(pd.Series(['2024-03-14T16:00:00.000Z','2020-01-05',pd.NA],dtype='string'))
        self.assertEqual(result.dt.strftime('%Y-%m-%d').iloc[0],'2024-03-14')
        self.assertTrue(pd.isna(result.iloc[2]))
        with self.assertRaises(ValueError): parsed_dates(pd.Series(['bad date'],dtype='string'))

    def test_export_retains_unresolved_and_reruns(self):
        with tempfile.TemporaryDirectory() as folder:
            p=Path(folder); source=p/'variants.csv.gz'; dest=p/'exact.csv.gz'
            pd.DataFrame({'VariationID':['1','2','3'],'snapshot':['2026-09']*3,
                'ClinicalSignificance':['Likely pathogenic','Benign/Likely benign','Benign']}).to_csv(source,index=False)
            for _ in range(2):
                report=export_exact_classes(source,dest,chunksize=2)
                self.assertEqual((report['records'],report['five_class_eligible'],report['unresolved']),(3,2,1))
                rows=pd.read_csv(dest,keep_default_na=False)
                self.assertEqual(rows['ClinicalSignificance'].iloc[1],'Benign/Likely benign')
                self.assertEqual(rows['exact_clinical_class'].iloc[1],'')

    def test_notebook_cleaner_chunking_and_rerun(self):
        notebook=json.loads((ROOT/'notebooks/02_data_cleaning.ipynb').read_text(encoding='utf-8'))
        setup=next(''.join(c['source']) for c in notebook['cells'] if 'def simplify_clinical_state' in ''.join(c.get('source',[])))
        clean=next(''.join(c['source']) for c in notebook['cells'] if 'def clean_snapshot' in ''.join(c.get('source',[])))
        with tempfile.TemporaryDirectory() as folder:
            p=Path(folder)
            columns=['VariationID','GeneSymbol','ClinicalSignificance','ReviewStatus','Assembly','Name']
            rows=[['1','GENE','Likely pathogenic','reviewed','GRCh38','NA'],
                  ['2','GENE','Benign','reviewed','GRCh37','x'],
                  ['3','-','Pathogenic','reviewed','GRCh38','x'],
                  ['4','GENE','Likely benign','reviewed','GRCh38','x']]
            source=p/'variant_summary_2022-01.txt.gz'
            pd.DataFrame(rows,columns=columns).to_csv(source,sep='\t',index=False)
            ns={'pd':pd,'Counter':Counter,'CLINVAR_FOLDER':p,'PROJECT_ROOT':p,'CORE_COLUMNS':columns}
            exec(setup,ns)
            # Keep fixture output separate from source; test the exact notebook function.
            exec(clean[:clean.index('\ncleaning_2022 =')].replace('chunksize=100_000','chunksize=2'),ns)
            for _ in range(2):
                report=ns['clean_snapshot'](source)
                self.assertEqual(report['saved_rows'],2)
                self.assertEqual(report['other_assembly_rows'],1)
                self.assertEqual(report['missing_field_rows'],1)
                saved=pd.read_csv(ns['processed_folder']/'clinvar_clean_2022-01.csv.gz',keep_default_na=False)
                self.assertEqual(saved['Name'].iloc[0],'NA')
                self.assertEqual(saved['ClinicalSignificance'].tolist(),['Likely pathogenic','Likely benign'])

    def test_source_date_from_metadata(self):
        with tempfile.TemporaryDirectory() as folder:
            p=Path(folder)/'clingen.csv'; p.write_text('"TITLE"\n"FILE CREATED: 2026-09-08",x\n')
            self.assertEqual(source_date(p),'2026-09-08')

class TimelinePersistenceTests(unittest.TestCase):
    def test_missing_and_stale_timeline_are_rebuilt(self):
        notebook=json.loads((ROOT/'notebooks/02_data_cleaning.ipynb').read_text(encoding='utf-8'))
        code=next(''.join(c['source']) for c in notebook['cells'] if 'temporary_timeline =' in ''.join(c.get('source',[])))
        with tempfile.TemporaryDirectory() as folder:
            p=Path(folder)
            timeline=pd.DataFrame({'2022-01':['vus','benign'],'2023-01':['pathogenic','benign']},index=pd.Index([1,2],name='VariationID'))
            ns={'processed_folder':p,'pd':pd,'timeline':timeline,'Counter':Counter,
                'snapshots':['2022-01','2023-01'],'state_names':['vus','benign','pathogenic'],'display':lambda x:None}
            for stale in (False,True):
                if stale: (p/'clinvar_state_timelines.csv.gz').write_bytes(b'old corrupt output')
                exec(code,ns)
                saved=pd.read_csv(p/'clinvar_state_timelines.csv.gz',index_col='VariationID')
                pd.testing.assert_frame_equal(saved,timeline)
                self.assertFalse((p/'clinvar_state_timelines.partial.csv.gz').exists())
