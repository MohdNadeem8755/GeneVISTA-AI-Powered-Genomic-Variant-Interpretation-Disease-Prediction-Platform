import gzip
import json
from pathlib import Path
import unittest
from src.biological_features import enrich, gene_context, load_development

ROOT=Path(__file__).resolve().parents[1]


class BiologicalFeatureTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.brca=json.load(gzip.open(ROOT/'data/processed/brca1_mavedb_v1.json.gz','rt'))
        cls.tp53=json.load(gzip.open(ROOT/'data/processed/tp53_evidence_v1.json.gz','rt'))

    def test_zero_frequency_and_missing_are_distinct(self):
        values=enrich('UNKNOWN','',{'gnomAD4.1_joint_AF':0},{},self.brca,self.tp53)
        self.assertEqual(values['gnomAD4.1_joint_AF'],0)
        self.assertEqual(values['gnomAD4.1_joint_AF__missing'],0)
        self.assertIsNone(values['gnomAD4.1_joint_POPMAX_AF'])
        self.assertEqual(values['gnomAD4.1_joint_POPMAX_AF__missing'],1)

    def test_negative_assay_preserved_and_no_clinical_codes(self):
        values=enrich('BRCA1','NM_007294.4(BRCA1):c.190T>G (p.Cys64Gly)',{}, {},self.brca,self.tp53)
        self.assertAlmostEqual(values['brca1_function_score'],-1.95113526060798)
        self.assertNotIn('PS3',json.dumps(values))
        self.assertFalse(any('label' in key or 'criterion' in key for key in values))

    def test_wrong_transcript_rejected(self):
        values=enrich('BRCA1','NM_007294.99(BRCA1):c.190T>G (p.Cys64Gly)',{}, {},self.brca,self.tp53)
        self.assertIsNone(values['brca1_function_score'])

    def test_gene_context_deduplicates_relations(self):
        row={'mondo_id':'MONDO:1','mode_of_inheritance':'AD','classification':'Definitive'}
        values=gene_context('GENE',{'GENE':[row,row]})
        self.assertEqual(values['clingen_disease_count'],1)
        self.assertEqual(values['clingen_definitive_count'],1)
        self.assertEqual(values['inheritance_AD'],1)
        self.assertEqual(gene_context('GENE;OTHER',{'GENE':[row]}),{'clingen_context_available':0})

    def test_test_labels_cannot_be_loaded(self):
        with self.assertRaises(ValueError):
            load_development('unused','unused','test')


if __name__=='__main__':
    unittest.main()
