import copy
import gzip
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from fastapi.testclient import TestClient
from backend.app.main import app
from backend.app.tp53_evidence import evidence
from src.tp53_evidence import assess_variant, computational_candidate, parse_variant, select_functional
from src.build_tp53_evidence import verify_transcripts

ROOT = Path(__file__).resolve().parents[1]


class Tp53EvidenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with gzip.open(ROOT/'data/processed/tp53_evidence_v1.json.gz', 'rt', encoding='utf-8') as handle:
            cls.package = json.load(handle)
        cls.record = {'name': 'NM_000546.6(TP53):c.742C>T (p.Arg248Trp)', 'gene_symbol': 'TP53'}

    def test_exact_transcript_and_protein_are_required(self):
        self.assertEqual(assess_variant(self.record, self.package)['status'], 'available')
        for name, expected in [
            (self.record['name'].replace('.6(', '.5('), 'transcript_mismatch'),
            (self.record['name'].replace('Arg248Trp', 'Arg248Gln'), 'mapping_conflict'),
            ('NM_000546.6(TP53):c.742_743del (p.Arg248fs)', 'unsupported_variant'),
        ]:
            self.assertEqual(assess_variant(dict(self.record, name=name), self.package)['status'], expected)
        self.assertIsNone(parse_variant(self.record['name']+' extra text'))
        self.assertEqual(assess_variant(dict(self.record, gene_symbol='BRCA1'), self.package)['status'], 'out_of_scope')

    def test_allele_specific_functional_results_are_not_mixed(self):
        rows = self.package['functional']['G117R']
        ga = select_functional(rows, 'c.349G>A')
        gc = select_functional(rows, 'c.349G>C')
        self.assertEqual(ga['published_code'], 'No evidence')
        self.assertEqual(gc['published_code'], 'BS3_Supporting')
        self.assertEqual(select_functional(rows, 'c.349G>T')['status'], 'allele_not_covered')
        exact = select_functional(self.package['functional']['G187D'], 'c.560G>A')
        self.assertEqual(len(exact['rows']), 1)
        self.assertEqual(exact['published_code'], 'BS3_Supporting')

    def test_duplicate_assays_are_one_criterion_and_conflicts_abstain(self):
        rows = self.package['functional']['F212L']
        result = select_functional(rows, 'c.634T>C')
        self.assertEqual(result['status'], 'matched')
        self.assertEqual(len(result['rows']), 2)  # provenance retained, one code
        contradictory = copy.deepcopy(rows)
        contradictory[0]['published_code'] = 'PS3'
        self.assertEqual(select_functional(contradictory, 'c.634T>C')['status'], 'ambiguous')

    def test_computational_thresholds_and_dependence(self):
        row = {'spliceai': 0.19, 'bayesdel': 0.16, 'align_gvgd': 'C65'}
        self.assertEqual(computational_candidate(row)['code'], 'PP3_Moderate')
        self.assertEqual(computational_candidate(dict(row, align_gvgd='C25'))['code'], 'PP3')
        self.assertIsNone(computational_candidate(dict(row, align_gvgd='C15'))['code'])
        self.assertEqual(computational_candidate(dict(row, bayesdel=-0.008, align_gvgd='C0'))['code'], 'BP4_Moderate')
        self.assertEqual(computational_candidate(dict(row, bayesdel=-0.0079, align_gvgd='C0'))['code'], 'BP4')
        self.assertIsNone(computational_candidate(dict(row, bayesdel=-0.1))['code'])
        splice = computational_candidate(dict(row, spliceai=0.2))
        self.assertEqual((splice['code'], splice['basis']), ('PP3', 'splicing'))
        self.assertIsNone(computational_candidate(row, pvs1_applied=True)['code'])
        for missing in [None, float('nan'), float('inf'), -0.1, 1.1]:
            self.assertIsNone(computational_candidate(dict(row, spliceai=missing))['code'])

    def test_splicing_and_pvs1_withhold_functional_criteria(self):
        record = dict(self.record, name='NM_000546.6(TP53):c.105G>T (p.Leu35Phe)')
        result = assess_variant(record, self.package)
        self.assertEqual(result['functional']['criterion_status'], 'withheld')
        self.assertEqual(result['computational']['candidate']['basis'], 'splicing')
        for context in [dict(observed_splicing=True), dict(pvs1_applied=True), dict(pvs1_for_splicing=True)]:
            result = assess_variant(self.record, self.package, **context)
            self.assertEqual(result['functional']['criterion_status'], 'withheld')
        self.assertEqual(assess_variant(self.record, self.package)['functional']['criterion_status'], 'review_required')

    def test_evidence_is_not_diagnosis_accuracy_or_added_to_model(self):
        first = assess_variant(self.record, self.package)
        with_label = assess_variant(dict(self.record, clinical_significance='Benign', prediction={'top_probability':0.999}), self.package)
        self.assertEqual(first, with_label)
        self.assertIsNone(first['classification'])
        self.assertIsNone(first['patient_disease_probability'])
        self.assertFalse(first['functional']['applied'])
        self.assertFalse(first['computational']['candidate']['applied'])
        self.assertNotIn('accuracy', first)

    def test_source_version_and_cds_proof_gate(self):
        for key, value in [('rule_version', '99'), ('schema_version', 2), ('table_transcript', 'NM_000546.5')]:
            package = dict(self.package, metadata=dict(self.package['metadata'], **{key: value}))
            self.assertEqual(assess_variant(self.record, package)['status'], 'incompatible_source')
        package = dict(self.package, metadata=dict(self.package['metadata'], transcript_verification={}))
        self.assertEqual(assess_variant(self.record, package)['status'], 'incompatible_source')

    def test_cds_comparison_rejects_sequence_differences(self):
        cds = 'ATG'+'GCT'*392+'TAA'
        def gb(accession, sequence):
            return f'VERSION     {accession}\n     CDS             1..1182\nORIGIN\n        1 {sequence.lower()}\n//\n'
        text = gb('NM_000546.5', cds)+gb('NM_000546.6', cds)
        _, proof = verify_transcripts(text)
        self.assertTrue(proof['cds_identical'])
        changed = gb('NM_000546.5', cds)+gb('NM_000546.6', cds[:5]+'A'+cds[6:])
        with self.assertRaises(ValueError):
            verify_transcripts(changed)

    def test_imported_counts_exclusions_and_missing_file(self):
        counts = self.package['metadata']['counts']
        self.assertEqual(counts['coding_alleles'], 2569)
        self.assertEqual(counts['excluded_functional_rows'], 6)
        self.assertTrue(all(r['protein_change'].startswith('R72') for r in self.package['metadata']['excluded_functional_rows']))
        self.assertNotIn('R72A', self.package['functional'])
        with patch.dict(os.environ, {'GENEVISTA_TP53_EVIDENCE':str(ROOT/'data/does-not-exist.json.gz')}):
            self.assertEqual(evidence(self.record)['status'], 'not_configured')

    def test_api_and_authentication(self):
        with patch.dict(os.environ, {'GENEVISTA_AUTH_ENABLED':'false'}), TestClient(app) as client:
            result = client.get('/api/variants/12347/tp53-evidence')
            self.assertEqual(result.status_code, 200)
            self.assertEqual(result.json()['status'], 'available')
            self.assertEqual(client.get('/api/variants/0/tp53-evidence').status_code, 422)
            self.assertEqual(client.get('/api/variants/9223372036854775807/tp53-evidence').status_code, 404)
        with patch.dict(os.environ, {'GENEVISTA_AUTH_ENABLED':'true'}), TestClient(app) as client:
            self.assertEqual(client.get('/api/variants/12347/tp53-evidence').status_code, 401)


if __name__ == '__main__':
    unittest.main()
