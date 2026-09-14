"""Conservative TP53 evidence matching. This module does not classify patients.

Published preliminary criteria stay separate from current-rule candidates.
Neither is combined with GeneVISTA model probabilities or clinical labels.
"""
import math
import re

TRANSCRIPT = 'NM_000546.6'
RULE_TRANSCRIPT = 'NM_000546.5'
RULE_VERSION = '2.4'
ARTICLE = 'https://doi.org/10.1186/s13073-025-01536-3'
CSPEC = 'https://cspec.genome.network/cspec/ui/svi/svi/GN009'
AA = dict(zip(
    'Ala Arg Asn Asp Cys Gln Glu Gly His Ile Leu Lys Met Phe Pro Ser Thr Trp Tyr Val'.split(),
    'ARNDCQEGHILKMFPSTWYV'))
NAME = re.compile(r'(NM_000546\.\d+)\(TP53\):(c\.([1-9]\d*)([ACGT])>([ACGT])) \(p\.([A-Z][a-z]{2})([1-9]\d*)([A-Z][a-z]{2})\)')


def parse_variant(name):
    match = NAME.fullmatch(name or '')
    if not match:
        return None
    transcript, change, position, ref, alt, aa_ref, aa_position, aa_alt = match.groups()
    if aa_ref not in AA or aa_alt not in AA or ref == alt or aa_ref == aa_alt:
        return None
    return {'transcript': transcript, 'coding_change': change, 'position': int(position),
            'ref': ref, 'alt': alt,
            'protein_change': f'{AA[aa_ref]}{aa_position}{AA[aa_alt]}'}


def select_functional(rows, coding_change):
    """Choose allele-specific rows before protein-level rows; retain provenance.

    Protein-equivalent alleles can have different Kotler/Funk assay results.
    Identical duplicate rows are one evidence item, not independent experiments.
    """
    substitution = re.fullmatch(r'c\.[1-9]\d*([ACGT]>[ACGT])', coding_change)
    if not substitution:
        return {'status': 'unsupported_allele', 'rows': []}
    exact = [r for r in rows if r['sequence_selector'] in (coding_change, substitution[1])]
    generic = [r for r in rows if r['sequence_selector'] in ('NA', 'AASub (no Kato data)')]
    specific = [r for r in rows if r not in generic]
    candidates = exact or (generic if not specific else [])
    if not candidates:
        return {'status': 'allele_not_covered' if rows else 'not_available', 'rows': []}
    signatures = {(r['published_code'], tuple(sorted(r['assays'].items()))) for r in candidates}
    if len(signatures) != 1:
        return {'status': 'ambiguous', 'rows': candidates}
    return {'status': 'matched', 'match_basis': 'DNA substitution and protein' if exact else 'protein-level assay',
            'rows': candidates, 'published_code': candidates[0]['published_code'],
            'assays': candidates[0]['assays']}


def computational_candidate(row, *, pvs1_applied=False):
    """Limited missense PP3/BP4 candidate per GN009 v2.4, never an applied code.

    All scores here come from the exact allele in published Table S2.
    Multiple computational predictors contribute to ONE criterion.
    """
    result = {'code': None, 'rule_version': RULE_VERSION, 'applied': False}
    if pvs1_applied:
        return dict(result, reason='PVS1 requires review of computational evidence dependencies.')
    splice, bayes, align = row.get('spliceai'), row.get('bayesdel'), row.get('align_gvgd')
    if not isinstance(splice, (int, float)) or not math.isfinite(splice) or not 0 <= splice <= 1:
        return dict(result, reason='Splicing prediction is missing or invalid.')
    if splice >= 0.2:
        return dict(result, code='PP3', basis='splicing', reason='Predicted splice impact requires RNA/splice review.')
    if not isinstance(bayes, (int, float)) or not math.isfinite(bayes) or align not in {'C0','C15','C25','C35','C45','C55','C65'}:
        return dict(result, reason='Missense predictors are missing or invalid.')
    code = None
    if bayes >= 0.16:
        if align == 'C65':
            code = 'PP3_Moderate'
        elif align in {'C25', 'C35', 'C45', 'C55'}:
            code = 'PP3'
    elif align != 'C65':
        code = 'BP4_Moderate' if bayes <= -0.008 else 'BP4'
    return dict(result, code=code, basis='missense', reason='Candidate for expert review.' if code else 'No computational criterion met.')


def assess_variant(record, bundle, *, pvs1_applied=False, pvs1_for_splicing=False, observed_splicing=None):
    if not record or record.get('gene_symbol') != 'TP53':
        return {'status': 'out_of_scope', 'reason': 'This pilot covers TP53 missense substitutions only.'}
    variant = parse_variant(record.get('name'))
    if not variant:
        return {'status': 'unsupported_variant', 'reason': 'An exact TP53 coding SNV and missense protein description are required.'}
    if variant['transcript'] != TRANSCRIPT:
        return {'status': 'transcript_mismatch', 'reason': f'The source tables use {TRANSCRIPT}; this transcript has not been matched.'}
    metadata = bundle.get('metadata', {})
    proof = metadata.get('transcript_verification', {})
    if (metadata.get('schema_version') != 1 or metadata.get('table_transcript') != TRANSCRIPT
            or metadata.get('rule_version') != RULE_VERSION or proof.get('cds_identical') is not True
            or proof.get('accessions') != [RULE_TRANSCRIPT, TRANSCRIPT]):
        return {'status': 'incompatible_source', 'reason': 'Evidence package or transcript verification needs rebuilding.'}
    row = bundle.get('computational', {}).get(variant['coding_change'])
    if row is None:
        return {'status': 'allele_not_covered', 'reason': 'This DNA substitution is not in the published computational table.'}
    if row['protein_change'] != variant['protein_change']:
        return {'status': 'mapping_conflict', 'reason': 'The recorded protein consequence differs from the source table; no evidence was assigned.'}
    computational = computational_candidate(row, pvs1_applied=pvs1_applied or pvs1_for_splicing)
    functional = select_functional(bundle.get('functional', {}).get(variant['protein_change'], []), variant['coding_change'])
    holds = []
    if row['spliceai'] >= 0.2:
        holds.append('Predicted splicing effect: protein-assay criteria are withheld pending splice review.')
    if pvs1_for_splicing or pvs1_applied:
        holds.append('PVS1 context: review functional-criterion dependencies before use.')
    if observed_splicing is True:
        holds.append('Observed splicing aberration: review RNA evidence before protein-assay criteria.')
    functional['criterion_status'] = 'withheld' if holds else 'review_required'
    functional['hold_reasons'] = holds
    functional['applied'] = False
    # Source functional codes are from the publication's v2 tables. They have
    # NOT been recalculated against the updated v2.4 functional flowchart.
    functional['rule_review_required'] = True
    direction_conflict = False
    fc = functional.get('published_code', '')
    cc = computational.get('code') or ''
    if functional['status'] == 'matched' and not holds:
        direction_conflict = ((cc.startswith('PP') and fc.startswith('BS'))
                              or (cc.startswith('BP') and fc.startswith('PS')))
    return {
        'status': 'available', 'variant': variant,
        'condition': {'name': 'Li-Fraumeni syndrome', 'id': 'MONDO:0018875', 'inheritance': 'Autosomal dominant'},
        'computational': dict(row, candidate=computational), 'functional': functional,
        'conflicting_evidence_directions': direction_conflict,
        'classification': None, 'patient_disease_probability': None,
        'review_required': [
            'Review published functional criteria against the current TP53 specification.',
            'Assess RNA/splicing evidence and confirm constitutional variant origin.',
            'Review population, phenotype, segregation and other applicable evidence before classification.',
        ],
        'sources': metadata['sources'], 'rule_version': RULE_VERSION,
        'limitations': 'Research evidence preview. No complete ACMG/AMP classification, diagnostic validation, or disease-risk estimate. Criteria are not added to model probabilities.',
    }
