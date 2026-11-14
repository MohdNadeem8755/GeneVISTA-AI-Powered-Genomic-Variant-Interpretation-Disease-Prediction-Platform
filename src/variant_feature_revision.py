"""Opt-in v2 inputs for NEW fits; preserve legacy saved-model preprocessing."""
import re
import pandas as pd
from src.train_annotation_candidates import INDEPENDENT, SEQUENCE, matrix
from src.variant_features import parse_change

PROFILE='variant_identity_v2'
CATEGORICAL=('gene_identity','genomic_base_change','transcript_base_change',
             'amino_change','amino_reference','amino_alternate')


def identity_features(gene, name, reference, alternate):
    change=parse_change(name)
    nt=change['nucleotide']; aa=change['protein']
    missing='__MISSING__'
    genes=sorted({g.strip() for g in re.split(r'[;,|]',str(gene)) if g.strip() not in ('','-')})
    values={'gene_identity':'|'.join(genes) if genes else missing,
            'genomic_base_change':f'{reference}>{alternate}' if reference in ('A','C','G','T') and alternate in ('A','C','G','T') and reference!=alternate else missing,
            'transcript_base_change':f"{nt['reference']}>{nt['alternate']}" if nt else missing,
            'amino_change':f"{aa['reference']}>{aa['alternate']}" if aa else missing,
            'amino_reference':aa['reference'] if aa else missing,
            'amino_alternate':aa['alternate'] if aa else missing}
    # Inspect HGVS change expressions only, never gene/transcript names.
    expressions=re.findall(r'\b[cgmnrp]\.\(?([^\s,;)]+)',name)
    sequence_text=' '.join(expressions).lower()
    for token in ('del','dup','ins','fs','inv'):
        values['event='+token]=float(token in sequence_text)
    return values


def revised_base_matrix(rows):
    """Rows: label, annotation JSON, gene, HGVS name, type, genomic ref, alt.

    Labels are passed through only for the legacy matrix's positional interface;
    neither the legacy matrix nor this transform uses them as feature inputs.
    """
    columns=tuple(INDEPENDENT)+SEQUENCE
    base=pd.DataFrame(matrix([r[:5] for r in rows],columns),columns=columns)
    revised=pd.DataFrame([identity_features(r[2],r[3],r[5],r[6]) for r in rows])
    for column in revised:
        base[column]=revised[column]
    return base
