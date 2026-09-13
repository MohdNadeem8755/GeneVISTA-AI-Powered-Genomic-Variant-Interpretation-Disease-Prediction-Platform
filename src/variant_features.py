"""Sequence-derived features shared by training and inference; no clinical labels."""
import math
import re

CLASSES = ['Pathogenic','Likely pathogenic','Uncertain significance','Likely benign','Benign']

def parse_change(name):
    nucleotide = sorted(set(re.findall(r'(?:c|g|m|n)\.([*\-]?\d+(?:[+-]\d+)?)([ACGT])>([ACGT])',name)))
    proteins = sorted(set(re.findall(r'p\.\(?([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|Ter|\*|=)',name)))
    change = {'nucleotide':None,'protein':None}
    if len(nucleotide)==1:
        position,reference,alternate=nucleotide[0]
        change['nucleotide']={'reference':reference,'alternate':alternate,'position':position,'label':reference+' → '+alternate}
    if len(proteins)==1:
        reference,position,alternate=proteins[0]
        change['protein']={'reference':reference,'alternate':alternate,'position':int(position),'label':reference+position+alternate}
    return change

def variant_features(row):
    name=str(row.get('name',''))
    genes=sorted(set(g.strip() for g in re.split(r'[;,|]',str(row.get('gene_symbol',''))) if g.strip() and g.strip()!='-'))
    features={'gene='+gene:1.0 for gene in genes}
    features['type='+str(row.get('variant_type','unknown')).lower()]=1.0
    change=parse_change(name)
    nucleotide=change['nucleotide'];protein=change['protein']
    if nucleotide:
        ref,alt=nucleotide['reference'],nucleotide['alternate']
        features['base='+ref+'>'+alt]=1.0
        features['transition']=float((ref,alt) in [('A','G'),('G','A'),('C','T'),('T','C')])
        position=nucleotide['position']
        features['intronic']=float('+' in position or '-' in position[1:])
        features['utr']=float(position.startswith(('*','-')))
        features['log_cdna_position']=math.log1p(int(re.search(r'\d+',position)[0]))/10
    if protein:
        ref,alt=protein['reference'],protein['alternate']
        features['amino='+ref+'>'+alt]=1.0
        features['amino_ref='+ref]=1.0;features['amino_alt='+alt]=1.0
        features['synonymous']=float(ref==alt or alt=='=')
        features['stop_gain']=float(alt in ('Ter','*'))
        features['log_protein_position']=math.log1p(protein['position'])/10
    lower=name.lower()
    for token in ('del','dup','ins','fs','inv'):
        features['event='+token]=float(token in lower)
    features['has_nucleotide_change']=float(nucleotide is not None)
    features['has_protein_change']=float(protein is not None)
    return features
