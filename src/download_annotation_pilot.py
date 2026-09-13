"""Download a bounded public annotation pilot; preserves source/version and assembly."""
import csv,gzip,json,urllib.request,urllib.parse,hashlib,datetime
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data/raw/external_annotations_pilot';OUT.mkdir(parents=True,exist_ok=True)
rows=[];seen=set()
with gzip.open(ROOT/'data/raw/clinvar/variant_summary_2026-09.txt.gz','rt',encoding='utf-8-sig') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  if r['Assembly']!='GRCh37' or r['OriginSimple'] not in ('germline','unknown'):continue
  if r['ClinicalSignificance'] not in ['Pathogenic','Likely pathogenic','Uncertain significance','Likely benign','Benign']:continue
  ref,alt=r['ReferenceAlleleVCF'],r['AlternateAlleleVCF'];pos=r['PositionVCF'];chrom=r['Chromosome']
  if len(ref)!=1 or len(alt)!=1 or ref not in 'ACGT' or alt not in 'ACGT' or ref==alt or not pos.isdigit() or int(pos)<=0:continue
  if chrom not in [str(i) for i in range(1,23)]+['X','Y']:continue
  key=f'chr{chrom}:g.{pos}{ref}>{alt}'
  if key in seen:continue
  # Deterministic small pilot; not an evaluation cohort.
  if int.from_bytes(hashlib.sha256(key.encode()).digest()[:4],'big')%101!=0:continue
  seen.add(key);rows.append({'variation_id':r['VariationID'],'query':key,'assembly':'GRCh37','gene':r['GeneSymbol']})
  if len(rows)>=100:break
fields='dbnsfp.revel,dbnsfp.cadd,dbnsfp.sift,dbnsfp.polyphen2,dbnsfp.gerp++,dbnsfp.phylop,gnomad_exome.af,gnomad_genome.af,cadd.phred,cadd.consequence'
payload=urllib.parse.urlencode({'ids':','.join(r['query'] for r in rows),'fields':fields}).encode()
req=urllib.request.Request('https://myvariant.info/v1/variant',data=payload,headers={'User-Agent':'GeneVISTA-research-pilot/1.0','Content-Type':'application/x-www-form-urlencoded'})
with urllib.request.urlopen(req,timeout=90) as response:data=json.load(response)
(OUT/'queries.json').write_text(json.dumps(rows,indent=2),encoding='utf-8')
(OUT/'annotations.json').write_text(json.dumps(data,indent=2),encoding='utf-8')
try:
 with urllib.request.urlopen('https://myvariant.info/v1/metadata',timeout=45) as response:metadata=json.load(response)
 (OUT/'source_metadata.json').write_text(json.dumps(metadata,indent=2),encoding='utf-8')
except Exception as e:print('Metadata unavailable',type(e).__name__)
summary={'retrieved_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'endpoint':'https://myvariant.info/v1/variant','assembly':'GRCh37/hg19','requested_fields':fields,'requested':len(rows),'returned':len(data),'matched':sum(not r.get('notfound',False) for r in data),'with_dbnsfp':sum(bool(r.get('dbnsfp')) for r in data),'with_population_frequency':sum(bool(r.get('gnomad_exome') or r.get('gnomad_genome')) for r in data),'with_cadd':sum(bool(r.get('cadd')) for r in data),'purpose':'Coverage pilot only. These scores are not independent truth labels. Upstream predictor training overlap and licensing must be audited before model evaluation or redistribution.'}
(OUT/'coverage.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
print(json.dumps(summary,indent=2))
