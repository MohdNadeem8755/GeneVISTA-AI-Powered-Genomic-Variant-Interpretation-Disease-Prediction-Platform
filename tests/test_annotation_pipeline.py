import csv
import gzip
import hashlib
import json
from pathlib import Path
import tempfile
import unittest

from src.annotation_pipeline import (CLASSES, FEATURES, assign_split, bucket,
    check_md5, compare, connect, merge, meta, number, prepare, snv_key)


class AnnotationPipelineTests(unittest.TestCase):
    def setUp(self):
        root = Path(__file__).resolve().parents[1] / '.runtime'
        root.mkdir(exist_ok=True)
        self.temp = tempfile.TemporaryDirectory(dir=root)
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.source = self.root / 'clinvar.tsv.gz'
        self.db = self.root / 'cohort.sqlite'

    def row(self, pos, label='Pathogenic', **changes):
        row = dict(Assembly='GRCh38', Chromosome='1', PositionVCF=str(pos),
            ReferenceAlleleVCF='A', AlternateAlleleVCF='G', OriginSimple='germline',
            VariationID=str(pos), GeneSymbol='GENE1', Name=f'NM_1:c.{pos}A>G',
            Type='single nucleotide variant', ClinicalSignificance=label,
            ReviewStatus='criteria provided, single submitter')
        row.update(changes)
        return row

    def cohort(self, rows):
        with gzip.open(self.source, 'wt', newline='', encoding='utf-8') as stream:
            writer=csv.DictWriter(stream,fieldnames=list(rows[0]),delimiter='\t')
            writer.writeheader(); writer.writerows(rows)
        return prepare(self.source,self.db,1)

    def annotations(self, rows):
        file=self.root/'dbNSFP5.4a_grch38.gz'
        header=['#chr','pos(1-based)','ref','alt',*FEATURES,'clinvar_clnsig']
        with gzip.open(file,'wt',encoding='utf-8',newline='') as stream:
            writer=csv.writer(stream,delimiter='\t'); writer.writerow(header)
            writer.writerows(rows)
        checksum=self.root/(file.name+'.md5')
        checksum.write_text(hashlib.md5(file.read_bytes()).hexdigest()+'  '+file.name)
        return file,checksum

    def test_coordinates_exact_labels_and_conflicts(self):
        rows=[self.row(100),self.row(100),self.row(200),self.row(200,'Benign'),
            self.row(300,'Pathogenic/Likely pathogenic'),self.row(400,Assembly='GRCh37'),
            self.row(500,AlternateAlleleVCF='GG'),self.row(600,OriginSimple='somatic'),
            self.row(700,VariationID='same'),self.row(800,VariationID='same')]
        audit=self.cohort(rows)
        self.assertEqual(audit['conflicting_keys'],3)
        self.assertEqual(audit['duplicate_source_rows_collapsed'],2)
        with connect(self.db) as con:
            self.assertEqual(con.execute("SELECT COUNT(*) FROM variants WHERE split<>'excluded'").fetchone()[0],1)
        self.assertIsNone(snv_key(self.row(100,Assembly='GRCh37')))
        self.assertIsNone(snv_key(self.row(100,ReferenceAlleleVCF='N')))
        with self.assertRaises(ValueError): prepare(self.source,self.db,1)

    def test_locus_groups_and_legacy_test_exclusion(self):
        self.cohort([self.row(p,CLASSES[p%5],AlternateAlleleVCF=a,VariationID=f'{p}{a}')
                     for p in range(100,200) for a in ('G','T')])
        with connect(self.db) as con:
            groups=con.execute('SELECT chrom,pos,COUNT(DISTINCT split) FROM variants GROUP BY chrom,pos').fetchall()
        self.assertTrue(all(count==1 for _,_,count in groups))
        for pos in range(100,1000):
            if bucket(f'genevista-v1|1|{pos}|{pos}',10)==0:
                self.assertNotEqual(assign_split('1',pos),'test')

    def test_merge_exact_alleles_missing_scores_and_allowlist(self):
        self.cohort([self.row(100),self.row(101)])
        rows=[['1',100,'A','G',*(['.']*len(FEATURES)),'Pathogenic'],
              ['1',100,'A','G',*(['0.2;0.7']*len(FEATURES)),'Benign'],
              ['1',101,'A','T',*(['0.9']*len(FEATURES)),'Pathogenic']]
        file,checksum=self.annotations(rows)
        report=merge(file,checksum,self.db)
        self.assertEqual(report['matched_variants'],1)
        with connect(self.db) as con:
            scores,n=con.execute('SELECT scores,source_rows FROM annotations').fetchone()
            self.assertEqual(meta(con,'annotation_status'),'ready')
        self.assertEqual(n,2)
        self.assertEqual(set(json.loads(scores)),set(FEATURES))
        self.assertEqual(json.loads(scores)['REVEL_score'],0.7)
        self.assertIsNone(number('.;NaN;inf'))
        self.assertEqual(number('0'),0)
        with self.assertRaises(ValueError): merge(file,checksum,self.db)

    def test_incomplete_download_and_no_annotations_fail_closed(self):
        self.cohort([self.row(100)])
        file,checksum=self.annotations([['1',100,'A','G',*(['.']*len(FEATURES)),'VUS']])
        checksum.write_text('0'*32)
        with self.assertRaisesRegex(ValueError,'MD5 mismatch'): merge(file,checksum,self.db)
        with self.assertRaisesRegex(ValueError,'Complete and verify'): compare(self.db,self.root/'comparison.json')
        with connect(self.db) as con:
            self.assertEqual(con.execute('SELECT COUNT(*) FROM annotations').fetchone()[0],0)

    def test_comparison_never_uses_test_labels(self):
        rows=[self.row(p,CLASSES[p%5]) for p in range(1000,1500)]
        self.cohort(rows)
        file,checksum=self.annotations([['1',p,'A','G',*([str((p%7)/7)]*len(FEATURES)),'not-an-input'] for p in range(1000,1500)])
        merge(file,checksum,self.db)
        with connect(self.db) as con:
            con.execute("UPDATE variants SET label='POISONED_TEST_LABEL' WHERE split IN ('test','calibration')")
        result=compare(self.db,self.root/'comparison.json',100)
        self.assertEqual(result['test_rows_read'],0)
        self.assertEqual(result['calibration_rows_read'],0)
        self.assertFalse(result['production_model_changed'])
        self.assertEqual(len(result['candidates']),2)

if __name__=='__main__': unittest.main()
