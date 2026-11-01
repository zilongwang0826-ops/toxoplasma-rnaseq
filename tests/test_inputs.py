"""Input integrity and reference invariants; no network or biological inference."""
import csv
import gzip
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'workflow/scripts'))
from common import read_samples
from download import fetch
from prepare_reference import prepare
from quant_qc import inspect_quant,check
from validate_inputs import validate

class Inputs(unittest.TestCase):
    def setUp(self):
        (ROOT/'results/validation').mkdir(parents=True,exist_ok=True)
        self.temp=tempfile.TemporaryDirectory(dir=ROOT/'results/validation')
        self.root=Path(self.temp.name)
        self.rows=read_samples(ROOT/'metadata/samples.tsv')
    def tearDown(self):
        self.temp.cleanup()  # Only this test's temporary synthetic fixtures.
    def sheet(self,rows):
        p=self.root/'samples.tsv'
        with p.open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
        return p
    def test_synthetic_fixture_is_deterministic(self):
        sys.path.insert(0,str(ROOT/'tests'))
        from make_smoke import generate
        generate(self.root/'a');generate(self.root/'b')
        self.assertEqual((self.root/'a/sample_1.fastq.gz').read_bytes(),(self.root/'b/sample_1.fastq.gz').read_bytes())
        with gzip.open(self.root/'a/sample_1.fastq.gz','rt') as f:
            self.assertEqual(sum(1 for _ in f),40000)

    def test_real_metadata(self):
        out=validate(ROOT/'metadata/samples.tsv')
        self.assertEqual(out['n_libraries'],12)
        self.assertEqual(out['archive_bytes'],70322230064)
    def test_swapped_condition_rejected_against_archive(self):
        self.rows[0]['condition'],self.rows[3]['condition']=self.rows[3]['condition'],self.rows[0]['condition']
        with self.assertRaises(ValueError):validate(self.sheet(self.rows), source_dir=ROOT/'metadata/source')

    def test_duplicate_biosample(self):
        self.rows[1]['biosample']=self.rows[0]['biosample']
        with self.assertRaises(ValueError): read_samples(self.sheet(self.rows))
    def test_duplicated_fastq(self):
        self.rows[1]['fastq_1']=self.rows[0]['fastq_1']
        with self.assertRaises(ValueError): read_samples(self.sheet(self.rows))
    def test_wrong_mate(self):
        self.rows[0]['fastq_1'],self.rows[0]['fastq_2']=self.rows[0]['fastq_2'],self.rows[0]['fastq_1']
        with self.assertRaises(ValueError): read_samples(self.sheet(self.rows))
    def test_duplicate_replicate(self):
        self.rows[1]['replicate']=self.rows[0]['replicate']
        with self.assertRaises(ValueError): read_samples(self.sheet(self.rows))
    def test_download_preserves_corruption(self):
        p=self.root/'input.gz';p.write_bytes(b'corrupted')
        with self.assertRaises(ValueError): fetch('https://example.org/input.gz',p,'0'*32,100)
        self.assertEqual(p.read_bytes(),b'corrupted')
    def test_existing_download_rebuilds_receipt(self):
        p=self.root/'input.gz';p.write_bytes(b'test')
        with patch('subprocess.run') as run:
            fetch('https://example.org/input.gz',p,hashlib.md5(b'test').hexdigest(),4)
            run.assert_not_called()
        self.assertTrue(Path(str(p)+'.source.json').exists())
    def test_non_https_rejected(self):
        with self.assertRaises(ValueError): fetch('ftp://example.org/file',self.root/'bad')
    def test_reference_mapping_and_decoys(self):
        for filename,text in {'cdna.gz':'>ENSMUST1.2 gene:ENSMUSG1\nACGTACGT\n',
                              'ncrna.gz':'>ENSMUST2.1\nTTTTGGGG\n',
                              'genome.gz':'>1 dna\nACGTACGTAAA\n',
                              'gtf.gz':'1\tx\ttranscript\t1\t8\t.\t+\t.\tgene_id "ENSMUSG1"; transcript_id "ENSMUST1"; gene_name "A";\n1\tx\ttranscript\t1\t8\t.\t+\t.\tgene_id "ENSMUSG2"; transcript_id "ENSMUST2";\n'}.items():
            with gzip.open(self.root/filename,'wt') as f:f.write(text)
        prepare(*(self.root/x for x in ['cdna.gz','ncrna.gz','gtf.gz']),self.root/'prepared',self.root/'genome.gz')
        data=json.loads((self.root/'prepared/reference.json').read_text())
        self.assertEqual(data['transcripts'],2);self.assertEqual(data['decoy_sequences'],1)
        self.assertEqual((self.root/'prepared/decoys.txt').read_text(),'1\n')
    def test_patch_annotation_scope_and_gene_identity(self):
        for name,text in {'cdna.gz':'>ENSMUST1.1 gene:ENSMUSG1.1\nACGTACGT\n',
                          'ncrna.gz':'>ENSMUST2.1 gene:ENSMUSG2.1\nTTTTGGGG\n',
                          'gtf.gz':'1\tx\ttranscript\t1\t8\t.\t+\t.\tgene_id "ENSMUSG1"; transcript_id "ENSMUST1";\n'}.items():
            with gzip.open(self.root/name,'wt') as f:f.write(text)
        inputs=[self.root/x for x in ['cdna.gz','ncrna.gz','gtf.gz']]
        with self.assertRaisesRegex(ValueError,'Transcript missing'):
            prepare(*inputs,self.root/'prepared')
        with gzip.open(self.root/'gtf.gz','at') as f:
            f.write('CHR_TEST_PATCH\tx\ttranscript\t1\t8\t.\t+\t.\tgene_id "ENSMUSG2"; transcript_id "ENSMUST2";\n')
        prepare(*inputs,self.root/'prepared')
        self.assertEqual(json.loads((self.root/'prepared/reference.json').read_text())['transcripts'],2)
        with gzip.open(self.root/'ncrna.gz','wt') as f:
            f.write('>ENSMUST2.1 gene:ENSMUSG_WRONG.1\nTTTTGGGG\n')
        with self.assertRaisesRegex(ValueError,'gene identity mismatch'):
            prepare(*inputs,self.root/'mismatched')

    def test_quantification_format(self):
        p=self.root/'quant.sf';p.write_text('Name\tLength\tEffectiveLength\tTPM\tNumReads\ntx1\t100\t70\t1000000\t100\n')
        self.assertEqual(inspect_quant(p),(['tx1'],100.0))
        p.write_text('garbage\nnotquant\n')
        with self.assertRaises(ValueError):inspect_quant(p)
    def test_unresolved_library_type_fails(self):
        for r in self.rows:
            folder=self.root/r['sample_id'];(folder/'aux_info').mkdir(parents=True)
            (folder/'aux_info/meta_info.json').write_text(json.dumps({'num_processed':100,'num_mapped':100,'library_types':['A']}))
            (folder/'quant.sf').write_text('Name\tLength\tEffectiveLength\tTPM\tNumReads\ntx1\t100\t70\t1000000\t100\n')
        with self.assertRaises(ValueError):check(ROOT/'metadata/samples.tsv',self.root,self.root/'qc.json',50,50)
        self.assertFalse(json.loads((self.root/'qc.json').read_text())['passed'])

if __name__=='__main__': unittest.main()
