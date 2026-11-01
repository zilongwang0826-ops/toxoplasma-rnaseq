"""Failure controls for automatic progression from a technical pilot."""
import copy
from pathlib import Path
import sys
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from pilot_gate import assess_metrics

class PilotGate(unittest.TestCase):
    def setUp(self):
        self.sample={'read_pairs':'10000'}
        self.meta={'num_processed':10000,'num_mapped':9000,'library_types':['IU']}
        self.qc={str(m):{'total_sequences':10000,'modules':{k:'PASS' for k in ['Per base sequence quality','Per base N content','Adapter Content']}} for m in (1,2)}
        self.thresholds={'min_mapping_percent':50,'min_assigned_fragments':100}
    def evaluate(self):return assess_metrics(self.sample,self.meta,9000,self.qc,self.thresholds)
    def test_passing_pilot(self):self.assertEqual(self.evaluate(),([],90.0,['IU']))
    def test_incomplete_quantification(self):
        self.meta['num_processed']=9999
        self.assertTrue(self.evaluate()[0])
    def test_adapter_warning_requires_review(self):
        self.qc['2']['modules']['Adapter Content']='WARN'
        self.assertTrue(self.evaluate()[0])
    def test_unresolved_library_type(self):
        self.meta['library_types']=['A']
        self.assertTrue(self.evaluate()[0])
    def test_missing_mate(self):
        del self.qc['1']
        self.assertTrue(self.evaluate()[0])
    def test_low_mapping(self):
        self.meta['num_mapped']=100
        self.assertTrue(self.evaluate()[0])
