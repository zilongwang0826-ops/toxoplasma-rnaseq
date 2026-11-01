"""Small fixtures verify the offline result checker's failure controls."""
import csv
import importlib.util
import tempfile
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('example_check', ROOT / 'scripts/check_example_results.py')
CHECKER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CHECKER)


def table(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


class ExampleResults(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='example-check-')
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        summaries = []
        for level in ['gene', 'transcript']:
            for day in ['11', '33']:
                summaries.append(dict(level=level, day=day, alpha='0.05', lfc_cutoff='1.5',
                                      features_input='2', significant_up='1', significant_down='0'))
                directory = self.root / 'analysis' / level / ('day' + day)
                rows = [dict(feature_id='feature1', padj='0.01', log2FoldChange='2', significant='TRUE', direction='up'),
                        dict(feature_id='feature2', padj='NA', log2FoldChange='NA', significant='FALSE', direction='not_significant')]
                table(directory / 'results.tsv', rows)
                table(directory / 'significant.tsv', rows[:1])
        table(self.root / 'analysis/summary.tsv', summaries)
        table(self.root / 'enrichment/summary.tsv', [dict(day='11', family='GO', direction='up', tested_terms='1', significant_terms='1', alpha='0.05')])
        table(self.root / 'enrichment/day11/GO/up.tsv', [dict(term_id='GO:test', padj='0.01')])
        table(self.root / 'comparison/summary.tsv', [dict(shared_selected='1', published_selected='2', published_list_recovery='0.5')])
        (self.root / 'report').mkdir()
        (self.root / 'report/index.html').write_text('<a href="../analysis/summary.tsv">Summary</a>')

    def test_consistent_snapshot_passes(self):
        self.assertEqual(CHECKER.check_results(self.root)['status'], 'PASS')

    def test_corrupt_selection_fails(self):
        path = self.root / 'analysis/gene/day11/results.tsv'
        rows = CHECKER.read_rows(path)
        rows[0]['significant'] = 'FALSE'
        table(path, rows)
        with self.assertRaisesRegex(ValueError, 'selection mismatch'):
            CHECKER.check_results(self.root)

    def test_missing_report_resource_fails(self):
        (self.root / 'report/index.html').write_text('<img src="missing.png">')
        with self.assertRaisesRegex(ValueError, 'Missing or external'):
            CHECKER.check_results(self.root)

    def test_wrong_recovery_ratio_fails(self):
        table(self.root / 'comparison/summary.tsv', [dict(shared_selected='1', published_selected='2', published_list_recovery='0.9')])
        with self.assertRaisesRegex(ValueError, 'recovery ratio mismatch'):
            CHECKER.check_results(self.root)
