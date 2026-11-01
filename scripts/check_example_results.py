#!/usr/bin/env python3
"""Check included result tables and HTML resources without raw data or model fitting."""
import argparse
import csv
import hashlib
from html.parser import HTMLParser
import json
import math
from pathlib import Path
from urllib.parse import unquote, urlsplit


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read_rows(path):
    with path.open(newline='', encoding='utf-8') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def number(value):
    return math.nan if value in {'NA', 'nan', ''} else float(value)


class ResourceParser(HTMLParser):
    def __init__(self):
        super().__init__()
        self.references = []

    def handle_starttag(self, tag, attrs):
        self.references.extend(value for key, value in attrs if key in {'src', 'href'} and value)


def check_results(results):
    results = Path(results).resolve()
    primary = []
    for summary in read_rows(results / 'analysis/summary.tsv'):
        directory = results / 'analysis' / summary['level'] / ('day' + summary['day'])
        ids, selected = set(), set()
        counts = {'up': 0, 'down': 0}
        with (directory / 'results.tsv').open(newline='', encoding='utf-8') as handle:
            for line, row in enumerate(csv.DictReader(handle, delimiter='\t'), 2):
                identity = row['feature_id']
                require(identity and identity not in ids, f'{directory}: duplicate/empty feature at line {line}')
                ids.add(identity)
                padj, effect = number(row['padj']), number(row['log2FoldChange'])
                expected = math.isfinite(padj) and padj < float(summary['alpha']) and abs(effect) > float(summary['lfc_cutoff'])
                require((row['significant'] == 'TRUE') == expected, f'{directory}: selection mismatch for {identity}')
                if expected:
                    direction = 'up' if effect > 0 else 'down'
                    require(row['direction'] == direction, f'{directory}: direction mismatch for {identity}')
                    selected.add(identity)
                    counts[direction] += 1
        require(len(ids) == int(summary['features_input']), f'{directory}: input count mismatch')
        for direction in counts:
            require(counts[direction] == int(summary['significant_' + direction]), f'{directory}: {direction} summary mismatch')
        listed = [row['feature_id'] for row in read_rows(directory / 'significant.tsv')]
        require(len(listed) == len(set(listed)) and set(listed) == selected, f'{directory}: significant table mismatch')
        primary.append({'level': summary['level'], 'day': summary['day'], 'rows': len(ids), **counts})
    require({(x['level'], x['day']) for x in primary} == {('gene', '11'), ('gene', '33'), ('transcript', '11'), ('transcript', '33')}, 'Expected four primary contrasts')
    go_count = 0
    for summary in read_rows(results / 'enrichment/summary.tsv'):
        path = results / 'enrichment' / ('day' + summary['day']) / summary['family'] / (summary['direction'] + '.tsv')
        rows = read_rows(path)
        require(len(rows) == int(summary['tested_terms']), f'{path}: tested-term count mismatch')
        significant = sum(number(row['padj']) < float(summary['alpha']) for row in rows)
        require(significant == int(summary['significant_terms']), f'{path}: significant-term count mismatch')
        go_count += 1
    comparison = read_rows(results / 'comparison/summary.tsv')
    for row in comparison:
        expected = int(row['shared_selected']) / int(row['published_selected'])
        require(math.isclose(expected, float(row['published_list_recovery']), rel_tol=1e-12), 'Published-list recovery ratio mismatch')
    report = results / 'report/index.html'
    parser = ResourceParser()
    parser.feed(report.read_text(encoding='utf-8'))
    local_resources = 0
    for reference in parser.references:
        parsed = urlsplit(reference)
        if parsed.scheme or parsed.netloc or not parsed.path:
            continue
        target = (report.parent / unquote(parsed.path)).resolve()
        require(target.is_relative_to(results) and target.exists(), f'Missing or external local report resource: {reference}')
        local_resources += 1
    return {'status': 'PASS', 'primary_tables': primary, 'GO_summary_rows': go_count,
            'comparison_summary_rows': len(comparison), 'local_report_resources': local_resources,
            'scope': 'Saved-table arithmetic and file references; no model refit or browser rendering.'}


def check_snapshot_hashes(repository, results):
    """Verify historical copied outputs using the included release source map."""
    records = read_rows(repository / 'docs/RELEASE_SOURCE_MAP.tsv')
    checked = 0
    for row in records:
        path = row['release_path']
        if not path.startswith('example-results/'):
            continue
        target = results / Path(path).relative_to('example-results')
        require(target.is_file(), f'Missing snapshot file: {path}')
        digest = hashlib.sha256(target.read_bytes()).hexdigest()
        require(digest == row['source_sha256'] == row['release_sha256'], f'Snapshot hash mismatch: {path}')
        checked += 1
    require(checked > 0, 'No snapshot entries found in source map')
    return checked


def main():
    repository = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results', type=Path, default=repository / 'example-results', help='Result snapshot to check')
    args = parser.parse_args()
    summary = check_results(args.results)
    summary['snapshot_files_hashed'] = check_snapshot_hashes(repository, args.results.resolve())
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
