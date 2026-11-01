#!/usr/bin/env python3
"""Check the integrity and basic schema of the retained small scientific sources."""
import csv
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent


def require(condition, message):
    if not condition:
        raise ValueError(message)


def main():
    manifest = json.loads((ROOT / 'source/sources.json').read_text())
    seen = set()
    for record in manifest:
        name = record['file']
        require(name not in seen, 'Duplicate source-manifest entry: ' + name)
        seen.add(name)
        path = ROOT / 'source' / name
        require(path.is_file(), 'Source missing: ' + name)
        data = path.read_bytes()
        require(len(data) == record['bytes'], 'Size mismatch: ' + name)
        require(hashlib.sha256(data).hexdigest() == record['sha256'], 'SHA-256 mismatch: ' + name)
        if 'md5' in record:
            require(hashlib.md5(data).hexdigest() == record['md5'], 'MD5 mismatch: ' + name)
    rows = list(csv.DictReader((ROOT / 'samples.tsv').open(), delimiter='\t'))
    required = {'sample_id','run_accession','biosample','condition','day','replicate','layout',
                'fastq_1','fastq_2','md5_1','md5_2','bytes_1','bytes_2'}
    require(rows and required <= rows[0].keys(), 'Sample schema incomplete')
    require(len(rows) == 12, 'Expected twelve study libraries')
    published = list(csv.DictReader((ROOT / 'published_dets.tsv').open(), delimiter='\t'))
    require(len(published) == 3017, 'Unexpected published DET row count')
    for day, count in [('11',936),('33',2081)]:
        selected = [row for row in published if row['day']==day]
        require(len(selected)==count, 'Published count mismatch for day ' + day)
        require(len({row['feature_id'] for row in selected})==count, 'Duplicate published transcript')
        require(all(row['level']=='transcript' and row['pvalue']=='NA' for row in selected),
                'Unexpected analysis level or fabricated nominal P value')
    print(json.dumps({'status':'passed','source_files_verified':len(seen),
                      'samples_verified':len(rows),'published_dets_verified':len(published)},indent=2))


if __name__ == '__main__':
    try:
        main()
    except (OSError,ValueError,KeyError) as error:
        sys.exit('Metadata verification failed: ' + str(error))
