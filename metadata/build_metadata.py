#!/usr/bin/env python3
"""Rebuild validated sample metadata from archived ENA, BioSample and article XML."""
import csv
import hashlib
import json
from pathlib import Path
import re
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / 'source'


def write_tsv(path, rows, columns):
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def main():
    runs = list(csv.DictReader((SOURCE / 'ena_runs.tsv').open(), delimiter='\t'))
    biosamples = {sample.attrib['accession']: {
        attribute.attrib['attribute_name']: attribute.text
        for attribute in sample.findall('.//Attribute')
    } for sample in ET.parse(SOURCE / 'biosamples.xml').getroot().findall('BioSample')}
    paper = ET.parse(SOURCE / 'epmc_article.xml').getroot()
    reported = {}
    for table in paper.findall('.//table-wrap'):
        if table.findtext('label') != 'TABLE 1':
            continue
        for row in table.findall('.//tr'):
            cells = [''.join(cell.itertext()).strip() for cell in row]
            indices = [index for index, cell in enumerate(cells) if re.fullmatch('[CT](11|33)[123]', cell)]
            if not indices:
                continue
            index = indices[0]
            reported[cells[index]] = dict(raw_reads=int(cells[index + 1].replace(',', '')),
                                         clean_reads=int(cells[index + 2].replace(',', '')))
    assert len(reported) == 12, 'Expected exactly twelve paper sample codes'
    rows, comparisons = [], []
    for run in runs:
        attrs = biosamples[run['sample_accession']]
        match = re.fullmatch('(acutely|chronically) (uninfected|infected) group ([123])', run['sample_alias'])
        assert match, 'Unexpected ENA alias'
        phase, status, replicate = match.groups()
        day = 11 if phase == 'acutely' else 33
        condition = 'control' if status == 'uninfected' else 'infected'
        sample_id = f"{'C' if condition == 'control' else 'T'}{day}{replicate}"
        assert attrs['biological_replicate'] == f'biological replicate {replicate}'
        assert attrs['sample_type'] == f"{'control' if condition == 'control' else 'test'} group on {day} day post-infection"
        assert run['library_layout'] == 'PAIRED'
        urls, md5s, sizes = [run[column].split(';') for column in ('fastq_ftp', 'fastq_md5', 'fastq_bytes')]
        assert len(urls) == len(md5s) == len(sizes) == 2
        assert int(run['base_count']) == int(run['read_count']) * 300
        assert reported[sample_id]['clean_reads'] == int(run['read_count']) * 2
        assert all(re.fullmatch('[0-9a-f]{32}', checksum) for checksum in md5s)
        rows.append(dict(sample_id=sample_id, run_accession=run['run_accession'],
                         biosample=run['sample_accession'], condition=condition, day=day,
                         replicate=replicate, layout='PAIRED', fastq_1='https://' + urls[0],
                         fastq_2='https://' + urls[1], md5_1=md5s[0], md5_2=md5s[1],
                         bytes_1=int(sizes[0]), bytes_2=int(sizes[1]),
                         experiment_accession=run['experiment_accession'], sample_alias=run['sample_alias'],
                         read_pairs=int(run['read_count']), base_count=int(run['base_count']),
                         read_length=150, instrument_model=run['instrument_model'],
                         library_selection=run['library_selection'], sex=attrs['sex'], strain=attrs['strain'],
                         tissue=attrs['tissue'], age_archive=attrs['age'],
                         batch='unknown', strandedness='unknown',
                         experimental_unit='reported biological replicate; pooling unresolved',
                         archived_processing='matches published clean-read counts'))
        comparisons.append(dict(sample_id=sample_id,run_accession=run['run_accession'],
                                archive_read_pairs=int(run['read_count']),
                                archive_individual_reads=2*int(run['read_count']),
                                paper_raw_reads=reported[sample_id]['raw_reads'],
                                paper_clean_reads=reported[sample_id]['clean_reads'],
                                clean_count_match=True))
    rows.sort(key=lambda row:(row['day'],row['condition'],int(row['replicate'])))
    comparisons.sort(key=lambda row:row['sample_id'])
    assert len(rows) == len({row['sample_id'] for row in rows}) == 12
    assert len({row['biosample'] for row in rows}) == len({row['run_accession'] for row in rows}) == 12
    assert all(sum(row['condition']==condition and row['day']==day for row in rows)==3
               for condition in ('control','infected') for day in (11,33))
    write_tsv(ROOT/'samples.tsv',rows,list(rows[0]))
    write_tsv(ROOT/'archive_paper_read_counts.tsv',comparisons,list(comparisons[0]))
    summary = dict(status='metadata_cross_checked',samples=len(rows),groups=4,replicates_per_group=3,
                   archived_fastq_bytes=sum(row['bytes_1']+row['bytes_2'] for row in rows),
                   archived_read_pairs=sum(row['read_pairs'] for row in rows),
                   archived_bases=sum(row['base_count'] for row in rows),
                   all_archive_read_counts_match_published_clean_counts=True,
                   unresolved=['pooling language','sequencer model discrepancy','batch information','library strandedness'],
                   source_sha256={name:hashlib.sha256((SOURCE/name).read_bytes()).hexdigest()
                                  for name in ('ena_runs.tsv','biosamples.xml','epmc_article.xml')})
    (ROOT/'validation.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))


if __name__ == '__main__':
    main()
