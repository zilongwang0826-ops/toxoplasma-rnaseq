#!/usr/bin/env python3
"""Extract published DETs without modifying the original XLSX supplement.

Only stored cell values are read. Formula cells fail explicitly. Python's standard
library is sufficient, so rebuilding the provenance tables needs no Excel engine.
"""
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import xml.etree.ElementTree as ET
import zipfile

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / 'source' / 'Table_1.XLSX'
NS = {'s': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
HEADERS = ['Transcript ID', 'Gene ID', 'Gene Symbol', 'Log2 (FoldChange)',
           'Catergery', 'Adjusted P-value', 'Function Description']


def read_sheets(path):
    with zipfile.ZipFile(path) as archive:
        strings = []
        if 'xl/sharedStrings.xml' in archive.namelist():
            strings = [''.join(item.itertext()) for item in
                       ET.fromstring(archive.read('xl/sharedStrings.xml')).findall('s:si', NS)]
        relations = {item.attrib['Id']: item.attrib['Target'] for item in
                     ET.fromstring(archive.read('xl/_rels/workbook.xml.rels'))}
        workbook = ET.fromstring(archive.read('xl/workbook.xml'))
        for sheet in workbook.findall('s:sheets/s:sheet', NS):
            rid = sheet.attrib['{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id']
            target = relations[rid]
            filename = target.lstrip('/') if target.startswith('/') else 'xl/' + target
            rows = []
            for row in ET.fromstring(archive.read(filename)).findall('s:sheetData/s:row', NS):
                values = [''] * len(HEADERS)
                for cell in row.findall('s:c', NS):
                    if cell.find('s:f', NS) is not None:
                        raise ValueError('Formula found; stored values cannot establish recalculation')
                    letters = re.match('[A-Z]+', cell.attrib['r']).group()
                    column = 0
                    for letter in letters:
                        column = column * 26 + ord(letter) - ord('A') + 1
                    value = cell.findtext('s:v', default='', namespaces=NS)
                    if cell.attrib.get('t') == 's':
                        value = strings[int(value)]
                    elif cell.attrib.get('t') == 'inlineStr':
                        value = ''.join(cell.find('s:is', NS).itertext())
                    if column > len(HEADERS):
                        if value:
                            raise ValueError('Unexpected additional nonempty column')
                    else:
                        values[column - 1] = value
                if any(values):
                    rows.append((int(row.attrib['r']), values))
            yield sheet.attrib['name'], rows


def write_tsv(path, rows, columns):
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def main():
    records = []
    audit = {'source': 'source/Table_1.XLSX', 'source_sha256': hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
             'source_doi': '10.3389/fmicb.2020.570903.s002',
             'interpretation': 'Published selected DET lists, not independently reproduced results',
             'nominal_p_values_available': False, 'sheets': {}}
    for name, sheet_rows in read_sheets(SOURCE):
        if name not in ('Acute', 'Chronic'):
            raise ValueError('Unknown supplementary sheet: ' + name)
        if sheet_rows[0][1] != HEADERS:
            raise ValueError('Source header changed')
        day = 11 if name == 'Acute' else 33
        raw = []
        group = []
        for source_row, values in sheet_rows[1:]:
            tx, gene, symbol, lfc, direction, padj, description = values
            if not re.fullmatch('ENSMUST[0-9]+', tx) or not re.fullmatch('ENSMUSG[0-9]+', gene):
                raise ValueError('Unexpected source identifier format')
            if not math.isfinite(float(lfc)) or not 0 <= float(padj) <= 1:
                raise ValueError('Invalid supplementary numeric value')
            expected_direction = 'Upregulated' if float(lfc) > 0 else 'Downregulated'
            if direction != expected_direction:
                raise ValueError('Direction contradicts LFC')
            raw.append(dict(zip(['source_row'] + HEADERS, [source_row] + values)))
            group.append(dict(level='transcript', feature_id=tx, day=day, log2FoldChange=lfc,
                              pvalue='NA', padj=padj, direction='up' if float(lfc)>0 else 'down',
                              gene_id=gene, gene_symbol=symbol,
                              gene_symbol_status='numeric_in_source' if re.fullmatch(r'[0-9]+(?:\.[0-9]+)?', symbol) else 'text',
                              source_sheet=name, source_row=source_row))
        if len(group) != len({row['feature_id'] for row in group}):
            raise ValueError('Duplicate published transcript identifiers within a contrast')
        write_tsv(ROOT / f'published_table1_{name.lower()}.tsv', raw, ['source_row'] + HEADERS)
        audit['sheets'][name] = dict(day=day,rows=len(group),unique_transcripts=len({row['feature_id'] for row in group}),
                                   unique_genes=len({row['gene_id'] for row in group}),
                                   upregulated=sum(row['direction']=='up' for row in group),
                                   downregulated=sum(row['direction']=='down' for row in group),
                                   min_abs_log2fc=min(abs(float(row['log2FoldChange'])) for row in group),
                                   max_padj=max(float(row['padj']) for row in group),
                                   exact_zero_padj=sum(float(row['padj'])==0 for row in group),
                                   numeric_gene_symbols=sum(row['gene_symbol_status']=='numeric_in_source' for row in group))
        records.extend(group)
    for direction in ('up','down'):
        a={row['feature_id'] for row in records if row['day']==11 and row['direction']==direction}
        b={row['feature_id'] for row in records if row['day']==33 and row['direction']==direction}
        audit[f'shared_{direction}_transcripts']=len(a & b)
    write_tsv(ROOT/'published_dets.tsv',records,list(records[0]))
    (ROOT/'published_dets_validation.json').write_text(json.dumps(audit,indent=2)+'\n')
    print(json.dumps(audit,indent=2))


if __name__ == '__main__':
    main()
