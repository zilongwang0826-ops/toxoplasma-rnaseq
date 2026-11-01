"""Bind successful synthetic execution evidence to source and container identities."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
from datetime import datetime,timezone


def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def source_fingerprints():
    files=[Path('workflow/Snakefile'),Path('tests/smoke.smk')]
    for directory in ['workflow/scripts','workflow/rules','tests','scripts']:
        files.extend(p for p in Path(directory).glob('*') if p.is_file() and p.suffix in {'.py','.R','.smk','.sh'})
    return {str(p):digest(p) for p in sorted(set(files))}

p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args()
marker=Path('results/container-validation/validated.json')
current={'image_lock_sha256':digest('containers/images.lock.json'),'source_sha256':source_fingerprints()}
if a.check:
    prior=json.loads(marker.read_text())
    if any(prior[k]!=v for k,v in current.items()):raise SystemExit('Container/source changed since tests; repeat container-test')
else:
    if not Path('results/tests/locus-sensitivity/verified.txt').read_text().startswith('PASS:'):
        raise SystemExit('Locus-counting sensitivity tests did not pass')
    go=json.loads(Path('results/tests/go-integration/verified.json').read_text())
    if go['status']!='PASS':raise SystemExit('Full GO integration did not pass')
    integrity=json.loads(Path('results/tests/python/verified.json').read_text())
    if integrity['status']!='PASS':raise SystemExit('Input integrity tests did not pass')
    quant=json.loads(Path('results/tests/salmon/verified.json').read_text())
    if quant['status']!='PASS':raise SystemExit('Synthetic quantification did not pass')
    with open('results/tests/container-analysis/integration-status.tsv') as f: rows=list(csv.DictReader(f,delimiter='\t'))
    if {r['component'] for r in rows}!={'tximport_file_integration','GO_database_integration'} or any(not r['status'].startswith('PASS') for r in rows):raise SystemExit('Required R integrations not validated')
    current['validated_at']=datetime.now(timezone.utc).isoformat()
    current['evidence']=['results/tests/salmon/verified.json','results/tests/container-analysis/integration-status.tsv']
    marker.parent.mkdir(parents=True,exist_ok=True);marker.write_text(json.dumps(current,indent=2)+'\n')
print('Container validation record verified' if a.check else 'Container validation recorded')
