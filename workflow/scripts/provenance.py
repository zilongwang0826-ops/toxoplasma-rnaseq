"""Record the source revision, configuration, references and resolved container identities."""
import argparse
import json
from pathlib import Path
import platform
import subprocess
from common import sha256, utc_now, write_json

if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',required=True)
    p.add_argument('--config-json',required=True)
    p.add_argument('--reference',required=True)
    p.add_argument('--qc',required=True)
    a=p.parse_args()
    effective=json.loads(a.config_json)
    commit=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()
    dirty=bool(subprocess.check_output(['git','status','--porcelain','--untracked-files=no'],text=True).strip())
    files=[]
    for folder in ['config','workflow','containers','metadata','profiles']:
        files.extend(x for x in Path(folder).rglob('*') if x.is_file() and '__pycache__' not in str(x) and x.suffix in {'.py','.R','.yaml','.smk','.json','.tsv','.sh'})
    launches=[]
    for path in sorted(Path(effective['output'],'launch').glob('*.txt')):
        lines=path.read_text().splitlines()
        if lines and len(lines[0])==40 and all(x in '0123456789abcdef' for x in lines[0]):
            launches.append({'record':str(path),'commit':lines[0],'command_and_settings':lines[1:]})
    write_json(a.output, {'completed_at':utc_now(),'commit':commit,'tracked_source_dirty':dirty,
        'commit_scope':'Source state when this provenance record was generated; earlier stages may use earlier launch revisions.',
        'launch_records':launches,
        'stage_provenance_note':'Retain launch logs and Snakemake per-output metadata. Do not attribute every earlier input or stage to the final report commit.',
        'platform':platform.platform(),'python':platform.python_version(),
        'input_and_source_sha256':{str(x):sha256(x) for x in sorted(files)},
        'effective_config':effective,
        'reference':json.loads(Path(a.reference).read_text()),
        'quantification_qc':json.loads(Path(a.qc).read_text()),
        'download_receipts':{str(x):json.loads(x.read_text()) for x in Path(effective['raw']).glob('*.source.json')},
        'containers':json.loads(Path('containers/images.lock.json').read_text()),
        'interpretation':'Execution provenance; scientific agreement must be assessed independently.'})
