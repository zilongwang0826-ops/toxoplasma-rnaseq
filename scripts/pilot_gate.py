"""Validate a complete real-library pilot before the approved full study run."""
import argparse
import json
from pathlib import Path
import sys
import zipfile

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'workflow/scripts'))
from common import read_samples,write_json,utc_now
from quant_qc import inspect_quant


def assess_metrics(sample,meta,assigned,qc,thresholds):
    errors=[]
    expected=int(sample['read_pairs'])
    for mate,report in qc.items():
        if report['total_sequences']!=expected:
            errors.append(f'Mate {mate}: FastQC count differs from archive')
        for module in ['Per base sequence quality','Per base N content','Adapter Content']:
            status=report['modules'].get(module)
            if status!='PASS':errors.append(f'Mate {mate}: review {module}: {status}')
    if set(qc)!={'1','2'}:errors.append('Both mate QC reports are required')
    processed,mapped=int(meta['num_processed']),int(meta['num_mapped'])
    percent=100*mapped/processed if processed else 0
    if processed!=expected:errors.append('Salmon pair count differs from archive')
    if mapped<0 or mapped>processed:errors.append('Invalid mapped fragment count')
    types=meta.get('library_types',[])
    if len(types)!=1 or types[0] not in {'IU','ISF','ISR','OU','OSF','OSR','MU','MSF','MSR'}:
        errors.append('Unresolved or invalid paired library type')
    if percent<thresholds['min_mapping_percent'] or assigned<thresholds['min_assigned_fragments']:
        errors.append('Pilot mapping or assigned depth below prespecified gate')
    return errors,percent,types


def gate(config_path,sample_id):
    import yaml
    cfg=yaml.safe_load(Path(config_path).read_text())
    samples=read_samples(cfg['samples'])
    sample=next(s for s in samples if s['sample_id']==sample_id)
    root=Path(cfg['output']);raw=Path(cfg['raw'])
    if cfg['quantification']['trim']:
        raise ValueError('The automatic pilot gate is defined for the untrimmed archive analysis')
    errors=[];qc={}
    for mate in (1,2):
        name=f'{sample_id}_{mate}'
        reads=raw/(name+'.fastq.gz')
        receipt=json.loads(Path(str(reads)+'.source.json').read_text())
        if reads.stat().st_size!=int(sample[f'bytes_{mate}']) or receipt['md5']!=sample[f'md5_{mate}']:
            errors.append(f'Mate {mate}: archive provenance mismatch')
        with zipfile.ZipFile(root/'qc/raw'/(name+'_fastqc.zip')) as z:
            summary=z.read(name+'_fastqc/summary.txt').decode()
            data=z.read(name+'_fastqc/fastqc_data.txt').decode()
        states={line.split('\t')[1]:line.split('\t')[0] for line in summary.splitlines()}
        totals=[line.split('\t')[1] for line in data.splitlines() if line.startswith('Total Sequences\t')]
        qc[str(mate)]={'modules':states,'total_sequences':int(totals[0]) if len(totals)==1 else 0}
    folder=root/'salmon'/sample_id
    meta=json.loads((folder/'aux_info/meta_info.json').read_text())
    ids,assigned=inspect_quant(folder/'quant.sf')
    metric_errors,percent,types=assess_metrics(sample,meta,assigned,qc,cfg['quantification'])
    errors.extend(metric_errors)
    processed,mapped=int(meta['num_processed']),int(meta['num_mapped'])
    output=root/'validation'/('pilot-'+sample_id+'.json')
    write_json(output,{'checked_at':utc_now(),'sample':sample_id,'passed':not errors,'errors':errors,
        'processed_pairs':processed,'mapped_pairs':mapped,'mapping_percent':percent,'assigned_fragments':assigned,
        'library_type':types,'transcripts':len(ids),'fastqc':qc,
        'interpretation':'Technical pilot QC; not a differential-expression or reproduction result.'})
    if errors:raise ValueError('; '.join(errors))
    print('Real-library technical pilot passed:',sample_id,percent,types)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config',default='config/config.yaml');p.add_argument('--sample',default='C111')
    a=p.parse_args();gate(a.config,a.sample)
