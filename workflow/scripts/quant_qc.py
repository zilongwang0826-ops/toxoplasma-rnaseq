"""Gate inference on matching transcript identities, library type and usable read depth."""
import argparse
import csv
import json
import math
from pathlib import Path
from common import read_samples, write_json, utc_now


def inspect_quant(path):
    with open(path,newline='') as handle:
        reader=csv.DictReader(handle,delimiter='\t')
        if reader.fieldnames != ['Name','Length','EffectiveLength','TPM','NumReads']:
            raise ValueError('Malformed Salmon quantification header: '+str(path))
        names=[]; total=0.0; tpm=0.0
        for row in reader:
            names.append(row['Name'])
            values=[float(row[k]) for k in ['Length','EffectiveLength','TPM','NumReads']]
            if not all(math.isfinite(x) and x>=0 for x in values):
                raise ValueError('Invalid quantification numbers')
            if values[0]<=0:
                raise ValueError('Nonpositive transcript length')
            total+=values[3];tpm+=values[2]
        if not names or any(not x for x in names) or len(names)!=len(set(names)) or total<=0:
            raise ValueError('Empty quantification or duplicated transcript identity')
        if abs(tpm-1e6)>100:
            raise ValueError('TPM total is inconsistent with a complete Salmon output')
    return names,total


def check(samples, root, out, min_percent, min_fragments, require_archive_counts=False):
    rows=[]; errors=[]; inferred=set(); reference_ids=None
    for sample in read_samples(samples):
        folder=Path(root)/sample['sample_id']
        meta=json.loads((folder/'aux_info/meta_info.json').read_text())
        types=meta.get('library_types',[])
        if len(types)!=1 or types[0] not in {'IU','ISF','ISR','OU','OSF','OSR','MU','MSF','MSR'}:
            errors.append(sample['sample_id']+': paired library type unresolved or invalid')
        inferred.update(types)
        processed,mapped=int(meta['num_processed']),int(meta['num_mapped'])
        if require_archive_counts and processed != int(sample['read_pairs']):
            errors.append(sample['sample_id']+': processed pair count differs from the complete archive')
        pct=100*mapped/processed if processed else 0
        if processed<=0 or mapped>processed or mapped<0 or pct<min_percent or mapped<min_fragments:
            errors.append(sample['sample_id']+': insufficient or invalid mapping/depth')
        ids,total=inspect_quant(folder/'quant.sf')
        if reference_ids is None:
            reference_ids=ids
        elif ids!=reference_ids:
            errors.append(sample['sample_id']+': inconsistent transcript reference/order')
        rows.append({'sample_id':sample['sample_id'],'processed_pairs':processed,'mapped_pairs':mapped,
                     'assigned_fragments':total,'mapping_percent':pct,'library_types':';'.join(types),'transcripts':len(ids)})
        if total<min_fragments:
            errors.append(sample['sample_id']+': insufficient assigned fragments')
    if len(inferred)>1:
        errors.append('Inconsistent inferred library types; inspect before inference')
    Path(out).parent.mkdir(parents=True,exist_ok=True)
    with open(str(out)+'.tsv','w',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(rows[0]),delimiter='\t');writer.writeheader();writer.writerows(rows)
    write_json(out,{'checked_at':utc_now(),'samples':rows,'errors':errors,'passed':not errors})
    if errors:
        raise ValueError('Quantification QC failed: '+'; '.join(errors))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--samples',required=True);p.add_argument('--quant-root',required=True)
    p.add_argument('--output',required=True);p.add_argument('--min-percent',type=float,default=50)
    p.add_argument('--min-fragments',type=int,default=1000000)
    p.add_argument('--require-archive-counts',action='store_true')
    a=p.parse_args();check(a.samples,a.quant_root,a.output,a.min_percent,a.min_fragments,a.require_archive_counts)
