"""Check complete FASTA target coverage and gene identity against an Ensembl GTF."""
import argparse
import collections
import csv
import gzip
import re
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'workflow/scripts'))
from prepare_reference import fasta
from common import sha256,write_json


def audit(fastas,gtf,output,table):
    annotations={};all_features=set();regions=set()
    with gzip.open(gtf,'rt') as handle:
        for line in handle:
            if line.startswith('#'):continue
            row=line.rstrip().split('\t')
            if len(row)!=9:continue
            attrs=dict(re.findall(r'(\w+) "([^"]*)"',row[8]))
            if 'transcript_id' in attrs:all_features.add(attrs['transcript_id'])
            if row[2]=='transcript':
                tid,gid=attrs['transcript_id'],attrs['gene_id']
                if tid in annotations and annotations[tid][0]!=gid:raise ValueError('Conflicting GTF gene mapping')
                annotations[tid]=(gid,row[0]);regions.add(row[0])
    targets=[]
    for path in fastas:
        for header,seq in fasta(path):
            tid=header.split()[0].split('.')[0]
            match=re.search(r'(?:^|\s)gene:([^\s]+)',header)
            gene=match.group(1).split('.')[0] if match else ''
            locus=re.search(r'(?:chromosome|scaffold):[^:]+:([^:]+):',header)
            annotation=annotations.get(tid)
            targets.append({'transcript_id':tid,'fasta_source':str(path),'fasta_gene_id':gene,
                'fasta_region':locus.group(1) if locus else '',
                'gtf_gene_id':annotation[0] if annotation else '',
                'in_transcript_rows':annotation is not None,'in_any_gtf_feature':tid in all_features,
                'gene_consistent':bool(annotation and gene and gene==annotation[0])})
    Path(table).parent.mkdir(parents=True,exist_ok=True)
    with open(table,'w',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(targets[0]),delimiter='\t');writer.writeheader();writer.writerows(targets)
    missing=[r for r in targets if not r['in_transcript_rows']]
    report={'targets':len(targets),'unique_targets':len({r['transcript_id'] for r in targets}),
        'gtf_transcripts':len(annotations),'gtf_regions':len(regions),
        'missing_targets':len(missing),'missing_in_all_features':sum(not r['in_any_gtf_feature'] for r in missing),
        'missing_by_fasta':dict(collections.Counter(r['fasta_source'] for r in missing)),
        'missing_by_region':dict(collections.Counter(r['fasta_region'] for r in missing)),
        'gene_mismatches_or_missing_header':sum(r['in_transcript_rows'] and not r['gene_consistent'] for r in targets),
        'inputs_sha256':{str(p):sha256(p) for p in [gtf]+fastas}}
    write_json(output,report)
    print({k:v for k,v in report.items() if k not in ['missing_by_region','inputs_sha256']})

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--fasta',nargs='+',required=True);p.add_argument('--gtf',required=True)
    p.add_argument('--output',required=True);p.add_argument('--table',required=True)
    a=p.parse_args();audit(a.fasta,a.gtf,a.output,a.table)
