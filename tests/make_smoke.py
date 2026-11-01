"""Deterministic synthetic reads for software validation, never biological inference."""
import argparse
import gzip
from pathlib import Path
import random


def reverse_complement(s): return s.translate(str.maketrans('ACGT','TGCA'))[::-1]

def generate(out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    rng=random.Random(483261)
    seqs=[''.join(rng.choices('ACGT',k=2500)) for _ in range(3)]
    (out/'transcripts.fa').write_text(''.join(f'>synthetic_tx{i+1}\n{s}\n' for i,s in enumerate(seqs)))
    # Fixed mixture: 90%, 8%, 2%; paired fragments with a 300-base insert.
    assignments=[0]*9000+[1]*800+[2]*200;rng.shuffle(assignments)
    with (out/'sample_1.fastq.gz').open('wb') as f1,(out/'sample_2.fastq.gz').open('wb') as f2:
        with gzip.GzipFile(fileobj=f1,mode='wb',mtime=0) as r1,gzip.GzipFile(fileobj=f2,mode='wb',mtime=0) as r2:
            for i,target in enumerate(assignments):
                start=rng.randrange(2201);fragment=seqs[target][start:start+300]
                if rng.random()<0.5:fragment=reverse_complement(fragment)
                for stream,sequence,mate in [(r1,fragment[:100],1),(r2,reverse_complement(fragment[-100:]),2)]:
                    stream.write(f'@synthetic_{i}/{mate}\n{sequence}\n+\n'.encode()+b'I'*100+b'\n')
    (out/'README.txt').write_text('SYNTHETIC SOFTWARE TEST ONLY. Three random transcripts; 10,000 pairs; mixture 90/8/2%.\n')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--outdir',required=True);a=p.parse_args();generate(a.outdir)
