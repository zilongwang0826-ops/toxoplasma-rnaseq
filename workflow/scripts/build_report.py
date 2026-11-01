"""Render a portable, dependency-free scientific report from executed result tables."""
import argparse
import csv
import html
import os
from pathlib import Path


def table(path):
    with open(path,newline='') as f:
        rows=list(csv.reader(f,delimiter='\t'))
    if not rows:
        raise ValueError('Empty expected table: '+str(path))
    header='<tr>'+''.join('<th>'+html.escape(x)+'</th>' for x in rows[0])+'</tr>'
    body=''.join('<tr>'+''.join('<td>'+html.escape(x)+'</td>' for x in row)+'</tr>' for row in rows[1:])
    return '<div class="table"><table>'+header+body+'</table></div>'


def render(results,output):
    root=Path(results); out=Path(output)
    required=[root/'analysis/summary.tsv',root/'enrichment/summary.tsv',root/'comparison/summary.tsv']
    for p in required:
        if not p.exists():
            raise FileNotFoundError('Report requires executed output: '+str(p))
    sections=['<h1>Mouse brain response to Toxoplasma gondii infection</h1>',
        '<p>Reanalysis of PRJNA483261; infected versus control at 11 and 33 days. '
        '<a href="https://doi.org/10.3389/fmicb.2020.570903">Hu et al. (2020)</a>.</p>',
        '<p>Three reported biological libraries per group. The publication’s pooling wording remains unresolved. '
        'Transcript and gene tests are separate; BH FDR and effect-size criteria are documented in the analysis parameters. '
        'No independent experimental validation is claimed.</p>']
    for title,path in zip(['Differential expression','Gene-level functional enrichment','Comparison with published transcripts'],required):
        sections += ['<h2>'+title+'</h2>',table(path)]
    sensitivity=root/'sensitivity/lilrb4a/summary.tsv'
    if sensitivity.exists():
        sections += ['<h2>Overlapping-locus counting sensitivity</h2>',
            '<p>Two overlapping Ensembl records at Lilrb4a are counted once in this sensitivity. '
            'Original differential-expression estimates are unchanged. This is a counting sensitivity, not a locus-level differential-expression test. '
            'Term-size filtering and BH correction are repeated for every day/direction family.</p>',table(sensitivity)]
    sections += ['<h2>Quality control and scientific visualizations</h2>']
    images=[]
    for path in sorted((root/'analysis').rglob('*.png')):
        refined=root/'figures'/path.relative_to(root/'analysis')
        images.append(refined if refined.exists() else path)
    for path in sorted((root/'enrichment').rglob('*.png')):
        refined=root/'figures'/path.relative_to(root)
        images.append(refined if refined.exists() else path)
    images+=sorted((root/'comparison').rglob('*.png'))
    for path in images:
        rel=os.path.relpath(path,out.parent)
        sections.append('<figure><img loading="lazy" src="'+html.escape(rel,quote=True)+'"><figcaption>'+html.escape(str(path.relative_to(root)))+'</figcaption></figure>')
    sections += ['<p><a href="../qc/multiqc_report.html">MultiQC report</a> · <a href="../provenance/run.json">Run provenance</a> · '
        '<a href="../validation/models/">Saved-model diagnostics</a> · <a href="../sensitivity/lilrb4a/">Locus sensitivity tables</a></p>',
        '<h2>Interpretation limits</h2><p>Whole-brain changes may reflect both cell composition and within-cell regulation. '
        'Unrecorded batch effects and pooling uncertainty cannot be resolved statistically. '
        'Enrichment is exploratory and does not establish pathway activation or causality. '
        'Published hit counts are reference observations, not targets used for optimization. '
        'GO is the primary functional analysis; inspect kegg-status.txt for KEGG execution status.</p>']
    page='<!doctype html><html lang="en"><meta charset="utf-8"><title>Toxoplasma RNA-seq report</title><style>body{font:16px system-ui;max-width:1100px;margin:40px auto;padding:0 20px;color:#172c3c}h1,h2{color:#123e50}table{border-collapse:collapse;font-size:13px}th,td{padding:8px;border-bottom:1px solid #ccd6df;text-align:left}.table{overflow:auto}img{max-width:100%}figure{margin:30px 0}figcaption{font-size:13px;color:#456}</style><main>'+''.join(sections)+'</main></html>'
    out.parent.mkdir(parents=True,exist_ok=True);out.write_text(page)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--results',required=True);p.add_argument('--output',required=True)
    a=p.parse_args();render(a.results,a.output)
