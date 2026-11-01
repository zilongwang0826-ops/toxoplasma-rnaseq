#!/usr/bin/env python3
"""Check frozen-table identity and every plotting transformation by stable keys."""
import hashlib,json,re
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from PIL import Image
import pymupdf as fitz
ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'report/figures/redesigned'
DOC=ROOT/'docs/figure-redesign'
checks=[]
def require(condition,name,**detail):
    if not condition: raise AssertionError(name+(': '+json.dumps(detail) if detail else ''))
    checks.append(dict(check=name,status='PASS',**detail))
def read(p):return pd.read_csv(ROOT/p,sep='\t',low_memory=False)
def close(a,b):return np.allclose(np.asarray(a,dtype=float),np.asarray(b,dtype=float),rtol=1e-12,atol=1e-14,equal_nan=True)
def compare(a,b,keys,columns,name):
    require(not a.duplicated(keys).any() and not b.duplicated(keys).any(),name+' unique keys')
    merged=a.merge(b,on=keys,how='left',suffixes=('_plot','_source'),validate='one_to_one',indicator=True)
    require((merged['_merge']=='both').all(),name+' aligned keys',rows=len(a))
    for c in columns:
        require(close(merged[c+'_plot'],merged[c+'_source']),name+' '+c)
baseline=json.loads((DOC/'baseline-tables.json').read_text())
# Release packaging is maintained separately; it is not a scientific result table.
baseline={p:h for p,h in baseline.items() if p!='docs/RELEASE_SOURCE_MAP.tsv'}
require(all(hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h for p,h in baseline.items()),'All frozen baseline tables unchanged',tables=len(baseline))
figbase=json.loads((DOC/'baseline-figures.json').read_text())
oldfig=[e for e in figbase if e['path']!='report/report.pdf']
require(all(hashlib.sha256((ROOT/e['path']).read_bytes()).hexdigest()==e['sha256'] for e in oldfig),'All original standalone figure assets unchanged',files=len(oldfig))
samples=read('metadata/samples.tsv'); require(len(samples)==12 and not samples.sample_id.duplicated().any(),'12 distinct sample identities')
require((samples.groupby(['day','condition']).size()==3).all(),'3 libraries per time and condition')
for level in ['gene','transcript']:
 for scope in ['all','11','33']:
  src=f'example-results/analysis/{level}/'+('qc' if scope=='all' else f'day{scope}/figures')
  stem='F01_sample_structure' if (level,scope)==('gene','all') else f'S_QC_{level}_{scope}'
  a=read(f'report/figures/redesigned/plot-data/{stem}_pca.tsv'); b=read(src+'/pca-coordinates.tsv')
  compare(a,b,['sample_id'],['PC1','PC2'],stem+' PCA')
  m=a.merge(samples[['sample_id','day','condition']],on='sample_id',suffixes=('_plot','_metadata'),validate='one_to_one')
  require((m.day_plot==m.day_metadata).all() and (m.condition_plot==m.condition_metadata).all(),stem+' group identity')
  dist=read(src+'/sample-distances.tsv').set_index('sample_id')
  dd=read(f'report/figures/redesigned/plot-data/{stem}_distances.tsv')
  require(len(dd)==len(dist)**2 and close(dd.distance,[dist.loc[r,c] for r,c in zip(dd.row,dd.column)]),stem+' full distance matrix')
  hm=read(src+'/heatmap-row-zscores.tsv').set_index('feature_id')
  hd=read(f'report/figures/redesigned/plot-data/S_heatmap_{level}_{scope}.tsv')
  require(len(hd)==hm.size and close(hd.z_score,[hm.loc[f,s] for f,s in zip(hd.feature_id,hd.sample_id)]),f'{level} {scope} all saved heatmap values',features=len(hm),libraries=len(hm.columns))
 for plot in ['volcano','MA']:
  stem=('F02_gene_differential_expression' if level=='gene' else 'S_transcript_differential_expression') if plot=='volcano' else f'S_MA_{level}'
  a=read(f'report/figures/redesigned/plot-data/{stem}.tsv')
  b=pd.concat([read(f'example-results/analysis/{level}/day{d}/results.tsv').assign(day=d) for d in [11,33]])
  compare(a,b,['day','feature_id'],['log2FoldChange','padj','baseMean'],stem)
  eligible=b[np.isfinite(b.log2FoldChange)&np.isfinite(b.padj)] if plot=='volcano' else b[np.isfinite(b.log2FoldChange)&np.isfinite(b.baseMean)&(b.baseMean>0)]
  require(len(a)==len(eligible),stem+' no eligible points omitted')
  if plot=='volcano':
   require(close(a.neglog10_fdr,-np.log10(np.maximum(a.padj,np.finfo(float).tiny))),stem+' adjusted-P transform')
  else:
   expect=np.where(~np.isfinite(a.padj),'Adjusted P unavailable',np.where(a.padj<.05,'BH-adjusted P < 0.05','Other tested'))
   require((a.fdr_status==expect).all(),stem+' original BH-only color meaning')
  mask=(b.padj<.05)&(b.log2FoldChange.abs()>1.5)
  require((mask==b.significant).all(),stem+' unchanged primary thresholds')
  observed=b[mask].groupby(['day','direction']).size().to_dict()
  expected={(11,'up'):496,(11,'down'):41,(33,'up'):1173,(33,'down'):18} if level=='gene' else {(11,'up'):724,(11,'down'):130,(33,'up'):1804,(33,'down'):101}
  require(observed==expected,stem+' original selection counts')
# GO figures preserve effect size, raw P, adjusted P, overlap and universe exactly.
cols=['overlap','term_size','selected_size','universe_size','enrichment_ratio','pvalue','padj']
go=pd.concat([read(f'example-results/enrichment/day{d}/GO/up.tsv').assign(day=d) for d in [11,33]])
a=read('report/figures/redesigned/plot-data/F03_immune_programs.tsv')
compare(a,go,['day','term_id'],cols,'Main GO themes')
require(len(a)==12 and a.groupby('day').term_id.apply(set).iloc[0]==a.groupby('day').term_id.apply(set).iloc[1],'Same six GO themes at both times')
for kind in ['primary','collapsed']:
 for day in [11,33]:
  src=f'example-results/enrichment/day{day}/GO/up.tsv' if kind=='primary' else f'example-results/sensitivity/lilrb4a/day{day}/up/collapsed.tsv'
  b=read(src);a=read(f'report/figures/redesigned/plot-data/S_GO_{kind}_{day}.tsv')
  compare(a,b,['term_id'],cols,f'GO {kind} {day}')
  expected=b[b.padj<.05].sort_values(['padj','pvalue','term_id']).head(15)
  require(list(a.term_id)==list(expected.term_id),f'GO {kind} {day} transparent top-15 selection')
down=read('report/figures/redesigned/plot-data/S_GO_downward_status.tsv')
for row in down.itertuples():
 src=f'example-results/enrichment/day{row.day}/GO/down.tsv' if row.counting=='primary' else f'example-results/sensitivity/lilrb4a/day{row.day}/down/collapsed.tsv'
 b=read(src)
 require(row.terms==len(b) and row.significant==int((b.padj<.05).sum())==0 and close([row.minimum_fdr],[b.padj.min()]),f'Downward GO {row.counting} {row.day}')
 require(row.selected_genes==b.selected_size.iloc[0] and row.universe==b.universe_size.iloc[0],f'Downward universe {row.counting} {row.day}')
pairs=read('report/figures/redesigned/plot-data/F04_published_agreement_pairs.tsv')
recovery=read('report/figures/redesigned/plot-data/F04_published_agreement_recovery.tsv')
for day in [11,33]:
 b=read(f'example-results/comparison/day{day}-published-transcript-comparison.tsv');a=pairs[pairs.day==day]
 compare(a,b,['feature_id'],['log2FoldChange_published','log2FoldChange_reanalysis'],'Published pairs '+str(day))
 bfinite=b[np.isfinite(b.log2FoldChange_published)&np.isfinite(b.log2FoldChange_reanalysis)]
 s=recovery[recovery.day==day].iloc[0]
 require(len(a)==len(bfinite)==s.finite_effect_pairs,'All finite published pairs '+str(day))
 require((a.log2FoldChange_published.abs()<=14).all() and (a.log2FoldChange_reanalysis.abs()<=14).all(),'Agreement limits include every point '+str(day))
 require(close([spearmanr(a.log2FoldChange_published,a.log2FoldChange_reanalysis).statistic],[s.effect_spearman]),'Conditional Spearman '+str(day))
 require(close([(np.sign(a.log2FoldChange_published)==np.sign(a.log2FoldChange_reanalysis)).mean()],[s.same_direction_fraction]),'Conditional direction agreement '+str(day))
 require(s.shared_selected==int(b.significant.sum()) and s.published_selected==len(b) and close([s.shared_selected/s.published_selected],[s.published_list_recovery]),'Published-list recovery denominator '+str(day))
# File-format and actual dimensions checks, including retained workflow.
manifest=read('report/figures/redesigned/figure-manifest.tsv')
require(len(manifest)==24 and not manifest.figure.duplicated().any(),'24 unique manifest entries')
for row in manifest.itertuples():
 for field in ['source_tables','generator','plot_data']:
  value=getattr(row,field)
  if pd.isna(value):continue
  require(all((ROOT/x).is_file() for x in value.split(';')),row.figure+' existing '+field)
pdfs=sorted(OUT.glob('*.pdf'));require(len(pdfs)==24,'24 complete figure sets')
for p in pdfs:
 require(p.with_suffix('.svg').is_file() and p.with_suffix('.png').is_file(),p.stem+' exports')
 with fitz.open(p) as doc:
  require(len(doc)==1,p.stem+' single canvas')
  page=doc[0]; im=Image.open(p.with_suffix('.png'))
  require(abs(im.width-page.rect.width*5)<2 and abs(im.height-page.rect.height*5)<2 and all(abs(x-360)<.1 for x in im.info['dpi']),p.stem+' 360 dpi dimensions')
  require(not page.get_images(),p.stem+' vector PDF geometry')
  text=page.get_text()
  require(len(text)>100,p.stem+' searchable vector text')
  # No non-whitespace glyph lies outside the canvas.
  outside=[ch['c'] for block in page.get_text('rawdict')['blocks'] if 'lines' in block for line in block['lines'] for span in line['spans'] for ch in span['chars'] if ch['c'].strip() and not page.rect.contains(fitz.Rect(ch['bbox']))]
  require(not outside,p.stem+' canvas text bounds')
with fitz.open(ROOT/'report/report.pdf') as doc:
 reporttext='\n'.join(p.get_text() for p in doc)
 require('Zilong Wang' in reporttext and 'zilongwang0826-ops' in reporttext,'Author and GitHub retained')
 require('## 8.' not in (ROOT/'report/report.md').read_text(),'Section 8 remains removed')
 require(all(not p.get_images() for p in doc),'Report embeds vector scientific figures')
 pages=len(doc)
 require(pages>0,'Report contains pages',pages=pages)
 # Use the actual placement transform and searchable PDF fonts, not nominal
 # source theme sizes. Notes and mathematical subscripts are recorded separately.
 placement_record=json.loads((DOC/'report-placements.json').read_text())
 placements=placement_record['placements']
 require(len(placements)==4,'Four main scientific figures embedded')
 note_prefixes=('Each point is','PCA uses','Primary selection:','Adjusted P =',
                'Six biological themes','FDR uses','Across-day FDR',
                'Correlation and direction','The dashed line','it is not genome')
 typography=[]
 for placement in placements:
  scale=placement['scale_x']
  require(abs(scale-placement['scale_y'])<1e-12,'Isotropic report placement '+placement['source_pdf'])
  report_page=doc[placement['report_page']-1]
  report_spans=[span for block in report_page.get_text('dict')['blocks'] for line in block.get('lines',[]) for span in line['spans'] if span['text'].strip()]
  entries=[]
  with fitz.open(ROOT/placement['source_pdf']) as source:
   for block in source[0].get_text('dict')['blocks']:
    for line in block.get('lines',[]):
     line_text=''.join(span['text'] for span in line['spans']).strip()
     for span in line['spans']:
      text=span['text'].strip()
      if not text: continue
      role='note' if line_text.replace('\u2212','-').startswith(note_prefixes) else 'essential'
      if text=='10' and 'log' in line_text and 'BH' in line_text: role='mathematical_subscript'
      effective=span['size']*scale
      # Confirm that the predicted size occurs in the actual final report text.
      require(any(other['text'].strip()==text and abs(other['size']-effective)<.02 for other in report_spans),'Final report font matches source transform',figure=Path(placement['source_pdf']).stem,text=text)
      if role!='mathematical_subscript':
       minimum=7.0 if role=='note' else 7.5
       require(effective>=minimum-0.01,'Final-size '+role+' typography',figure=Path(placement['source_pdf']).stem,text=text,effective_pt=effective,minimum_pt=minimum)
      entries.append(dict(text=text,role=role,source_pt=span['size'],effective_pt=effective))
  typography.append(dict(figure=Path(placement['source_pdf']).stem,report_page=placement['report_page'],scale=scale,minimum_essential_pt=min(x['effective_pt'] for x in entries if x['role']=='essential'),minimum_note_pt=min(x['effective_pt'] for x in entries if x['role']=='note'),text=entries))
 (DOC/'typography.json').write_text(json.dumps(dict(status='PASS',method='Source PDF span size multiplied by actual placement scale, cross-checked against final report spans.',mathematical_subscripts='Relative mathematical typesetting; base axis label subject to essential-label threshold.',figures=typography),indent=2)+'\n')
output=dict(status='PASS',baseline_tables=len(baseline),original_assets_unchanged=len(oldfig),figure_sets=len(pdfs),report_pages=pages,checks=checks,
 limitations=['Libraries are observed units; pooling and one-mouse-to-one-library independence unresolved.','No raw processing, model fitting, hypothesis testing or inference rerun.','Narrative GO selection is illustrative; full rankings and nonsignificant downward results retained.'])
(DOC/'validation.json').write_text(json.dumps(output,indent=2)+'\n')
print(json.dumps({k:v for k,v in output.items() if k not in ['checks','limitations']},indent=2));print('Checks:',len(checks))
