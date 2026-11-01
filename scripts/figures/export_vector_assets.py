#!/usr/bin/env python3
"""Export PNG/SVG from the source-rendered PDF geometry; no plot-data changes."""
from pathlib import Path
import pymupdf as fitz
from reportlab.graphics import renderPDF
import sys
import csv
from xml.etree import ElementTree as ET
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from build_workflow_diagram import build_drawing

root=Path(__file__).resolve().parents[2]
out=root/'report/figures/redesigned'
drawing=build_drawing()
drawing.width+=12; drawing.height+=12; drawing.translate(6,6)
renderPDF.drawToFile(drawing,str(out/'F00_analysis_workflow.pdf'))
for path in sorted(out.glob('*.pdf')):
    with fitz.open(path) as document:
        assert len(document)==1
        page=document[0]
        svg=ET.fromstring(page.get_svg_image(text_as_path=True))
        ns='http://www.w3.org/2000/svg'
        ET.register_namespace('',ns)
        svg.insert(0,ET.Element('{'+ns+'}rect',width='100%',height='100%',fill='white'))
        path.with_suffix('.svg').write_text(ET.tostring(svg,encoding='unicode'))
        pixels=page.get_pixmap(matrix=fitz.Matrix(5,5),alpha=False)
        pixels.set_dpi(360,360)
        pixels.save(path.with_suffix('.png'))
print(f'Exported {len(list(out.glob("*.pdf")))} figure sets at 360 dpi')

manifest=out/'figure-manifest.tsv'
with manifest.open() as stream:
    reader=csv.DictReader(stream,delimiter='\t'); fields=reader.fieldnames; rows=[r for r in reader if r['figure']!='F00_analysis_workflow']
rows.append(dict(figure='F00_analysis_workflow',category='methods',width_in=str(drawing.width/72),height_in=str(drawing.height/72),dpi='360',source_tables='report/report.md;workflow/Snakefile',plot_data='',purpose='Retained documented gene/transcript workflow',generator='scripts/build_workflow_diagram.py;scripts/figures/export_vector_assets.py'))
with manifest.open('w',newline='') as stream:
    writer=csv.DictWriter(stream,fieldnames=fields,delimiter='\t');writer.writeheader();writer.writerows(rows)
