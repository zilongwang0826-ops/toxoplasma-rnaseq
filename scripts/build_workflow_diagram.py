#!/usr/bin/env python3
"""Draw the executed RNA-seq workflow as SVG and embeddable vector artwork.

The report renderer embeds the same Drawing directly, preserving vector text
and arrows. This is a methods schematic, not a new scientific analysis.
"""
from pathlib import Path

from reportlab.graphics import renderSVG
from reportlab.graphics.shapes import Drawing, Line, Polygon, Rect, String
from reportlab.lib.colors import HexColor, white


WIDTH, HEIGHT = 507, 288
INK = HexColor('#173342')
EDGE = HexColor('#78919f')
FILL = HexColor('#edf2f4')


def build_drawing():
    drawing = Drawing(WIDTH, HEIGHT)

    def label(x, y, text, size=8.4, bold=False, color=INK):
        drawing.add(String(x, y, text, textAnchor='middle', fontSize=size,
                           fontName='Helvetica-Bold' if bold else 'Helvetica',
                           fillColor=color))

    def box(x, y, width, height, title, detail, emphasized=False):
        drawing.add(Rect(x, y, width, height, rx=3, ry=3,
                         fillColor=INK if emphasized else FILL,
                         strokeColor=INK if emphasized else EDGE,
                         strokeWidth=0.65))
        color = white if emphasized else INK
        label(x + width / 2, y + height / 2 + 3, title, bold=True, color=color)
        label(x + width / 2, y + height / 2 - 9, detail, size=7.5, color=color)

    def arrow(points):
        for a, b in zip(points, points[1:]):
            drawing.add(Line(*a, *b, strokeColor=EDGE, strokeWidth=0.9))
        x, y = points[-1]
        px, py = points[-2]
        if y < py:
            tip = [x, y, x-2.5, y+5, x+2.5, y+5]
        elif x > px:
            tip = [x, y, x-5, y-2.5, x-5, y+2.5]
        else:
            tip = [x, y, x+5, y-2.5, x+5, y+2.5]
        drawing.add(Polygon(tip, strokeColor=EDGE, fillColor=EDGE, strokeWidth=0.5))

    label(WIDTH / 2, 279, 'Snakemake orchestration  |  Locked Docker images', 9, True)
    box(3, 225, 238, 35, 'ENA paired FASTQ', 'Archive checks and FastQC')
    box(266, 225, 238, 35, 'Ensembl 100 / GRCm38.p6', 'Transcripts, annotation and genome decoys')
    box(134, 172, 239, 35, 'Salmon quantification + QC', 'Decoy-aware index; bias correction', True)
    arrow([(122,225),(122,216),(204,216),(204,207)])
    arrow([(385,225),(385,216),(303,216),(303,207)])
    box(134, 121, 239, 33, 'tximport', 'Separate gene and transcript estimates')
    arrow([(253.5,172),(253.5,154)])
    box(3, 65, 238, 38, 'Gene DESeq2 models', 'Counts + length offsets; 11 and 33 dpi')
    box(266, 65, 238, 38, 'Transcript DESeq2 models', 'lengthScaledTPM counts; 11 and 33 dpi')
    arrow([(253.5,121),(253.5,113),(122,113),(122,103)])
    arrow([(253.5,113),(385,113),(385,103)])
    box(3, 10, 238, 38, 'GO enrichment', 'Tested universe; overlapping-locus sensitivity')
    box(266, 10, 238, 38, 'Published transcript comparison', 'Selected-list recovery and effect agreement')
    arrow([(122,65),(122,48)])
    arrow([(385,65),(385,48)])
    return drawing


def main():
    folder = Path(__file__).resolve().parents[1] / 'report/figures'
    folder.mkdir(parents=True, exist_ok=True)
    renderSVG.drawToFile(build_drawing(), str(folder / 'analysis_workflow.svg'))


if __name__ == '__main__':
    main()
