#!/usr/bin/env python3
"""Build the research PDF from report/report.md and source-driven local vector figures.

Requires reportlab, Pillow and pypdf. This renderer reads only the supported Markdown
blocks used by the report and never fetches remote images or executes HTML.
"""
import argparse
import json
from html import escape
from pathlib import Path
import re

from PIL import Image as PILImage
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.enums import TA_LEFT, TA_RIGHT, TA_CENTER
from reportlab.pdfbase.pdfmetrics import stringWidth
from reportlab.platypus import Flowable, Image, KeepTogether, PageBreak, Paragraph, SimpleDocTemplate, Spacer, Table, TableStyle

from build_workflow_diagram import build_drawing
from pypdf import PdfReader, PdfWriter, Transformation

VECTOR_PLACEMENTS = []

class VectorFigure(Flowable):
    """Record a rectangle for lossless insertion of an existing PDF artwork."""
    def __init__(self, path, width, height):
        super().__init__()
        self.path, self.width, self.height = path, width, height
        self.hAlign = "CENTER"
    def draw(self):
        x, y = self.canv.absolutePosition(0, 0)
        VECTOR_PLACEMENTS.append((self.canv.getPageNumber()-1, self.path, x, y, self.width, self.height))

def insert_vector_artwork(path):
    reader = PdfReader(path)
    writer = PdfWriter()
    for page in reader.pages:
        writer.add_page(page)
    for page_number, figure_path, x, y, width, height in VECTOR_PLACEMENTS:
        artwork = PdfReader(figure_path).pages[0]
        transform = Transformation().scale(width/float(artwork.mediabox.width), height/float(artwork.mediabox.height)).translate(x, y)
        writer.pages[page_number].merge_transformed_page(artwork, transform)
    if reader.metadata:
        writer.add_metadata({str(k): str(v) for k,v in reader.metadata.items()})
    temporary = path.with_suffix(".vector-tmp.pdf")
    with temporary.open("wb") as stream:
        writer.write(stream)
    temporary.replace(path)


def record_vector_placements(root, output):
    """Record the actual affine scale for final-size typography verification."""
    placements = []
    for page, source, x, y, width, height in VECTOR_PLACEMENTS:
        artwork = PdfReader(source).pages[0]
        source_width = float(artwork.mediabox.width)
        source_height = float(artwork.mediabox.height)
        placements.append({
            'report_page': page + 1,
            'source_pdf': str(source.relative_to(root)),
            'x_pt': x, 'y_pt': y,
            'width_pt': width, 'height_pt': height,
            'source_width_pt': source_width, 'source_height_pt': source_height,
            'scale_x': width / source_width, 'scale_y': height / source_height,
        })
    record = root / 'docs/figure-redesign/report-placements.json'
    record.write_text(json.dumps({
        'report': str(output.relative_to(root)) if output.is_relative_to(root) else output.name,
        'units': 'PDF points; y measured from page bottom',
        'placements': placements,
    }, indent=2) + '\n')


INK = colors.black
MUTED = colors.black
LINE = colors.black
WIDTH = A4[0] - 88
STYLES = {
    'title': ParagraphStyle('title', fontName='Times-Bold', fontSize=20, leading=23, textColor=INK, spaceAfter=12),
    'heading': ParagraphStyle('heading', fontName='Times-Bold', fontSize=13, leading=16, textColor=INK, spaceBefore=10, spaceAfter=7, keepWithNext=True),
    'subheading': ParagraphStyle('subheading', fontName='Times-Bold', fontSize=11, leading=14, textColor=INK, spaceBefore=7, spaceAfter=5, keepWithNext=True),
    'body': ParagraphStyle('body', fontName='Times-Roman', fontSize=11, leading=14.1, textColor=INK, spaceAfter=8),
    'caption': ParagraphStyle('caption', fontName='Times-Roman', fontSize=9.5, leading=12, textColor=MUTED, spaceAfter=9),
    'cell': ParagraphStyle('cell', fontName='Times-Roman', fontSize=10, leading=12.1, textColor=INK),
    'header': ParagraphStyle('header', fontName='Times-Bold', fontSize=10, leading=12.1, textColor=INK),
}


def inline(text):
    text = escape(text)
    text = re.sub(r'`([^`]+)`', r'<font name="Courier" size="8.5">\1</font>', text)
    text = re.sub(r'\*\*([^*]+)\*\*', r'<b>\1</b>', text)
    text = re.sub(r'(?<!\*)\*([^*]+)\*(?!\*)', r'<i>\1</i>', text)
    text = re.sub(r'\[([^\]]+)\]\(([^)]+)\)', r'<link href="\2" color="#000000">\1</link>', text)
    return text


def make_table(lines):
    rows = [[cell.strip() for cell in line.strip().strip('|').split('|')] for line in lines]
    separator = next((row for row in rows if all(re.fullmatch(r':?-+:?', cell) for cell in row)), None)
    rows = [row for row in rows if not all(re.fullmatch(r':?-+:?', cell) for cell in row)]
    count = len(rows[0])
    # Fit columns to their actual text rather than stretching sparse tables.
    widths = [max(stringWidth(row[j], 'Times-Bold' if i == 0 else 'Times-Roman', 10)
                  for i, row in enumerate(rows)) + 14 for j in range(count)]
    if sum(widths) > WIDTH:
        widths = [width * WIDTH / sum(widths) for width in widths]
    alignments = [TA_CENTER if cell.startswith(':') and cell.endswith(':') else
                  TA_RIGHT if cell.endswith(':') else TA_LEFT
                  for cell in (separator or ['---'] * count)]
    cells = []
    for i, row in enumerate(rows):
        cells.append([Paragraph(inline(cell), ParagraphStyle(
            f'table-{i}-{j}', parent=STYLES['header' if i == 0 else 'cell'],
            alignment=alignments[j])) for j, cell in enumerate(row)])
    result = Table(cells, colWidths=widths, repeatRows=1, hAlign='LEFT')
    result.setStyle(TableStyle([
        ('LINEABOVE', (0, 0), (-1, 0), .7, LINE),
        ('LINEBELOW', (0, 0), (-1, 0), .4, LINE),
        ('LINEBELOW', (0, -1), (-1, -1), .7, LINE),
        ('VALIGN', (0, 0), (-1, -1), 'TOP'),
        ('LEFTPADDING', (0, 0), (-1, -1), 7), ('RIGHTPADDING', (0, 0), (-1, -1), 7),
        ('TOPPADDING', (0, 0), (-1, -1), 4), ('BOTTOMPADDING', (0, 0), (-1, -1), 4),
    ]))
    return result


def read_blocks(source):
    lines = source.read_text(encoding='utf-8').splitlines()
    blocks, index = [], 0
    while index < len(lines):
        line = lines[index].strip()
        if not line:
            index += 1
            continue
        if line == '<!-- pagebreak -->':
            blocks.append(PageBreak()); index += 1; continue
        if line.startswith('|'):
            selected = []
            while index < len(lines) and lines[index].strip().startswith('|'):
                selected.append(lines[index]); index += 1
            blocks.extend([KeepTogether([make_table(selected)]), Spacer(1, 11)]); continue
        if line.startswith('#'):
            match = re.fullmatch(r'(#{1,3})\s+(.+)', line)
            if not match:
                raise ValueError('Unsupported heading: ' + line)
            name = {1: 'title', 2: 'heading', 3: 'subheading'}[len(match[1])]
            blocks.append(Paragraph(inline(match[2]), STYLES[name])); index += 1; continue
        picture = re.fullmatch(r'!\[([^\]]*)\]\(([^)]+)\)', line)
        if picture:
            path = (source.parent / picture[2]).resolve()
            if not path.is_relative_to(source.parent.parent.resolve()):
                raise ValueError('Figures must be local repository files')
            if path == (source.parent / 'figures/analysis_workflow.svg').resolve():
                figure = build_drawing()
                figure.hAlign = 'CENTER'
                group = [figure, Spacer(1, 5)]
                index += 1
                while index < len(lines) and not lines[index].strip():
                    index += 1
                if index < len(lines) and lines[index].startswith('**Workflow overview.'):
                    group.append(Paragraph(inline(lines[index]), STYLES['caption']))
                    index += 1
                blocks.append(KeepTogether(group))
                continue
            # Markdown keeps portable PNG previews; print uses sibling vector art.
            if path.parent.name == 'redesigned' and path.with_suffix('.pdf').is_file():
                path = path.with_suffix('.pdf')
            if path.suffix == '.pdf':
                page = PdfReader(path).pages[0]
                w, h = float(page.mediabox.width), float(page.mediabox.height)
            else:
                with PILImage.open(path) as image:
                    w, h = image.size
            scale = min(WIDTH / w, 390 / h)
            figure = VectorFigure(path, w * scale, h * scale) if path.suffix == '.pdf' else Image(str(path), width=w * scale, height=h * scale)
            figure.hAlign = 'CENTER'
            group = [figure, Spacer(1, 5)]
            index += 1
            while index < len(lines) and not lines[index].strip():
                index += 1
            if index < len(lines) and lines[index].startswith('**Figure '):
                group.append(Paragraph(inline(lines[index]), STYLES['caption']))
                index += 1
            blocks.append(KeepTogether(group))
            continue
        paragraph = [line]
        index += 1
        while index < len(lines) and lines[index].strip() and not lines[index].strip().startswith(('#', '|', '<!--', '![')):
            paragraph.append(lines[index].strip()); index += 1
        text = ' '.join(paragraph)
        style = STYLES['caption'] if text.startswith('**Figure ') else STYLES['body']
        blocks.append(Paragraph(inline(text), style))
    return blocks


def frame(canvas, doc):
    canvas.saveState()
    canvas.setFont('Times-Roman', 10)
    canvas.setFillColor(INK)
    canvas.drawCentredString(A4[0] / 2, 26, str(doc.page))
    canvas.restoreState()


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=root / 'report/report.md')
    parser.add_argument('--output', type=Path, default=root / 'report/report.pdf')
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    title = next((line[2:].strip() for line in args.source.read_text(encoding='utf-8').splitlines()
                  if line.startswith('# ')), args.source.stem)
    doc = SimpleDocTemplate(str(args.output), pagesize=A4, rightMargin=44, leftMargin=44,
                            topMargin=45, bottomMargin=43, title=title,
                            author='Zilong Wang', subject='PRJNA483261 computational RNA-seq reconstruction')
    VECTOR_PLACEMENTS.clear()
    doc.build(read_blocks(args.source.resolve()), onFirstPage=frame, onLaterPages=frame)
    insert_vector_artwork(args.output)
    record_vector_placements(root, args.output.resolve())
    print(args.output)


if __name__ == '__main__':
    main()
