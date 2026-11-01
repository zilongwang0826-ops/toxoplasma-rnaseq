# Research report

Project period: October 2026. Historical execution and report-preparation dates are distinguished in the [project timeline](../docs/PROJECT_TIMELINE.md).

Read the [complete report](report.md) on GitHub or download the [PDF](report.pdf). The report presents the study design, methods, results, comparison with the publication, interpretation and references. Scientific tables and source figures are reused without numerical changes.

The [HTML diagnostic appendix](../example-results/report/index.html) and [MultiQC report](../example-results/qc/multiqc_report.html) provide additional plots and QC information after downloading the repository. Keep their relative directory structure intact.

## Supplementary vector exports

Two large transcript-level SVGs are supplied as a separate [vector attachment](../docs/RELEASE_ASSETS.md). Their PDF/PNG versions, source tables and generation scripts remain in this directory tree. Reading or rebuilding the main report does not require those two SVGs. Restore them before running the full saved-figure validator described below; its complete-format checks are unchanged.

## Rebuild the PDF

The renderer uses local Markdown, four redesigned vector scientific figures and a retained vector workflow schematic, built-in PDF fonts, ReportLab, Pillow and pypdf. Body text and headings use black serif type; tables use top, header and bottom rules, content-sized columns and Markdown column alignment. The schematic is generated from `scripts/build_workflow_diagram.py` as SVG and embedded directly as vector artwork in the report. Its two branches preserve the gene and transcript analysis units. The scientific PDF artwork is merged directly without rasterization. The report runtime is recorded in `../scripts/figures/environment-versions.json`; dependencies are pinned in `requirements.txt`.

```bash
python3.12 -m venv .report-venv
.report-venv/bin/python -m pip install -r report/requirements.txt
.report-venv/bin/python scripts/build_workflow_diagram.py
.report-venv/bin/python scripts/build_research_report.py
```

No scientific pipeline or external image request is run by the renderer. After changing report text, rebuild the PDF and inspect every rendered page. `<!-- pagebreak -->` markers in the Markdown control print pagination and are hidden in GitHub's reading view.

The report builder records actual figure placement scales in `../docs/figure-redesign/report-placements.json`. After rebuilding, run `python scripts/figures/validate_redesign.py` from the repository root to check final-size labels against the source PDF and report text. The target is at least 7.5 pt for essential figure text and 7 pt for supporting figure notes. Page count follows the content; retain all four main figures and inspect every rendered page for readable, complete placement.

## Rebuild redesigned figures

See [the frozen-result figure layer](../scripts/figures/README.md), [figure inventory](figures/redesigned/figure-manifest.tsv) and [validation](../FIGURE_SCIENTIFIC_VALIDATION.md). Run `bash scripts/figures/build_all.sh` with the documented R/Python dependencies. This rebuilds all 24 figure sets and the report from saved downstream tables.
