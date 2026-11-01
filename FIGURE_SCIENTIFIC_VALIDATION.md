# Toxoplasma scientific validation

## Phase 2 validation — 31 October 2026

PASS for saved-result and presentation checks. [Recorded validation](docs/figure-redesign/validation.json) records 791 passed checks, including the original numerical checks and final-size typography checks. All 106 baseline scientific/metadata/audit TSV tables match their SHA-256 baseline hashes. The release source-map TSV is packaging metadata maintained separately and is excluded from the scientific hash check. All 121 original standalone figure assets remain byte-identical. No raw processing, normalization, PCA recomputation, model fitting, enrichment testing or new inference was run.

The unit remains the reported biological library: 12 libraries and three per time/condition group. Different animals at 11 and 33 dpi are not treated as paired. Reference assembly GRCm38.p6, Ensembl release 100, annotation provenance, Salmon and DESeq2 methods remain as documented in the preserved analysis snapshot.

## Verified invariants

| Evidence | Check |
|---|---|
| All PCA views | Coordinate values and sample IDs aligned by key; time and condition joined to metadata; saved explained-variance labels retained; equal axis units. |
| All distances | Every matrix cell retained and compared by row/column sample ID; symmetric input with zero diagonal; complete-linkage ordering changes only position. |
| All heatmaps | Every value in each saved 40-feature matrix matched by feature/sample key; feature IDs and row Z scores retained. |
| Volcanoes | Every finite FDR-eligible estimate retained; `-log10(max(padj, machine-positive floor))` verified; missing adjusted P is not interpreted as nonsignificance. |
| Primary DE selection | Strict `padj < 0.05` and `abs(log2FoldChange) > 1.5` unchanged. Gene up/down counts: 496/41 at day 11, 1,173/18 at day 33. Transcript counts: 724/130 and 1,804/101. |
| MA plots | Mean normalized count shown on a log10 axis; all finite effects retained; original BH-only color criterion explicitly distinguished from the volcano effect filter. |
| GO | Stable term IDs align enrichment ratio, overlap, selected/universe sizes, raw P and adjusted P to source; top-15 selection follows padj/pvalue/ID; all nonsignificant downward families retained. |
| Main GO comparison | Same six narrative themes at both times; full names/IDs; exact counts and source adjusted P. Selection is illustrative and not a new test or claim of independent pathways. |
| Published agreement | 760 and 1,743 finite pairs; conditional Spearman 0.878021867474803 and 0.932812261968955 independently reproduced from matched effects; direction fractions preserved. |
| Recovery | Shared selected counts 696/936 and 1,624/2,081 independently checked using the full published-list denominator; rounded display 74.36% and 78.04%. |
| Exports/report | 24 complete SVG/PDF/PNG sets, 360-dpi dimensions, existing input paths, vector PDF geometry/text, canvas bounds and retained attribution. Essential main-figure labels are at least 7.5 pt after placement in the report; supporting notes are at least 7 pt. |

Numerical alignment uses explicit identity keys, not implicit row order, with relative tolerance 1e-12 and absolute tolerance 1e-14 for numeric serialization. Hash identity protects the unmodified full tables. Display-only changes are point/label placement, term/feature ordering, matrix reshaping, descriptive counts, log-axis transforms and finite underflow floors. No confidence intervals or new statistical uncertainty were invented.

## Tests actually run

- `Rscript scripts/figures/redesign_figures.R`: all 23 scientific compositions generated successfully from saved tables; no warnings in the final build log.
- `python scripts/figures/export_vector_assets.py`: 24 complete format sets, including the retained workflow.
- `python scripts/build_research_report.py`: eight-page report generated with vector artwork.
- `python scripts/figures/validate_redesign.py`: 791 checks passed, including final report placement and typography checks.
- `./run.sh test`: 54 tests passed, including 14 environment-lock/control-flow tests. Transfer-failure messages are deliberate synthetic recovery fixtures, not failed tests. Docker and mamba operations in the new control-flow tests are stubs.
- `python scripts/check_example_results.py`: PASS; four full primary result tables, four GO-summary rows, four comparison-summary rows, 36 local diagnostic-report resources and 231 snapshot hashes checked.

The inference-oriented R suites were not rerun because the statistical model code and scientific tables were unchanged. Existing documented inference validation remains historical evidence. This verification concerns saved-result arithmetic, provenance and presentation; it is not an independent reproduction of the raw-data analysis. See [release checks](docs/RELEASE_CHECKS.md) for current execution boundaries.

## Unresolved scientific limitations preserved

Pooling language leaves one-mouse-to-one-library independence unresolved. Batch/cage assignments are unavailable. Whole-brain expression mixes cell populations; lymphocyte-associated patterns can reflect infiltration, within-cell regulation or both. Conditional publication-list correlation is not genome-wide accuracy. Across-day GO adjusted-P differences are not tests of pathway-activity differences. Parent/child GO redundancy remains explicit in full rankings. No new scientific error was discovered or silently corrected during this presentation upgrade.
