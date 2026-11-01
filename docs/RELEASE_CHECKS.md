# Verification scope

## Current prose and layout revision — 1 November 2026

The report now uses black serif text and compact three-line tables. A two-metric comparison table was converted to prose with all four values retained. Repeated qualifications were reduced and the confirmed workflow-building, virtual-machine execution and paper-comparison tasks were integrated into the project aim. The report flows to six pages without reducing figure size; body, table and caption text are larger than in the preceding copy. All six rendered pages were inspected.

The saved-result checker passed again with 231 snapshot hashes and four primary tables containing 388,956 rows. The unchanged figure/report validator passed all 791 checks in a disposable copy with the two SVG attachments restored: 106 baseline tables, 121 original assets, 24 figure sets and four main report figures. Scientific files and figure assets were preserved; no analysis was rerun.

## Earlier local editorial revision — 1 November 2026

The report and README were revised for a clearer account of the M2 reanalysis. The scientific tables, figure assets, analysis code, parameters and environment locks were preserved. The PDF was rebuilt from the revised Markdown using the existing vector figures; all eight rendered pages were inspected.

The saved-result checker passed with 231 snapshot hashes, four primary tables containing 388,956 rows, four GO-summary rows, four comparison-summary rows and 36 local HTML resources. In a disposable copy with the two SVG attachments restored, the figure/report validator passed 791 checks across 106 baseline tables, 121 original assets and 24 figure sets. Only its arbitrary eight-to-ten-page constraint was removed; figure count, typography, geometry and scientific checks remain. No full analysis, quantification, model fit or enrichment calculation was rerun.

## Earlier application package checks — 1 November 2026

The following checks completed in a disposable copy with two compute threads. The two separately packaged SVGs were restored before the full figure check. No validator output overwrote the delivered scientific snapshot.

| Check | Result |
|---|---|
| `python3 -B -m unittest discover -s tests -p 'test_*.py' -v` | 54 tests passed; input checks, synthetic failure/recovery fixtures and stubbed environment-control tests. |
| `python3 -B scripts/check_example_results.py` | Passed: four primary tables containing 388,956 rows, four GO-summary rows, four comparison-summary rows, 36 local HTML resources and 231 snapshot hashes. |
| `python -B scripts/figures/validate_redesign.py` | 791 checks passed across 106 baseline scientific/metadata tables, 121 original assets and all 24 restored figure sets; the report remains eight pages. |

At that packaging checkpoint, the report Markdown/PDF, source records, scientific tables, plotted values and figure files were byte-identical to the accepted package. The later editorial revision is described above. The [source map](RELEASE_SOURCE_MAP.tsv) records the current file hashes; [current selection](PUBLIC_FILE_SELECTION.tsv) identifies the ordinary files and separate vector attachment. [Figure validation](../FIGURE_SCIENTIFIC_VALIDATION.md) explains the protected numerical and display invariants.

## Execution boundaries

These are local engineering and saved-result checks. They do not constitute a new FASTQ-to-result run, container rebuild, Salmon quantification, DESeq2 refit, GO recomputation or independent biological replication. Diagnostic HTML resource paths were checked; interactive browser behavior was not tested. The earlier packaging checks did not rebuild the report or figures. The editorial revision rebuilt and visually checked the report; the figure assets remain unchanged, with their original rendering records retained.

Use [the running guide](RUNNING.md) for a fresh production run and [figure build instructions](../scripts/figures/README.md) for presentation regeneration. The full saved-figure validator requires restoration of the [two supplementary SVGs](RELEASE_ASSETS.md).
