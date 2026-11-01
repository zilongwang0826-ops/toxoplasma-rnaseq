# Historical execution and numerical checks — 29 October 2026

This file summarizes the archived 12-library execution and subsequent computational checks. The retained records report successful execution, result/model consistency checks, report-link checks and an unchanged no-work invocation. These are historical records, not new executions during release preparation or independent biological validation. The release-specific checks and remaining verification limits are listed in [RELEASE_CHECKS.md](RELEASE_CHECKS.md).

| Check | Executed evidence | Coverage and limits |
| --- | --- | --- |
| Archive/sample metadata | 12 runs, 12 BioSamples, four groups of three; source documents and supplementary lists verified | Pooling, batch and sequencer discrepancies remain documented |
| Raw FASTQ integrity | All 24 files independently reread for MD5 and SHA-256; exact archive bytes/MD5 and receipt SHA-256 matched | Raw content was verified at audit time; subsequent external changes are not monitored |
| Read quality | 24 FastQC record counts match archive pairs; base-quality, N-content and adapter modules pass | Composition/GC/duplication flags retained and assessed; not automatic grounds for trimming |
| Quantification | All 12 libraries; exact processed-pair counts, finite estimates, target identity/order, count/TPM totals, common transcript/decoy hashes and IU layout | Salmon 1.10.3 reconstruction; no inferential replicate uncertainty propagation |
| Python tests | 36 tests pass locally and in the approved tools container | Includes transport recovery, corrupt inputs, protected evidence, delivery links and no-overwrite audit behavior |
| Container smoke workflow | 11 jobs pass, two threads and one job; includes Salmon, DESeq2, tximport, full GO integration and locus-counting checks | Synthetic data verify implementation, not biological conclusions |
| Four primary DE contrasts | Independent result/filter/list consistency, sample order, normalized counts, all 14 input MD5s; BH/Wald/baseMean recalculations agree within numerical tolerance | Three reported independent libraries per group, conditional on pooling ambiguity |
| Six saved models | Exact design/colData/features, full rank, positive finite length offsets or size factors, all coefficients converged, valid dispersions | No sample removal, refit or count replacement; Cook flags remain inspectable |
| GO and locus sensitivity | Four complete tested-background families; frozen annotation snapshot; documented overlapping pair collapsed separately | Original tables retained; broad selected themes stable, some small terms sensitive |
| Published transcript comparison | All published identifiers represented; executed selected-list overlap, direction and effect comparisons | Published selected lists are not known genome-wide positives/negatives; correlations are selection-conditioned |
| Primary-result preservation | SHA-256 of 69 primary TSV files identical before and after supplementary validation | Presentation files and final provenance/report are separately regenerated |
| Full delivery audit | PASS: 24 raw files, 12 quantifications, 24 FastQC records, required outputs, input hashes, report links and stage metadata | The original study exit precedes supplementary jobs; their own completion records are retained |
| Presentation and no-work resume | PASS: corrected GO labels inspected; representative figures of every type visually reviewed; final artifact/link audit passed; actual `./run.sh run` reports no work | Primary 69 TSV hashes remain unchanged; no browser-wide pixel regression or independent second-machine run claimed |
| Independent environment rebuild | Not executed | Validated images ran successfully; resolved locks/recipes are supplied but not independently rebuilt on a second machine |

## Actual environments and provenance

Production used Linux x86_64, Snakemake 9.27.1 in its isolated host environment and the two content-addressed Docker images recorded in `containers/images.lock.json`. Scientific jobs enforce two CPUs and bounded memory; indexing remained below the approved 24 decimal GB cap after documented recovery. Runtime statistics: R 4.4.2, DESeq2 1.46.0, tximport 1.34.0, AnnotationDbi 1.68.0, org.Mm.eg.db 3.20.0 and GO.db 3.20.0. Salmon 1.10.3, FastQC 0.12.1 and MultiQC 1.27.1 were executed.

The end-of-report Git revision is not assigned retroactively to every earlier stage. Launch records and per-output Snakemake timestamps/commands are retained; the audit explicitly labels stage-to-launch revision matching as temporal association rather than a per-job Git attestation. Public compact evidence is under `docs/validation/`; detailed operational/development records remain outside this repository.

The primary study exited with status 0 after report generation. A separate six-job validation pass produced model diagnostics, locus sensitivity and labelled QC figures, then regenerated provenance/report. Primary input/result tables were not recomputed. Source changes to presentation require their own validated rerun and final artifact check.

## Consequential review findings

T331 host assignment is 73.46%, below its infected-day 33 comparators, but its PCA/distance and Cook diagnostics do not justify removal. T111 has a weaker acute response; nearest-neighbour statements are limited to the stated analysis scope. All 12 samples remain in the declared design.

The overlapping Lilrb4a/Gm49339 identifiers affect individual GO terms. The frozen annotation sets are identical (467 terms). After collapsing this pair, day 11 upward significance remains 929 terms with one lost/one gained; day 33 changes 1129→1116. Eleven day 33 significant terms leave the minimum-size family; two retained terms lose significance. Downward families remain empty. See `enrichment-identifier-audit.md` for the complete interpretation and limitations.

Four primary models contain 13/34/69/125 rows without P values after count filtering. Saved max-Cook flags do not correspond one-to-one with P-value removal: the installed DESeq2 1.46.0 two-level-design code retains a flagged row when at least three counts exceed the count in the sample with maximum Cook distance. The 28 retained-P exceptions are all nonsignificant at BH 0.05. A few day 33 P values underflow to numerical zero; plots apply a finite lower bound.

## Recovery and scope limits

The initial ordinary GTF omitted 2,040 selected targets; comprehensive same-release annotation resolved the mismatch without filtering targets. Initial index and quantifier loading attempts experienced confirmed OOM and were retried within approved limits. Three interrupted HTTPS acquisitions were resumed from preserved partial bytes and passed archive validation. These were real diagnosed failures, not evidence-free restart attempts. See `reference-scope.md`, `resource-validation.md` and `download-recovery.md`.

The approved dataset is complete; no KEGG snapshot, qPCR, predicted host–parasite network validation or second full independent run is claimed. Statistical inference remains conditional on biological replication, unknown batches and bulk-tissue interpretation limits.
