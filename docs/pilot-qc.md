# Real-data pilot QC observations

Status: C111 technical pilot passed, and the subsequent full cohort and scientific validation completed. This page preserves pilot-stage observations; see [results.md](results.md) and [validation-status.md](validation-status.md) for the complete-study assessment.

## C111, read 2

The complete archive file was acquired on 27 October 2026 and verified against ENA MD5 `6d63ddb6087225ea7807b452cd2b1f95` and the expected 3,359,193,762 bytes. FastQC reports 37,526,096 sequences, each 150 bases, matching the archive's expected paired-fragment count; GC content is 49%.

Per-base quality, per-sequence quality, N content, read-length distribution, overrepresented sequences and adapter content passed. Per-base sequence content failed, with the imbalance confined mainly to positions 1–9; subsequent positions are close to balanced. GC distribution and duplication produced warnings. FastQC estimated a deduplicated fraction of 55.259%; this is not an estimate of PCR artefacts.

The start-position pattern is consistent with known RNA-seq priming/fragmentation bias described in the [FastQC sequence-content documentation](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/4%20Per%20Base%20Sequence%20Content.html). [FastQC duplication documentation](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/8%20Duplicate%20Sequences.html) explains that expression abundance can generate repeated RNA-seq reads and that raw sequence identity does not distinguish technical from biological duplicates.

Decision: preserve the archived reads without extra trimming or deduplication. This is a technical QC decision, conditional on the companion mate, Salmon mapping/assignment and cohort-level checks. It does not establish sample suitability by itself or validate the study's biological conclusions. At this initial single-mate check, full pilot and cohort QC were pending; their subsequent completion is recorded below and in the final validation record.

Evidence stays under `results/qc/raw/`, `data/fastq/*.source.json` and the run logs on the compute server; only small QC artifacts are returned locally.

## Completed paired-library pilot

C111 quantification completed at 06:10:35 UTC on 28 October 2026; the Snakemake job finished at 06:10:39 with exit 0 after 54 minutes 7 seconds. Both mates contain 37,526,096 reads, matching the complete archive. Salmon processed exactly 37,526,096 paired fragments and assigned 31,339,012 (83.5126%). It inferred `IU` (inward, unstranded paired-end); the strand-mapping bias was 0.500625, consistent with balanced orientations. This does not establish the library type of the remaining samples.

Both mates passed per-base quality, per-sequence quality, N content, overrepresented sequences and adapter checks. R1 has the same start-localized composition pattern as R2, with imbalance mainly at positions 1–9 and near-balanced composition thereafter. Composition, GC and duplication flags remain documented rather than hidden. No additional trimming, deduplication or sample exclusion was performed.

The completed quant.sf contains exactly 141,450 distinct target IDs, matching the prepared mapping. Numeric entries were finite and nonnegative. NumReads sums to 31,339,012.007 (rounding difference 0.007 from the mapped count); TPM sums to 1,000,000.00001. All 46 short graph-excluded targets are retained with zero TPM and NumReads. Sequence/name hashes match the completed index. This resolves the short-target output check described in `reference-scope.md`.

The automatic gate passed at 06:11:09 UTC, and the existing driver began the remaining 63 workflow jobs at 06:11:12. The completed index and C111 results were reused. Evidence is in `docs/validation/pilot-C111.json`, `C111-salmon-meta.json`, `C111-library-format.json` and `C111-quant-invariants.json`. The approximately 6.7 MiB quant.sf was returned for local verification but remains excluded from Git; raw FASTQ was not returned.

The pilot alone does not validate infection-associated changes. Subsequent cohort QC, gene/transcript inference, enrichment sensitivity and published-results comparison have now executed and been reviewed; their evidence is recorded separately in the final scientific report.
