# Statistical analysis specification

## Experimental unit and estimands

The primary experimental unit is the reported biological replicate, represented by one archived brain RNA-seq library. The design contains 12 libraries: three infected and three control libraries at each of 11 and 33 days post-infection. These are cross-sectional groups; no pairing or longitudinal reuse of a mouse is assumed. The article's reference to pooled replicates remains unresolved despite the separate sample records. Consequently, independence of biological replicates is an explicit assumption requiring author clarification. Batch information is unavailable and no unsupported batch covariate is introduced.

Primary contrasts are infected versus control **within each time point**, fitted separately with `~ condition`. This permits comparison with the publication without imposing a common temporal effect. An optional secondary model, `~ day + condition + day:condition`, estimates `(infected33 - control33) - (infected11 - control11)` on the log2 scale. It is exploratory, has limited power, and does not describe within-animal change.

## Distinct analysis units and count import

The published supplementary tables contain 936 acute and 2,081 chronic **transcript** hits, corresponding to 701 and 1,448 unique annotated genes. These quantities are not interchangeable.

Two separate analysis tracks use the same verified Salmon files and Ensembl release 100 mapping:

- **Gene track:** `tximport(type="salmon", tx2gene=..., countsFromAbundance="no")`, followed by `DESeqDataSetFromTximport`. DESeq2 receives the supported sample-specific, abundance-weighted transcript-length offsets. This track supports the primary gene-level functional interpretation.
- **Transcript track:** `tximport(type="salmon", txOut=TRUE, countsFromAbundance="lengthScaledTPM")`, followed by `DESeqDataSetFromTximport`. These library-scaled estimated counts are an explicit tximport transformation, not a TPM matrix supplied as raw counts. This track supports comparison with published transcript identities. It does not propagate Salmon inferential uncertainty and is not a differential transcript usage analysis.

The implementation checks sample identity and order, transcript identity/order across files, finite nonnegative estimates, nonempty libraries, and a unique transcript-to-gene mapping. Trailing numeric Ensembl versions are removed only for mapping/comparison; collisions are rejected. Quantification quality and minimum assigned-fragment checks are enforced upstream by the workflow. Every sample is retained unless an explicit, documented redesign is made; no post-hoc exclusion based on significance is implemented.

## Filtering, model fitting and thresholds

Within each six-library comparison, retain features with at least 10 estimated counts in at least three libraries (configurable). DESeq2 estimates normalization factors, dispersions and unshrunken log2 fold changes; performs two-sided Wald tests of zero effect; applies its default Cook's-distance handling and independent filtering; and calculates Benjamini–Hochberg adjusted P values. Automatic outlier replacement is unavailable at three replicates per condition under DESeq2's default minimum replication requirement.

The primary classification follows the publication's results section: `padj < 0.05` and `abs(log2FoldChange) > 1.5`. The effect-size cutoff is a descriptive filter after a zero-effect test; it is **not** a test that the effect exceeds 1.5. The methods-section definition is retained as a separately labelled exploratory sensitivity column: nominal `pvalue < 0.05` and `abs(log2FoldChange) >= 1.5`. Nominal-P lists have no FDR control and cannot replace the primary findings. Adjustments are separate for each feature level and time point, not study-wide across all analyses.

All input features remain in `results.tsv`. Status values distinguish count-filtered features, Cook's-distance/untestable features, independent-filtered features, and features with defined adjusted P values. Fitted R objects preserve offsets, model coefficients and dispersions. No package version is inferred from the paper's invalid `DEseq2 v.3.6.1` specification; container and runtime versions are recorded separately.

## Diagnostic figures

PCA and sample distances use variance-stabilized values from count-filtered features; PCA uses up to the 500 most variable features. The variance-stabilizing transformation uses `blind=FALSE` with the declared design, but selection does not depend on significance. Heatmaps show up to 40 variable features as row Z scores; hierarchical clustering determines row and column order, and a shared symmetric color scale is shown. They are exploratory diagnostics, not independent evidence of treatment effects. MA plots highlight FDR status; volcano plots additionally display the prespecified effect-size cutoff. All plots are written as vector PDF and 300 dpi PNG, and the underlying plotting values are retained where appropriate.

Bulk brain differences may reflect cell-composition changes, immune infiltration, or within-cell regulation. These data do not distinguish those mechanisms. Modest sample size, unresolved pooling, unreported batch effects and BALB/c versus reference-strain differences constrain interpretation.

## Functional enrichment

Only the **gene-level** primary lists feed functional analysis. GO annotations are extracted from the analysis container's `org.Mm.eg.db` and `GO.db`; database metadata, package versions and the exact association snapshot are saved. `GOALL` includes ancestor terms. Ensembl IDs remain the counting unit, and duplicate gene-term association rows are removed. This prevents repeated associations but does not prove that distinct Ensembl IDs represent distinct biological loci.

For each day and direction, the eligible universe consists of genes with nonmissing BH-adjusted P values in that contrast, intersected with genes annotated in the relevant database. This matches genes that could have entered the selected list. Mapping loss is reported for both universe and selected lists. Gene sets with 10–500 eligible genes are tested using a one-sided hypergeometric test. Zero-overlap eligible terms remain in the multiple-testing family. BH adjustment spans BP, CC and MF together within each time point and direction; families across time/direction are separate. GO redundancy and expression/length-related selection bias remain limitations of over-representation analysis.

KEGG runs only when a frozen, appropriately licensed `term_id, term_name, gene_id` TSV snapshot is provided. Its checksum is recorded and it is corrected as a separate family. Without a snapshot, the workflow writes an explicit not-executed status; it does not substitute an uncontrolled live database response. No enriched term is evidence of causality or biochemical pathway activation.

### Overlapping-locus counting sensitivity

The documented overlapping chr10 records `ENSMUSG00000062593` and `ENSMUSG00000112148` share the modern Lilrb4a mapping. Both are FDR-eligible at both time points; only the second is significant at 11 dpi, and both are significantly upregulated at 33 dpi. `locus_sensitivity.R` preserves the primary per-Ensembl results and maps only this pair to one counting unit. Membership in each selected direction is the union of its originally selected members. Annotation associations are also unioned and deduplicated; equality of the two frozen GO term sets is checked and recorded.

All four time/direction families are recalculated using their actual annotated universe, the same 10–500 term-size limits and a complete BH correction including zero-overlap terms. The script first independently reproduces the primary ORA probabilities from the frozen snapshot. It reports changes in background/selected counts, eligible term families and significant terms, including size-boundary changes. This is a sensitivity to identifier/locus counting, **not** a new locus-level differential-expression test or locus-level FDR guarantee. Other many-to-one Entrez groups are audited but are not automatically merged.

### Saved-model and presentation checks

`model_diagnostics.R` reads the saved DESeq2 objects without refitting. It verifies sample and feature identity, design rank, length offsets or size factors, coefficient convergence, dispersion validity and Cook's-distance diagnostics. Flagged features and per-sample diagnostics remain inspectable; they do not automatically remove samples or replace counts. With only three libraries per group, a high Cook statistic does not by itself determine exclusion: the pinned DESeq2 result procedure also has documented count-based exceptions.

`publication_figures.R` renders directly labelled PCA and scaled sample-distance figures from the original saved plotting tables. It performs no normalization, PCA refit, model change or sample exclusion. The primary figures and numeric coordinate/distance tables remain available alongside the presentation versions.

## Published-results comparison

`compare_published.R` joins executed transcript-level estimates to the verified Supplementary Table 1 using transcript IDs. It records availability, testability, FDR eligibility, selected-list overlap, Jaccard index, direction concordance and Spearman effect-size correlation. Some supplementary gene-symbol cells contain numeric values consistent with spreadsheet conversion; comparison therefore uses stable Ensembl IDs and does not repair symbols speculatively. The published table contains significant hits only, so absent transcripts are not known negatives and specificity cannot be estimated. Correlation on the published selected set is selection-biased and descriptive. Nominal P values were not supplied by the paper.

The comparison does not assert a successful reproduction from equal hit counts or visual similarity. Computational assessment requires acceptable read/quantification QC, verified library identities and contrast assignments, documented model-execution checks, effects and pathways consistent with the evidence, and an explanation of residual differences. These checks do not establish mouse-level independence: biological inference remains conditional on the unresolved pooling assumption stated above. Exact numerical equality is not a prespecified success requirement given software, preprocessing and threshold ambiguities.

## Execution and validation boundaries

`tests/test_analysis.R` uses deterministic simulated negative-binomial counts with known positive and negative effects and a time interaction. It tests design and sample-order guards, count filtering, effect direction, the DESeq2 tximport offset interface, interaction orientation, mathematical enrichment probabilities, BH correction, empty-result handling, and figure creation. Its test outputs are explicitly synthetic and excluded from scientific interpretation. It can exercise the supported DESeq2 import interface using constructed tximport objects when the `tximport` package itself is unavailable; this does not validate Salmon file import or the complete pipeline.

Real import, GO database lookup and container integration require the declared container dependencies. Passing local synthetic tests alone does not meet the project completion criteria.

## Sources

- Hu et al. (2020), [original publication](https://doi.org/10.3389/fmicb.2020.570903), methods, Table 1, Figure 2 and Supplementary Table 1.
- [DESeq2 Bioconductor reference and vignette](https://bioconductor.org/packages/DESeq2/).
- [tximport Bioconductor reference and vignette](https://bioconductor.org/packages/tximport/).
