# Reproduction targets and acceptance evidence

The publication is [Hu et al. (2020)](https://doi.org/10.3389/fmicb.2020.570903). Its [Supplementary Table 1](https://doi.org/10.3389/fmicb.2020.570903.s002) is the reference list for feature-level comparison. The table contains Ensembl **transcript** identifiers, corresponding gene identifiers, log2 fold changes, direction labels and adjusted P values. It contains selected significant transcripts only, with no raw counts, full tested universe or nominal P values.

## Published targets verified from supplementary data

| Measure | 11 dpi | 33 dpi |
| --- | ---: | ---: |
| Published differentially expressed transcripts | 936 | 2,081 |
| Upregulated transcripts | 789 | 1,908 |
| Downregulated transcripts | 147 | 173 |
| Distinct associated Ensembl genes | 701 | 1,448 |
| Smallest absolute published log2 fold change | 1.501389933 | 1.5000992 |
| Largest published adjusted P value | 0.0498975495963185 | 0.049931754 |

There are 618 transcripts upregulated in both comparisons and six downregulated in both. These counts were computed from the downloaded supplementary table and agree with the paper's totals and overlap statements. They validate extraction and interpretation of the **published reference**, not independent reproduction of the RNA-seq analysis.

One chronic-infection adjusted P value is stored as zero in the source spreadsheet. The unmodified value is preserved. Any negative-log visualization must use an explicitly reported finite floor for plotting only. Nominal P values are unavailable and represented as `NA`, never inferred from adjusted P values.

The source has five numeric cells in the Gene Symbol column (one acute, four chronic), potentially resulting from spreadsheet date conversion. This explanation is unconfirmed. The extract preserves these values and flags them as `numeric_in_source`; matching uses Ensembl identifiers and display labels should be obtained from the pinned reference annotation. No gene symbol has been silently repaired.

## Comparison definitions

| Target | Required inputs | Evaluation | Interpretation limit |
| --- | --- | --- | --- |
| Dataset identity | Twelve paired libraries, archive checksums and sample assignment evidence | Verify all MD5 sums, mate integrity and sample mappings before quantification | Archive metadata validation does not validate downloaded FASTQ contents |
| Transcript DET results | Full Salmon quantification and transcript-level DESeq2 output per day | Count significant up/down transcripts using BH adjusted P < 0.05 and absolute log2FC > 1.5; report count differences from the table above | Matching counts alone do not establish reproduction |
| Feature recovery | Version-normalized Ensembl transcript IDs from both result sets | Report intersection, union, Jaccard similarity and published-list recovery; include denominator and unmatched IDs | Annotation changes and absent historical transcript versions can reduce overlap |
| Effect direction and magnitude | Matched transcript log2FC estimates | Report direction agreement and effect correlation for matched features with finite effects, including sample size | The published list is selected on significance; its correlations are not a genome-wide replication estimate |
| Gene-level response | tximport gene aggregation and gene-level DESeq2 output | Analyze as a complementary modern interpretation; compare mapped published gene sets descriptively | The paper's 936/2,081 transcript counts are not gene-level targets; multiple isoforms can map to one gene |
| Sample structure | Full library quantification and transformed counts | Examine unsupervised PCA, correlations and technical QC; report unanticipated structure | PCA separation is descriptive; bulk tissue cell mixtures and unknown batch factors remain possible explanations |
| Enrichment themes | Gene results, versioned GO annotations and tested gene universe | Evaluate interferon/antigen-presentation themes at 11 dpi and lymphocyte/T-cell activation themes at 33 dpi with corrected enrichment tests | Do not substitute expected themes for computed results; modern annotation/background can change individual term rankings |
| Sensitivity to reporting gaps | Identical inputs under prespecified alternative thresholds or options | Separate nominal-P from BH-FDR DET summaries; assess declared methodological deviations | Threshold sensitivity does not repair unavailable preprocessing details or uncertain pooling |

Gene-level enrichment should use a gene universe derived from the actual analysis eligibility and annotation coverage. A transcript-derived gene set and a gene-level DE set answer related but different questions; their results must be labeled separately. The historical clusterProfiler/GO/KEGG snapshots and background universe are unavailable from the paper. GO-based biological-theme concordance is therefore more defensible than claiming exact historical pathway reconstruction. KEGG-specific conclusions require a separately recorded and legally accessible KEGG snapshot and mapping policy.

The original study also includes qRT-PCR, transcription-factor networks and predicted host–parasite interactions. Those results are background context. This project's required core is the RNA-seq differential response and functional enrichment. It does not imply new wet-lab validation or claim experimental validation of computationally predicted interactions.

## Success criteria and current evidence

Completion requires all intended real libraries to pass identity checks, a successful containerized end-to-end workflow, valid count import and statistical contrasts, executed enrichment with a documented universe, figures generated from real outputs, and a documented comparison with the published reference. Missing resources or unresolved critical failures must remain visible in the project status. Exact numerical agreement is not required for a scientifically useful reconstruction, but disagreement must be quantified and its plausible causes separated from verified explanations.

The complete real-data analysis has executed. See [results.md](results.md) for the measured differences and [validation-status.md](validation-status.md) for acceptance evidence. `metadata/validation.json` and `metadata/published_dets_validation.json` remain source-metadata checks, distinct from the later real-data execution and model/artefact validation.

## Rebuilding reference tables

The source workbook remains unchanged. From the repository root:

```bash
python3 metadata/build_metadata.py
python3 metadata/normalize_supplement.py
python3 metadata/verify_sources.py
```

These are source-evidence maintenance commands, not RNA-seq analysis steps. They require Python's standard library and the retained small source files. `metadata/published_dets.tsv` normalizes column names and adds analysis level, day and source-row provenance. `metadata/published_table1_acute.tsv` and `metadata/published_table1_chronic.tsv` preserve all source columns. No nominal P values or absent whole-transcriptome results are manufactured.
