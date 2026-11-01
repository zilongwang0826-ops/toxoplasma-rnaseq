# Reproducibility-gap register

| Gap | Evidence and impact | Implemented handling | Resolution status |
| --- | --- | --- | --- |
| Biological pooling | Paper wording differs from 12-library / biological-replicate description | Three libraries per group; mouse independence is conditional | Author clarification unavailable |
| Sequencer | Paper HiSeq 2500; ENA HiSeq X Ten | Preserve both claims | Unresolved |
| Archived preprocessing | Every archive count equals the published clean count | No automatic retrimming; retain archive MD5 | Exact original trimming cannot be reconstructed |
| Historical reference file set | GRCm38.p6 / Ensembl 100 reported without exact file list | Freeze cDNA, ncRNA, primary-assembly genome and comprehensive chr_patch_hapl_scaff GTF URLs/hashes; retain the original GTF scope audit | Selected reconstruction, not proven identical |
| Historical quantifier | Salmon 1.2.0 reported; missing full command | Salmon 1.10.3 with explicit selective alignment/decoys/bias options | Modernized method |
| DESeq2 version | Reported 3.6.1 does not identify the package environment | Pin DESeq2 1.46.0 in a Bioconductor 3.20 container | Historical environment unresolved |
| Count import | TPM mentioned before count-based testing | Supported tximport gene offsets and transcript abundance-derived counts | Reconstruction explicitly differs where original is unclear |
| Significance thresholds | Methods nominal P and ≥1.5; results adjusted P and >1.5 | BH primary; nominal-P sensitivity retained separately | Both definitions documented |
| Statistical design/batches | Per-comparison models and batch metadata omitted | Day-specific contrasts; no invented covariates or pairing | Unknown confounding cannot be corrected |
| Transcript/gene unit | Supplementary table explicitly lists transcript and gene IDs | Separate inference and comparison by transcript ID | Unit resolved; published numerical symbol cells flagged |
| GO/KEGG background | Annotation date, universe and mapping policy omitted | Versioned GO snapshot; tested/mappable genes; controlled duplicate mapping | Historical enrichment not exactly recoverable |
| Transcript estimation uncertainty | No documented inferential replicates | Report count-model inference without uncertainty propagation | Limitation; no isoform-usage claims |
| Complete workflow execution | No published executable workflow or exact package lock | Independent implementation with Git, content-addressed images and execution provenance | Executed on all 12 libraries; independent numerical/model review and locus sensitivity completed |
| Host-parasite networks and qPCR | Additional experimental/database analyses in original paper | Outside selected RNA-seq reconstruction targets | Not claimed as reproduced |

See `study-design.md` and `reproduction-targets.md` for source-linked evidence and prespecified comparison criteria. No implementation decision is tuned to force the published hit counts.
