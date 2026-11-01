# Study design and metadata evidence

## Scientific question and scope

How does the mouse brain transcriptional response to oral infection with *Toxoplasma gondii* Pru oocysts differ from the contemporaneous uninfected state at 11 and 33 days post-infection (dpi)? This project reconstructs the public RNA-seq analysis of [Hu et al. (2020), Frontiers in Microbiology](https://doi.org/10.3389/fmicb.2020.570903), rather than claiming exact re-execution of unpublished code.

The primary comparisons are infected versus control at each time point. The study is cross-sectional: different mice were killed at each time point. Replicate number is a label within each group, not a pairing or longitudinal subject identifier. The intended experimental unit is one independently sampled mouse brain; the archived unit is one reported biological replicate/library. The pooling ambiguity below prevents independently proving a one-mouse-to-one-library mapping.

## Verified sample organization

[PRJNA483261](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA483261) has 12 runs, 12 experiments and 12 distinct BioSamples. ENA sample aliases and NCBI BioSample `sample_type` and `biological_replicate` attributes agree on all assignments. `metadata/samples.tsv` is the executable sample sheet. `metadata/build_metadata.py` checks these assignments against archived source responses and paper Table 1.

| Time | Condition | Paper sample IDs | Runs in sample order | Replicates |
| --- | --- | --- | --- | --- |
| 11 dpi | Control | C111, C112, C113 | SRR7657877, SRR7657876, SRR7657879 | 3 |
| 11 dpi | Infected | T111, T112, T113 | SRR7657878, SRR7657881, SRR7657880 | 3 |
| 33 dpi | Control | C331, C332, C333 | SRR7657883, SRR7657882, SRR7657875 | 3 |
| 33 dpi | Infected | T331, T332, T333 | SRR7657874, SRR7657873, SRR7657872 | 3 |

The paper sample codes are reconstructed from verified condition, time and replicate labels, and are corroborated by exact agreement with each paper sample's clean-read count. They are not original SRA accession fields. All records describe female BALB/c mouse brain. BioSample age is 6–8 weeks; the paper describes seven-week-old mice. The paper reports oral challenge with 100 Pru oocysts and PBS controls. No sequencing-batch or cage assignments are available in the retrieved records. They must not be invented or introduced as known covariates.

## Archive identity and preprocessing

All twelve runs are paired-end. The archive lists 409,862,759 read pairs and 122,958,827,700 bases; the ratio is 300 bases per pair for every library, consistent with 2 × 150 bp. These are archive metadata, not an independently completed FASTQ inspection.

For **every** sample, twice the ENA read-pair count equals the paper's clean-read count exactly, and differs from its raw-read count. The checks are recorded in `metadata/archive_paper_read_counts.tsv`. This is strong evidence that the archived reads are already filtered. Consequently, primary processing preserves the downloaded files and does not automatically trim them again. Any additional trimming requires an observed QC problem and a documented sensitivity analysis. The exact preprocessing that produced the archived files cannot be recovered from counts alone.

ENA lists 70,322,230,064 bytes of compressed FASTQ across both mates of all libraries (70.322 GB or 65.493 GiB). Full files and their MD5 digests are listed in the sample sheet. No complete FASTQ has been downloaded as part of metadata verification.

## Material uncertainties and contradictions

| Topic | Evidence | Consequence |
| --- | --- | --- |
| Replicate pooling | The paper describes twelve mice in four groups and three biological replicates per group, but also contains the phrase “of 3 pooled replicates” in library preparation. The archive explicitly labels three distinct biological replicates per group. | Proceeding with three independent libraries per group is supported by archive labels and Figure 2, but independence at the mouse level remains conditional. Author clarification would resolve this. |
| Sequencer | Paper: HiSeq 2500. All ENA runs: HiSeq X Ten. | Preserve both records; do not silently choose one as historically certain. |
| Selection/strandedness | Paper describes mRNA libraries with NEBNext Ultra; ENA lists RT-PCR selection. Strandedness is not stated in retrieved records. | Use automatic Salmon library-type inference and retain its evidence; do not infer strandedness from the kit name alone. |
| DESeq2 version | Paper reports DEseq2 v3.6.1. | This does not identify a usable DESeq2 package version and may conflate R and package versions. Pin and report the actual replacement environment. |
| Significance | Methods: unadjusted P < 0.05 and absolute log2 fold change ≥ 1.5. Results: adjusted P < 0.05 and absolute log2 fold change > 1.5. | Primary analysis controls BH FDR; record a nominal-P sensitivity table separately. Never call nominal-P results FDR-controlled. |
| Statistical input | Paper describes TPM normalization before count-based DESeq2, without specifying the import object. | Use supported Salmon count import and length handling. TPM must not be treated as raw counts. |
| Reference | Mouse GRCm38.p6, Ensembl release 100 is reported. Exact FASTA/GTF files, transcript selection, index command and decoys are omitted. | Record exact reconstructed reference URLs, checksums and index settings; assess differences explicitly. |
| Preprocessing | Trimmomatic 0.39 is named, but adapters and trimming options are omitted. Archive counts match clean counts. | Preserve archived data and avoid unmotivated repeat trimming. |
| Supplement symbols | Five supplementary Gene Symbol cells contain numeric values. Ensembl gene/transcript IDs remain present. | Preserve the source and flag the values; use Ensembl IDs for matching and pinned annotation for labels. |
| Enrichment | clusterProfiler 3.6.2 and BH adjusted P < 0.05 are reported. Database snapshot, eligible background, mapping policy and gene-set size bounds are unspecified. | Freeze available annotations and use the actually tested, mappable gene universe. Treat historical pathway identities as imperfectly recoverable. |

## Analysis commitments

Each time-point comparison uses only its three infected and three control libraries, with an intercept and condition term. The sign of every reported log2 fold change is infected relative to control. This avoids confounding infection with age/time and does not impose a shared infection effect across time points. A combined time-by-condition interaction would answer an additional question and must not be substituted silently.

Gene-level analysis sums transcript abundance through an explicit transcript-to-gene map and retains the import method's length correction. Transcript-level results must be labeled separately and are the closest unit to the paper's stated DET analysis. Gene counts and transcript counts are not numerically interchangeable. Multiple-testing correction is applied within each prespecified comparison and analysis level; a claim spanning both comparisons requires its own multiplicity definition.

Sample QC will assess library size, Salmon mapping rate and inferred library type, transformed-count PCA, sample correlations and the relation between technical metrics and condition. No sample is excluded merely to improve significance or reproduce the published DET count. With three libraries per group, uncertainty and outlier sensitivity remain substantial. Bulk tissue effects can reflect changing cell composition as well as regulation within cells.

## Source and licensing

Primary sources are the [publisher article](https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2020.570903/full), [Europe PMC XML](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC7604304/fullTextXML), [NCBI BioProject](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA483261), ENA run metadata and NCBI BioSample XML. Retrieval URLs, UTC times and SHA-256 hashes are recorded in `metadata/source/sources.json`. The paper and its supplementary material are distributed under CC BY 4.0, attributed to Hu and colleagues (2020). Public source snapshots are scientific evidence, not workflow execution results.
