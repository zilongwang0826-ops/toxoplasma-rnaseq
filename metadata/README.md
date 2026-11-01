# Dataset and publication metadata

`sample_id` uses the paper's C/T + day + replicate notation, cross-checked against the archive's group/replicate labels and exact clean-read counts. `run_accession`, `biosample`, `experiment_accession`, FASTQ URLs, MD5 digests and compressed sizes come from the ENA run report. Condition, day and replicate are independently verified against NCBI BioSample attributes. FASTQ URLs use HTTPS equivalents of the ENA-provided FTP host and path.

The twelve libraries are archived as three biological replicates in each of four condition/time groups. The paper's pooling language remains unresolved. `replicate` is not a subject identifier and does not establish pairing across groups. Unknown batch and strandedness are explicitly recorded as `unknown`. `read_pairs`, `base_count` and `read_length` describe archive metadata; full-file checks occur during workflow execution.

## Files

- `samples.tsv`: executable sample manifest, including both mates and checksums.
- `archive_paper_read_counts.tsv`: archive-versus-publication read-count checks.
- `validation.json`: metadata validation outcome and input SHA-256 hashes.
- `published_dets.tsv`: normalized published significant **transcript** list, for comparison only.
- `published_table1_acute.tsv`, `published_table1_chronic.tsv`: source-column-preserving workbook extracts.
- `published_dets_validation.json`: extracted reference counts, effect bounds and overlaps.
- `source/`: immutable small public source responses and original supplementary workbooks.
- `source/sources.json`: retrieval URLs, acquisition timestamps, SHA-256 hashes, and supplemental licensing information.
- `build_metadata.py`, `normalize_supplement.py`, `verify_sources.py`: reproducible source-evidence maintenance utilities.

No full sequencing reads are stored here. Source records do not represent an independently executed RNA-seq analysis. Do not feed the selected published DET table into DESeq2 or enrichment as if it were newly generated results.

## Supplementary-table schema

`published_dets.tsv` contains `level`, `feature_id`, `day`, `log2FoldChange`, `pvalue`, `padj`, `direction`, `gene_id`, `gene_symbol`, `gene_symbol_status`, `source_sheet`, and `source_row`. All entries have `level=transcript`; the paper reports ENSMUST identifiers without version suffixes. `direction` is `up` or `down` relative to infected versus control. `pvalue=NA` means no nominal P value was provided. The source's one exact zero adjusted P value is retained. `source_row` is the one-based Excel row, including the header row. Gene associations are published annotation, not newly inferred mappings.

## Attribution

Hu R-S, He J-J, Elsheikha HM, Zou Y, Ehsan M, Ma Q-N, Zhu X-Q, Cong W (2020). Transcriptomic Profiling of Mouse Brain During Acute and Chronic Infections by *Toxoplasma gondii* Oocysts. *Frontiers in Microbiology* 11:570903. [Article DOI](https://doi.org/10.3389/fmicb.2020.570903). The paper and supplementary tables are licensed CC BY 4.0. The [publisher's Figshare collection](https://api.figshare.com/v2/collections/5178827) supplies the original supplementary files. Derived metadata tables preserve attribution and identify their transformations.

Five published symbol cells contain numeric values (one acute and four chronic). These values are preserved as strings and flagged `gene_symbol_status=numeric_in_source`; they may reflect spreadsheet date conversion, but that cause is unconfirmed. They must not be used as valid gene symbols. Ensembl identifiers remain the comparison keys; current labels should come from the pinned reference annotation.
