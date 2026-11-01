# Ensembl release-100 reference scope

The initial real-data pilot stopped during reference preparation because the Ensembl `cdna.all` and `ncrna` target files contain transcripts from patch/alternate regions absent from the ordinary `Mus_musculus.GRCm38.100.gtf.gz` annotation.

The actual downloaded ordinary GTF contained 142,699 transcript IDs across 45 regions. Of the selected FASTA targets, 1,686 cDNA and 354 ncRNA IDs were absent from **all** GTF feature rows, not merely its transcript rows. The first rejected ID, `ENSMUST00000180575`, is located on `CHR_MG153_PATCH` in the FASTA header. This identifies an annotation-scope mismatch rather than a transcript-row parser omission.

The reconstruction now selects the same assembly/release's comprehensive annotation, [Mus_musculus.GRCm38.100.chr_patch_hapl_scaff.gtf.gz](https://ftp.ensembl.org/pub/release-100/gtf/mus_musculus/Mus_musculus.GRCm38.100.chr_patch_hapl_scaff.gtf.gz). The original GTF and acquisition receipt are retained under a timestamped `resources/ensembl100/source-history/` directory. No FASTA targets are silently dropped. Reference preparation also checks FASTA-header versus GTF gene identity when a header gene field is present.

Target FASTA scope includes patch and alternate loci; the decoy genome remains the selected primary assembly. Alternative transcript representations may increase ambiguous mappings. This reconstructed file set is explicit and must not be claimed to be the exact unpublished reference selection used by Hu et al.

`scripts/audit_reference_scope.py` records target coverage, gene identity, input hashes and differences between annotation scopes. The comprehensive GTF contains 144,778 transcript IDs across 118 regions. The measured audit covered all 141,450 unique FASTA targets, with zero missing transcript IDs and zero FASTA/GTF gene-ID mismatches. Reference preparation then completed with 53,028 gene IDs and 66 primary-assembly decoy sequences. Compact audit and checksum evidence is in `docs/validation/reference-*-scope.json` and `real-reference-prepared.json`.

## Completed index target accounting

The completed sparse index name table contains all 141,450 target transcript IDs and all 66 decoy names, with no missing or unexpected target identifiers. All 3,017 published transcript-by-day entries have IDs represented in that table.

Forty-six targets have original lengths of 9–30 nt and cannot enter a k=31 graph, although Salmon retains their names and lengths in index metadata. None occurs in the published DET lists. `first_decoy_index=141404` uses an internal compact identifier; applying the stored short-reference correction restores the full name-table decoy start to 141450. It must not be interpreted as loss of 46 metadata records. Evidence and the exact IDs are in `validation/index-target-audit.json`.

Salmon reports clipping poly-A tails from 747 derived-reference sequences and substituting 78,088,277 non-ATCG nucleotides; original downloaded sequences remain unchanged. `keep_duplicates=true` retains duplicate target identities. These are recorded quantifier preprocessing operations, not silent project-level exclusions.

Source interpretation follows the pinned pufferfish [FixFasta short-reference handling](https://github.com/COMBINE-lab/pufferfish/blob/salmon-v1.10.3/src/FixFasta.cpp#L414-L458), [reference-name restoration](https://github.com/COMBINE-lab/pufferfish/blob/salmon-v1.10.3/src/PufferfishBinaryGFAReader.cpp#L327-L384), and [sparse loader](https://github.com/COMBINE-lab/pufferfish/blob/salmon-v1.10.3/src/PufferfishSparseIndex.cpp#L183). The completed C111 quant.sf retains all 46 short targets with zero TPM and NumReads; this was checked directly. Metadata and single-library technical consistency do not establish valid biological conclusions.
