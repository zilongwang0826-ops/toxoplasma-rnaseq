# Source provenance and verification

The project combines a public-data workflow with the compact snapshot of its executed results in `example-results/`. New workflow runs write to the separate, Git-ignored `results/` directory. Scientific inputs, result tables and figure assets are preserved. Local editorial and layout revisions change report text, its renderer/PDF and reader-facing documentation; they do not rerun or alter the analysis. Current file hashes are recorded below.

## Trace the evidence

- [Source metadata](../metadata/README.md) records accessions, retrieval records, checksums and the published comparison lists.
- [Supplementary computational information](../report/supplementary_reproduction.md) describes the completed analysis and its environment.
- [Execution checks](validation-status.md) retain the scientific and computational scope of the original run.
- [Source mapping](RELEASE_SOURCE_MAP.tsv) retains original source identities and hashes while recording current distributed bytes. The map excludes itself and the root checksum index to avoid circular hashes.
- [Current file selection](PUBLIC_FILE_SELECTION.tsv) lists shipped project files and the two unchanged SVGs supplied in the separate [vector attachment](RELEASE_ASSETS.md).

The repository-level `CHECKSUMS.sha256` verifies the current delivered files. Original source hashes and timestamps refer to their original execution, rather than implying a new raw-data run. Internal editorial records, obsolete presentation copies and historical packaging manifests are retained in the local development archive outside the application package.

## Reproduce and maintain

Follow [the running guide](RUNNING.md) for environment restoration, source acquisition and bounded production execution. Preserve raw data, reference identities, software locks and scientific outputs. Run the appropriate [checks](RELEASE_CHECKS.md) after relevant edits. Rebuild changed reports and inspect their final pages; full-data execution, saved-result consistency and independent biological confirmation remain distinct claims.
