# Supplementary vector figures

The two large transcript-level SVGs are distributed as `toxoplasma-supplementary-vectors.zip`. PDF and PNG views for all 24 figures remain in this repository. The SVG bytes and scientific figure contents are unchanged.

| Figure | Included views |
|---|---|
| Transcript MA plots | [PDF](../report/figures/redesigned/S_MA_transcript.pdf) · [PNG](../report/figures/redesigned/S_MA_transcript.png) |
| Transcript differential expression | [PDF](../report/figures/redesigned/S_transcript_differential_expression.pdf) · [PNG](../report/figures/redesigned/S_transcript_differential_expression.png) |

## Download and restore

Download the [SVG archive](https://github.com/zilongwang0826-ops/toxoplasma-rnaseq/releases/download/supplementary-vectors/toxoplasma-supplementary-vectors.zip) and its [SHA-256 sidecar](https://github.com/zilongwang0826-ops/toxoplasma-rnaseq/releases/download/supplementary-vectors/toxoplasma-supplementary-vectors.zip.sha256) from the [supplementary vector release](https://github.com/zilongwang0826-ops/toxoplasma-rnaseq/releases/tag/supplementary-vectors). The [asset record](../RELEASE_ASSETS.json) identifies the archive and both original files.

Save the ZIP and checksum sidecar in the same download directory, then verify the archive:

```bash
shasum -a256 -c toxoplasma-supplementary-vectors.zip.sha256
```

The ZIP SHA-256 is `57459957e4f59e38e6d7a7f7970ec7ed37e513367a23b7a8a688720a67950fb6`. From the **toxoplasma-rnaseq repository root**, restore the SVGs without overwriting existing files:

```bash
unzip -n /path/to/downloads/toxoplasma-supplementary-vectors.zip -d .
shasum -a256 -c RELEASE_ASSETS.sha256
```

The archive paths start with `report/figures/redesigned/`. The two restored SVGs are Git-ignored; retain them locally for the complete figure check. If an existing file has a different hash, preserve it and investigate before replacing it.

## Verification

`shasum -a256 -c CHECKSUMS.sha256` verifies ordinary repository files. `RELEASE_ASSETS.sha256` verifies the two restored SVGs. The full `scripts/figures/validate_redesign.py` check requires all figure formats and should be run only after restoring the attachment. The report and saved-result checks do not require these SVGs. [Current file selection](PUBLIC_FILE_SELECTION.tsv) records both distribution routes.
