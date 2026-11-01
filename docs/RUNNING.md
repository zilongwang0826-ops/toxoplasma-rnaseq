# Running the workflow

## Read the report and inspect the results

The [research report](../report/report.md) and [PDF](../report/report.pdf) provide the complete narrative. Rebuild the PDF using the instructions in [report/README.md](../report/README.md).

For diagnostic details, open `example-results/report/index.html` or `example-results/qc/multiqc_report.html` locally. Preserve the entire `example-results/` tree so figures and tables resolve. These are historical outputs, not results of the release checks. The MultiQC file retains its upstream interface and software attribution.

## Local checks

From the repository root with Python 3.9 or later:

```bash
./run.sh test
./run.sh check
python3 scripts/check_example_results.py
```

The first command runs synthetic engineering tests, including explicit environment-lock comparisons and stubbed preparation control flow. The stubs do not run Docker or mamba and do not establish Linux integration. The second rebuilds curated metadata from included source records and checks sample identity. The third verifies the published snapshot: table arithmetic, selected-list consistency, local diagnostic-report references and all 231 copied-file hashes. These checks make no network request. New output goes to `results/`, leaving `example-results/` unchanged. None of these commands runs the RNA-seq analysis.

## Restore the production environment

Use a Linux x86_64 host with Docker Engine, mamba, about 32 GB RAM and 200 GB free storage. The package locks target this architecture. Prepare an isolated project environment:

```bash
# If mamba is not on PATH, set MAMBA_BIN to its executable path first.
bash scripts/prepare_remote.sh
export SNAKEMAKE="$PWD/.workflow-env/bin/snakemake"
./run.sh dry-run
./run.sh dag > workflow.dot
```

Preparation installs the project-local orchestration environment, builds pinned Docker recipes and runs the synthetic container workflow. It requires network access and substantial storage. It does not download the full sequencing dataset. `containers/build.sh` prefers resolved package locks and records new local image IDs; image IDs can differ between equivalent rebuilds. Keep the resulting locks and test evidence. Existing historical validation markers cannot validate a changed release: run `./run.sh container-test` after any relevant source or image change.

Before container building, preparation exports the installed orchestration environment to `results/environment/workflow-linux-64.explicit.txt` and compares it with `config/locks/workflow-linux-64.explicit.txt`. This check applies both to newly created and reused environments. It compares complete package URLs and hashes, ignoring comments, blank lines, surrounding whitespace and entry order. Missing, extra or changed packages fail; empty or malformed explicit lists also fail. Merely finding a working Snakemake executable is insufficient.

On a mismatch, preparation stops before building containers or running their tests. It preserves the existing environment and the export. Inspect the differences first, copy the export and diagnostics to a dated evidence directory, and determine whether the installed environment or the intended lock is wrong. If restoring the committed lock, retain the existing environment under a unique backup directory before deliberately creating a replacement from that lock; rerun preparation to verify it. Do not overwrite the committed lock with the observed export just to pass the check. An intentional dependency change requires a separately reviewed lock update and renewed validation.

The original YAML fallback remains available when the explicit lock is absent or empty, with an explicit **unlocked, not lock-verified** warning. This path does not establish conformity with a pinned environment. A distributed locked release should retain its committed explicit lock. To inspect an export without preparing containers:

```bash
python3 scripts/check_workflow_environment.py \
  config/locks/workflow-linux-64.explicit.txt \
  results/environment/workflow-linux-64.explicit.txt
```

The archived run used Snakemake 9.27.1; Salmon 1.10.3; R 4.4.2; DESeq2 1.46.0; tximport 1.34.0; org.Mm.eg.db and GO.db 3.20.0. See the explicit locks and runtime records for complete packages. Do not install production scientific dependencies into a global Mac environment.

## Execute a new analysis

Review `config/config.yaml` and resource capacity before starting. The default result directory is `results/`, independent of the included `example-results/`. Initialize or clone a Git repository, commit the source/configuration and environment locks, then run:

```bash
export SNAKEMAKE="$PWD/.workflow-env/bin/snakemake"
./run.sh run
```

The driver requires clean tracked source, current container validation and accessible Docker. Each scientific rule runs in the pinned container; the profile limits scheduling to two cores and one job. This downloads approximately 65.493 GiB of FASTQ plus references. Source reads and references are protected outputs; preserve them on interrupted downloads and follow [download recovery](download-recovery.md). Avoid force flags that replace validated inputs.

A complete integrity audit needs the full workspace, including reads, quantifications, fitted models and annotation snapshot. The compact example alone cannot pass that audit:

```bash
python3 workflow/scripts/container_exec.py --image tools --threads 2 --memory-mb 1000 -- \
  python scripts/audit_completed_run.py --project . --results results \
  --output results/audits/new-verification --rehash-raw
```

Choose a new output directory each time. Compare results by stable feature IDs and declared thresholds, not only hit counts. For a later rerun in the same workspace, an unchanged `./run.sh run` should report no work; that expectation was not re-tested against the full data during release preparation.
