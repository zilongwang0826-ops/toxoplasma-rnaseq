#!/usr/bin/env bash
# Orchestration runs on the host; every scientific job enters Docker explicitly.
set -euo pipefail
cd "$(dirname "$0")"
MODE=${1:-help}
if (($#)); then shift; fi
SMK=${SNAKEMAKE:-snakemake}
PROFILE=${PROFILE:-profiles/local}
CONFIG=${CONFIG:-config/config.yaml}
mkdir -p .cache
SOURCE_CACHE=$(mktemp -d "$PWD/.cache/snakemake.XXXXXX")
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
case "$MODE" in
    help)
        cat <<'EOF'
Usage: ./run.sh {check|dry-run|dag|test|container-test|run} [Snakemake options]
  check           Validate curated metadata locally (no download).
  dry-run         Validate the full DAG without Docker or data downloads.
  dag             Print the workflow DAG as DOT.
  test            Run lightweight Python unit tests locally.
  container-test  Run synthetic indexing/quantification and statistics in Docker.
  run             Execute the complete containerized study after resource approval.
Set CONFIG, PROFILE and SNAKEMAKE to select configurations/executable.
EOF
        ;;
    check)
        python3 metadata/build_metadata.py
        python3 workflow/scripts/validate_inputs.py --samples metadata/samples.tsv --output results/validation/local-metadata.json
        ;;
    test)
        python3 -m unittest discover -s tests -p 'test_*.py' -v
        ;;
    dry-run|dag)
        extra=(--dry-run)
        if [[ "$MODE" == dag ]]; then extra=(--dag); fi
        "$SMK" --runtime-source-cache-path "$SOURCE_CACHE" --snakefile workflow/Snakefile --configfile "$CONFIG" --profile "$PROFILE" "${extra[@]}" "$@"
        ;;
    container-test)
        docker info >/dev/null
        "$SMK" --runtime-source-cache-path "$SOURCE_CACHE" --snakefile tests/smoke.smk --cores 2 --jobs 1 --scheduler greedy --printshellcmds --forceall "$@"
        python3 scripts/record_container_validation.py
        ;;
    run)
        docker info >/dev/null
        test -f containers/images.lock.json || { echo 'Build and validate pinned images first' >&2; exit 2; }
        python3 scripts/record_container_validation.py --check
        git rev-parse HEAD >/dev/null
        if [[ -n "$(git status --porcelain --untracked-files=no)" ]]; then
            echo 'Commit tracked source/configuration changes before a scientific run' >&2; exit 2
        fi
        mkdir -p results/launch
        stamp=$(date -u +%Y%m%dT%H%M%SZ)
        { git rev-parse HEAD; printf 'CONFIG=%s\nPROFILE=%s\n' "$CONFIG" "$PROFILE"; printf '%q ' "$0" run "$@"; printf '\n'; } > "results/launch/$stamp.txt"
        "$SMK" --runtime-source-cache-path "$SOURCE_CACHE" --snakefile workflow/Snakefile --configfile "$CONFIG" --profile "$PROFILE" "$@" \
            2>&1 | tee "results/launch/$stamp.log"
        ;;
    *) echo "Unknown command: $MODE" >&2; exit 2;;
esac
