#!/usr/bin/env bash
# Prepare the isolated Linux production environment and validate its containers.
set -euo pipefail
cd "$(dirname "$0")/.."
[[ $(uname -s) == Linux && $(uname -m) == x86_64 ]] || { echo 'Linux x86_64 required' >&2; exit 2; }
MAMBA_BIN=${MAMBA_BIN:-$(command -v mamba || true)}
[[ -n "$MAMBA_BIN" && -x "$MAMBA_BIN" ]] || { echo 'Install mamba or set MAMBA_BIN to its executable path' >&2; exit 2; }
docker info >/dev/null
python3 - <<'PY'
import shutil
if shutil.disk_usage('.').free < 30*1024**3: raise SystemExit('Need 30 GiB free for preparation and reserve')
PY
export CONDA_PKGS_DIRS="$PWD/.cache/conda"
mkdir -p "$CONDA_PKGS_DIRS" results/environment
spec=config/workflow-environment.yaml
lock=config/locks/workflow-linux-64.explicit.txt
if [[ -s "$lock" ]]; then spec="$lock"; fi
if [[ ! -x .workflow-env/bin/snakemake ]]; then
    "$MAMBA_BIN" create --yes --prefix "$PWD/.workflow-env" \
        --strict-channel-priority --file "$spec"
fi
.workflow-env/bin/snakemake --version
resolved=results/environment/workflow-linux-64.explicit.txt
"$MAMBA_BIN" env export --prefix "$PWD/.workflow-env" --explicit > "$resolved"
if [[ "$spec" == "$lock" ]]; then
    python3 scripts/check_workflow_environment.py "$lock" "$resolved"
else
    echo 'WARNING: no nonempty workflow explicit lock; YAML fallback is unlocked and not lock-verified.' >&2
fi
BASE_IMAGE=$(python3 -c 'import json; d=json.load(open("containers/base-image.json")); print(d["repository"]+"@"+d["manifest_digest"])')
export BASE_IMAGE
bash containers/build.sh
SNAKEMAKE="$PWD/.workflow-env/bin/snakemake" ./run.sh container-test
