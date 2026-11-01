#!/usr/bin/env bash
# Continue an already approved real-data phase after the existing pilot finishes.
# This script does not grant approval or launch/restart a missing pilot.
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p results/launch
exec 9>results/launch/study-driver.lock
flock -n 9 || { echo 'Another study driver already owns this workspace' >&2; exit 2; }
trap 'status=$?; printf "%s\n" "$status" > results/launch/study.exit' EXIT
printf '%s\n' "$$" > results/launch/study-driver.pid
export SNAKEMAKE="$PWD/.workflow-env/bin/snakemake"
while [[ ! -f results/launch/pilot.exit ]]; do
    if ! pgrep -f '[s]nakemake .*results/salmon/C111/quant.sf$' >/dev/null; then
        sleep 5
        [[ -f results/launch/pilot.exit ]] || { echo 'Pilot process missing without completion status; inspect before resuming' >&2; exit 3; }
    fi
    sleep 30
done
[[ $(cat results/launch/pilot.exit) == 0 ]] || { echo 'Pilot failed; full study not started' >&2; exit 4; }
python3 workflow/scripts/container_exec.py --image tools --threads 2 --memory-mb 1000 -- \
    python scripts/pilot_gate.py --sample C111
python3 - <<'PY'
import csv,json,shutil
from pathlib import Path
rows=list(csv.DictReader(open('metadata/samples.tsv'),delimiter='\t'))
archive=sum(int(r[f'bytes_{m}']) for r in rows for m in (1,2))
receipts=list(Path("resources/ensembl100").rglob("*.source.json"))
if not receipts:raise SystemExit("Reference acquisition receipts are missing")
reference_bytes=sum(json.loads(path.read_text())["bytes"] for path in receipts)
if archive+reference_bytes>75*1024**3:raise SystemExit('Approved 75 GiB download budget exceeded')
remaining=sum(int(r[f'bytes_{m}']) for r in rows for m in (1,2) if not Path(f"data/fastq/{r['sample_id']}_{m}.fastq.gz").exists())
if shutil.disk_usage('.').free<remaining+35*1024**3:
    raise SystemExit('Insufficient storage for remaining compressed reads, 10 GiB derived outputs and 25 GiB reserve')
PY
printf 'PILOT GATE PASSED; continuing the approved full study\n'
./run.sh run
