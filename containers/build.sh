#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
: "${BASE_IMAGE:?Set BASE_IMAGE to mambaorg/micromamba@sha256:<verified manifest digest>}"
case "$BASE_IMAGE" in mambaorg/micromamba@sha256:*) ;; *) echo 'A digest-pinned official micromamba base is required' >&2; exit 2;; esac
mkdir -p results/container-validation containers/locks
printf "image\trecipe\tspec\n" > results/container-validation/build-inputs.tsv
for image in tools statistics; do
    spec="containers/$image.yaml"
    recipe="containers/Dockerfile.$image"
    if [[ ${RESOLVE_NEW:-0} != 1 && -s "containers/locks/$image-linux-64.explicit.txt" ]]; then
        spec="containers/locks/$image-linux-64.explicit.txt"
        recipe="containers/Dockerfile.locked"
    fi
    printf "%s\t%s\t%s\n" "$image" "$recipe" "$spec" >> results/container-validation/build-inputs.tsv
    docker build --build-arg "ENVIRONMENT_FILE=$spec" --build-arg "LOCK_FILE=$spec" --platform linux/amd64 --build-arg "BASE_IMAGE=$BASE_IMAGE" \
        --file "$recipe" --tag "projet2-$image:2026.1" . \
        2>&1 | tee "results/container-validation/build-$image.log"
    resolved="results/container-validation/$image-resolved.explicit.txt"
    docker run --rm --pull=never "projet2-$image:2026.1" micromamba env export --name base --explicit > "$resolved"
    if [[ "$recipe" == containers/Dockerfile.locked ]]; then
        python3 - "$spec" "$resolved" <<'CHECK'
import sys
from pathlib import Path
def packages(path):
    return {line for line in Path(path).read_text().splitlines() if line.startswith('https://')}
if packages(sys.argv[1]) != packages(sys.argv[2]):
    raise SystemExit('Rebuilt environment differs from the explicit input lock')
CHECK
    else
        cp "$resolved" "containers/locks/$image-linux-64.explicit.txt"
    fi
done
python3 - <<'PY'
import csv,hashlib,json,subprocess,os
from pathlib import Path
from datetime import datetime,timezone
images={}
for row in csv.DictReader(open('results/container-validation/build-inputs.tsv'),delimiter='\t'):
    name=row['image']
    tag='projet2-'+name+':2026.1'
    info=json.loads(subprocess.check_output(['docker','image','inspect',tag],text=True))[0]
    files=[Path(row['recipe']),Path(row['spec'])]
    images[name]={'tag':tag,'id':info['Id'],'architecture':info['Architecture'],
                  'recipe_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}}
Path('containers/images.lock.json').write_text(json.dumps({'built_at':datetime.now(timezone.utc).isoformat(),'base_image':os.environ['BASE_IMAGE'],'images':images},indent=2)+'\n')
PY
