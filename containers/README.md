# Container specifications and validation

Two Linux amd64 images separate read processing from R analysis. The first build uses `tools.yaml` and `statistics.yaml` with pinned application versions and strict channel priority. The official micromamba base's registry manifest digest is recorded in `base-image.json`; image layers were not pulled during the audit.

On the selected Linux production host:

```bash
export BASE_IMAGE=$(python3 -c 'import json; d=json.load(open("containers/base-image.json")); print(d["repository"]+"@"+d["manifest_digest"])')
bash containers/build.sh
./run.sh container-test
```

The build writes manager-generated explicit package specifications and a JSON image lock. The execution wrapper requires a local `sha256:` image ID and matching recipe hashes; Snakemake also records the image ID in rule parameters so an image change is not invisible to its rerun logic. A separate successful-test marker is bound to the image lock and relevant source hashes. A failed build or test is not sufficient to enable a scientific run.

The original image build retains a diagnostic package list internally; executable locks are generated with `micromamba env export --explicit`, including `@EXPLICIT` and package hashes. For an exact package rebuild, `build.sh` uses `Dockerfile.locked` and prefers an existing nonempty `containers/locks/<image>-linux-64.explicit.txt` over a YAML solve. To deliberately resolve new dependencies, preserve the prior locks, update the specifications and select `RESOLVE_NEW=1` after review. Docker image IDs may differ across builds even with identical package URLs, so retain both the explicit package list and image digest. Offline image export and registry publication require a separate storage/distribution decision.

Both recipes were built on the remote VM on 27 October 2026. The initial nine-job smoke workflow passed; the later validation record covers the expanded eleven-job workflow. The resolved MultiQC patch version is 1.27.1 (the initial YAML constraint was 1.27); explicit package locks, rather than the initial constraints alone, define the tested environment. See ../docs/validation/ for compact execution evidence. Biological execution uses two threads initially and a bounded container memory limit. No host daemon/socket is mounted into scientific containers. Container builds retain package/build caches; cleanup requires a deliberate review.
