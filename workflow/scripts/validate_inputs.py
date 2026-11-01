"""Validate the curated experimental design before any bulk download."""
import argparse
import csv
import json
import re
import shutil
from pathlib import Path
from common import read_samples, sha256, write_json, utc_now


def validate(samples, min_free_gib=0, max_download_gib=80, disk_path=".", source_dir=None):
    rows = read_samples(samples)
    source = Path(source_dir) if source_dir else Path(samples).parent / "source"
    archive = source / "ena_runs.tsv"
    manifests = json.loads((source / "sources.json").read_text())
    expected = next(x["sha256"] for x in manifests if x["file"] == "ena_runs.tsv")
    if sha256(archive) != expected:
        raise ValueError("Archived ENA metadata checksum mismatch")
    with archive.open(newline="") as handle:
        runs = {r["run_accession"]: r for r in csv.DictReader(handle, delimiter="\t")}
    for row in rows:
        archived = runs[row["run_accession"]]
        phase, status, replicate = re.fullmatch(r"(acutely|chronically) (uninfected|infected) group ([123])", archived["sample_alias"]).groups()
        identity = ("11" if phase == "acutely" else "33", "control" if status == "uninfected" else "infected", replicate)
        if (row["day"], row["condition"], row["replicate"]) != identity or row["biosample"] != archived["sample_accession"]:
            raise ValueError("Sample labels disagree with archived biological identity")
        for m in (1, 2):
            for source_column, manifest_column, prefix in [("fastq_ftp", "fastq", "https://"), ("fastq_md5", "md5", ""), ("fastq_bytes", "bytes", "")]:
                if row[f"{manifest_column}_{m}"] != prefix + archived[source_column].split(";")[m-1]:
                    raise ValueError("Sample file details disagree with archived source")
    total = sum(int(r[f"bytes_{m}"]) for r in rows for m in (1, 2))
    if total > max_download_gib * 1024**3:
        raise ValueError("Archive size exceeds configured download budget")
    free = shutil.disk_usage(disk_path).free
    if free < min_free_gib * 1024**3:
        raise ValueError("Insufficient free disk space for the configured execution budget")
    return {"validated_at": utc_now(), "samples_sha256": sha256(samples), "archive_metadata_sha256": expected, "n_libraries": len(rows),
            "archive_bytes": total, "disk_free_bytes": free,
            "experimental_unit": "library attributed to an individual mouse; pooling wording remains unresolved",
            "design": "independent infected versus control libraries within each day; n=3 per group",
            "status": "metadata_structure_validated; biological independence remains an assumption"}


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--samples", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--min-free-gib", type=float, default=0)
    p.add_argument("--max-download-gib", type=float, default=80)
    a = p.parse_args()
    write_json(a.output, validate(a.samples, a.min_free_gib, a.max_download_gib))
