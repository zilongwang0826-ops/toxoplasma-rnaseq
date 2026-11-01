"""Small, dependency-free provenance and input validation utilities."""
import csv
import hashlib
import json
import os
import re
from datetime import datetime, timezone
from pathlib import Path


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    os.replace(tmp, path)


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def read_samples(path, complete=True):
    with open(path, newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {"sample_id", "run_accession", "biosample", "condition", "day", "replicate", "layout",
                "fastq_1", "fastq_2", "md5_1", "md5_2", "bytes_1", "bytes_2"}
    if not rows or not required.issubset(rows[0]):
        raise ValueError("Sample sheet is empty or missing required columns: " + str(sorted(required)))
    for field in ("sample_id", "run_accession", "biosample"):
        if len({r[field] for r in rows}) != len(rows):
            raise ValueError("Expected one distinct biological library/run/BioSample per row: " + field)
    urls = [r[f"fastq_{m}"] for r in rows for m in (1, 2)]
    if len(urls) != len(set(urls)):
        raise ValueError("Duplicate FASTQ URL: potential pseudoreplication")
    counts = {}
    replicate_labels = set()
    for r in rows:
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]*", r["sample_id"]):
            raise ValueError("Unsafe sample identifier")
        if r["condition"] not in {"control", "infected"} or r["day"] not in {"11", "33"}:
            raise ValueError("Unexpected condition or time point")
        if r["layout"].upper() != "PAIRED":
            raise ValueError("This study workflow requires paired-end libraries")
        if not re.fullmatch(r"[SED]RR[0-9]+", r["run_accession"]):
            raise ValueError("Invalid run accession")
        for mate in (1, 2):
            url = r[f"fastq_{mate}"]
            expected = r["run_accession"] + f"_{mate}.fastq.gz"
            if url.rsplit("/", 1)[-1] != expected:
                raise ValueError("FASTQ URL does not match its run accession and mate")
            if not url.startswith("https://ftp.sra.ebi.ac.uk/"):
                raise ValueError("Expected verified ENA HTTPS download URL")
            if not re.fullmatch(r"[0-9a-f]{32}", r[f"md5_{mate}"]):
                raise ValueError("Missing/invalid source checksum")
            if int(r[f"bytes_{mate}"]) <= 0:
                raise ValueError("Invalid source byte count")
        label = (r["day"], r["condition"], r["replicate"])
        if r["replicate"] not in {"1", "2", "3"} or label in replicate_labels:
            raise ValueError("Invalid or repeated replicate label")
        replicate_labels.add(label)
        key = (r["day"], r["condition"])
        counts[key] = counts.get(key, 0) + 1
    if complete and (len(rows) != 12 or counts != {(d,c):3 for d in ("11","33") for c in ("control","infected")}):
        raise ValueError("Expected 12 libraries: three per day-by-condition cell")
    return rows
