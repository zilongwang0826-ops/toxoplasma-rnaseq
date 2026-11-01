#!/usr/bin/env python3
"""Audit a completed study without modifying its inputs or execution state.

Usage: python3 scripts/audit_completed_run.py --output results/audits/analysis
The output directory must not exist. Add --rehash-raw for streaming MD5/SHA-256
verification of all compressed reads; otherwise raw verification is receipt/size
consistency only. No network, containers, workflow jobs, or Git writes are used.
"""
import argparse
import base64
import csv
from datetime import datetime, timezone
import hashlib
from html.parser import HTMLParser
import json
import math
from pathlib import Path
import re
import shlex
import subprocess
import sys
from urllib.parse import unquote, urlsplit
import zipfile

sys.dont_write_bytecode = True  # The audit must not create project __pycache__ files.
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "workflow/scripts"))
from common import read_samples

CRITICAL_QC = ("Per base sequence quality", "Per base N content", "Adapter Content")
HASH_FIELDS = {"index_seq_hash": "SeqHash", "index_name_hash": "NameHash",
               "index_seq_hash512": "SeqHash512", "index_name_hash512": "NameHash512",
               "index_decoy_seq_hash": "DecoySeqHash", "index_decoy_name_hash": "DecoyNameHash"}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read_json(path):
    return json.loads(Path(path).read_text())


def read_tsv(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def digests(path):
    md5, sha = hashlib.md5(), hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            md5.update(block)
            sha.update(block)
    return {"md5": md5.hexdigest(), "sha256": sha.hexdigest()}


def project_path(root, value):
    path = (root / value).resolve()
    require(path.is_relative_to(root), "Path escapes project: " + str(value))
    return path


def audit_raw(path, sample, mate, rehash=False):
    receipt = read_json(str(path) + ".source.json")
    expected_bytes = int(sample[f"bytes_{mate}"])
    expected_md5 = sample[f"md5_{mate}"]
    require(path.stat().st_size == expected_bytes, "Raw file byte count mismatch")
    require(receipt.get("bytes") == expected_bytes, "Receipt byte count mismatch")
    require(receipt.get("url") == sample[f"fastq_{mate}"], "Receipt URL mismatch")
    require(receipt.get("md5") == expected_md5, "Receipt MD5 differs from manifest")
    require(receipt.get("source_md5") == expected_md5, "Receipt source MD5 differs from manifest")
    require(bool(re.fullmatch(r"[0-9a-f]{64}", receipt.get("sha256", ""))), "Invalid receipt SHA-256")
    result = {"sample_id": sample["sample_id"], "mate": mate, "bytes": expected_bytes,
              "url": receipt["url"], "receipt_md5": expected_md5,
              "receipt_sha256": receipt["sha256"], "rehash_performed": rehash}
    if rehash:
        actual = digests(path)
        require(actual["md5"] == expected_md5, "Raw content MD5 mismatch")
        require(actual["sha256"] == receipt["sha256"], "Raw content SHA-256 mismatch")
        result["computed"] = actual
    return result


def audit_quant(path, meta, sample, index, expected_ids, previous_ids=None):
    """Read each quantification once, including content hashes and numeric checks."""
    require(not meta.get("quant_errors"), "Salmon reports quantification errors")
    require(meta.get("salmon_version") == "1.10.3", "Salmon version differs from the pinned 1.10.3 environment")
    md5, sha = hashlib.md5(), hashlib.sha256()
    names, seen = [], set()
    counts, tpm = 0.0, 0.0
    with Path(path).open("rb") as handle:
        def lines():
            for line in handle:
                md5.update(line)
                sha.update(line)
                yield line.decode("utf-8")
        reader = csv.DictReader(lines(), delimiter="\t")
        require(reader.fieldnames == ["Name", "Length", "EffectiveLength", "TPM", "NumReads"],
                "Unexpected quant.sf columns")
        for row in reader:
            name = row["Name"]
            require(name and name not in seen and None not in row, "Duplicate/empty ID or malformed quant row")
            values = [float(row[k]) for k in ("Length", "EffectiveLength", "TPM", "NumReads")]
            require(all(math.isfinite(x) and x >= 0 for x in values), "Invalid quantification number")
            require(values[0] > 0, "Nonpositive transcript length")
            names.append(name)
            seen.add(name)
            tpm += values[2]
            counts += values[3]
    require(seen == expected_ids, "Quantification target IDs differ from prepared reference")
    require(previous_ids is None or names == previous_ids, "Quantification target order differs across samples")
    require(int(meta["num_processed"]) == int(sample["read_pairs"]), "Incomplete archive pair processing")
    mapped = int(meta["num_mapped"])
    require(0 < mapped <= int(meta["num_processed"]), "Invalid mapped-pair count")
    percent = float(meta["percent_mapped"])
    require(math.isfinite(percent) and abs(percent - 100 * mapped / int(meta["num_processed"])) <= 0.001,
            "Reported mapping percentage disagrees with processed/mapped counts")
    # Salmon 1.10 writes estimated counts to three decimal places. The sum of
    # individual rounding errors is bounded by half a unit per target.
    count_tolerance = max(1.0, len(names) * 0.0005 + 1e-6)
    require(abs(counts - mapped) <= count_tolerance, "Assigned counts disagree with mapped pairs")
    require(abs(tpm - 1e6) <= 100, "TPM total differs from 1e6 by more than 100")
    for key, index_key in HASH_FIELDS.items():
        require(meta.get(key) == index.get(index_key) and bool(index.get(index_key)),
                "Index signature mismatch or missing: " + key)
    types = meta.get("library_types", [])
    require(len(types) == 1 and types[0] in {"IU", "ISF", "ISR", "OU", "OSF", "OSR", "MU", "MSF", "MSR"},
            "Unresolved paired library type")
    result = {"sample_id": sample["sample_id"], "targets": len(names), "processed_pairs": int(meta["num_processed"]),
              "mapped_pairs": mapped, "assigned_fragments": counts, "count_rounding_tolerance": count_tolerance,
              "tpm_sum": tpm, "library_type": types[0], "salmon_version": meta.get("salmon_version"),
              "md5": md5.hexdigest(), "sha256": sha.hexdigest()}
    return result, names


def audit_fastqc(path, name, pairs):
    with zipfile.ZipFile(path) as archive:
        states = {row[1]: row[0] for row in csv.reader(
            archive.read(name + "_fastqc/summary.txt").decode().splitlines(), delimiter="\t")}
        data = archive.read(name + "_fastqc/fastqc_data.txt").decode().splitlines()
    counts = [int(line.split("\t")[1]) for line in data if line.startswith("Total Sequences\t")]
    require(counts == [pairs], "FastQC read count differs from archive")
    require(all(states.get(module) == "PASS" for module in CRITICAL_QC), "Critical FastQC module is not PASS")
    return {"sample_mate": name, "total_sequences": pairs, "modules": states}


class LocalLinks(HTMLParser):
    def __init__(self):
        super().__init__()
        self.links = []

    def handle_starttag(self, tag, attrs):
        attributes = dict(attrs)
        for key in ("href", "src"):
            if attributes.get(key):
                self.links.append((tag, attributes[key]))


def audit_links(report, root):
    root = Path(root).resolve()
    report = Path(report).resolve()
    parser = LocalLinks()
    parser.feed(report.read_text())
    checked, external = [], []
    for tag, value in parser.links:
        url = urlsplit(value)
        require(url.scheme != "file", "Nonportable file-URL report target: " + value)
        if url.scheme or url.netloc:
            external.append(value)
            continue
        if not url.path:
            continue
        target = (report.parent / unquote(url.path)).resolve()
        require(target.is_relative_to(root), "Report link escapes delivered project: " + value)
        require(target.is_dir() or (target.is_file() and target.stat().st_size > 0), "Missing/empty report target: " + value)
        checked.append({"tag": tag, "url": value, "target": str(target.relative_to(root)),
                        "kind": "directory" if target.is_dir() else "file"})
    return {"local_targets": checked, "external_urls_not_requested": sorted(set(external)),
            "scope": "Local target existence only; no browser rendering or external URL checks"}


def launch_records(root, results):
    records = []
    for path in sorted((results / "launch").glob("*.txt")):
        if not re.fullmatch(r"\d{8}T\d{6}Z", path.stem):
            continue
        lines = path.read_text().splitlines()
        if not lines or not re.fullmatch(r"[0-9a-f]{40,64}", lines[0]):
            continue
        timestamp = datetime.strptime(path.stem, "%Y%m%dT%H%M%SZ").replace(tzinfo=timezone.utc).timestamp()
        records.append({"path": str(path.relative_to(root)), "timestamp": timestamp, "revision": lines[0],
                        "recorded_command": lines[-1] if len(lines) > 1 else "", "sha256": digests(path)["sha256"]})
    return records


def decode_output(filename):
    try:
        return base64.urlsafe_b64decode(filename + "=" * (-len(filename) % 4)).decode("utf-8")
    except (ValueError, UnicodeError):
        return None


def stage_evidence(root, requested, launches):
    """Associate actual job timestamps with launches, explicitly as an inference."""
    records, errors, groups = {}, [], {}
    directory = root / ".snakemake/metadata"
    for path in sorted(directory.rglob("*")):
        if not path.is_file():
            continue
        output = decode_output(str(path.relative_to(directory)))
        if output not in requested:
            continue
        data = read_json(path)
        records[output] = data
        if data.get("incomplete") or not isinstance(data.get("endtime"), (int, float)):
            errors.append("Incomplete execution metadata: " + output)
        key = (data.get("rule"), data.get("starttime"), data.get("endtime"), data.get("shellcmd"))
        group = groups.setdefault(key, {"outputs": [], "metadata_records": [], "data": data})
        group["outputs"].append(output)
        group["metadata_records"].append(str(path.relative_to(root)))
    for output in sorted(requested - records.keys()):
        errors.append("Missing execution metadata: " + output)
    stages = []
    for group in groups.values():
        data = group["data"]
        start = data.get("starttime")
        candidates = [x for x in launches if isinstance(start, (float, int)) and x["timestamp"] <= start]
        launch = max(candidates, key=lambda x: x["timestamp"]) if candidates else None
        command = data.get("shellcmd", "")
        parts = shlex.split(command)
        image_id = parts[parts.index("--expected-id") + 1] if "--expected-id" in parts else ""
        stages.append({"rule": data.get("rule", ""), "outputs": json.dumps(sorted(group["outputs"])),
                       "start_epoch": start, "end_epoch": data.get("endtime"), "incomplete": data.get("incomplete", False),
                       "image_id": image_id, "shell_command": command,
                       "associated_launch": launch["path"] if launch else "",
                       "associated_launch_revision": launch["revision"] if launch else "unknown",
                       "revision_evidence": "latest recorded launch preceding job start; temporal association, not per-job Git attestation" if launch else "no preceding launch record",
                       "metadata_records": json.dumps(sorted(group["metadata_records"])),
                       "recorded_input_checksums": json.dumps(data.get("input_checksums", {}), sort_keys=True)})
    markers = []
    incomplete = root / ".snakemake/incomplete"
    for path in sorted(incomplete.rglob("*")):
        if not path.is_file():
            continue
        output = decode_output(str(path.relative_to(incomplete)))
        metadata = records.get(output, {})
        if output not in requested:
            classification = "historical_or_unrelated_output"
        elif (not metadata.get("incomplete", True) and isinstance(metadata.get("endtime"), (float, int))
              and metadata["endtime"] >= path.stat().st_mtime and (root / output).exists()):
            classification = "superseded_by_later_complete_record"
        else:
            classification = "unresolved_requested_output"
            errors.append("Unresolved incomplete marker for requested output: " + str(output))
        markers.append({"path": str(path.relative_to(root)), "output": output, "classification": classification})
    return sorted(stages, key=lambda x: (x["start_epoch"] or 0, x["rule"])), markers, errors


def declared_artifacts(cfg, samples):
    """Files and job-output identities of this study's declared workflow."""
    out, raw, ref = cfg["output"], cfg["raw"], cfg["reference_dir"]
    required = {f"{out}/{name}" for name in ["analysis/summary.tsv", "analysis/execution-status.txt", "analysis/parameters.tsv",
        "analysis/session-info.txt", "analysis/input-checksums.tsv", "analysis/samples-used.tsv", "enrichment/summary.tsv",
        "enrichment/execution-status.txt", "enrichment/parameters.tsv", "enrichment/session-info.txt", "enrichment/input-checksums.tsv",
        "enrichment/go-annotation-snapshot.tsv", "enrichment/go-database-metadata.tsv", "comparison/summary.tsv",
        "comparison/input-checksums.tsv", "comparison/session-info.txt", "comparison/interpretation.txt", "provenance/run.json",
        "report/index.html", "qc/multiqc_report.html", "qc/quantification.json", "qc/quantification.json.tsv", "launch/study.exit"]}
    requested = {f"{out}/{name}" for name in ["validation/metadata.json", "analysis", "enrichment", "comparison",
        "provenance/run.json", "report/index.html", "qc/multiqc_report.html", "qc/quantification.json", "qc/quantification.json.tsv",
        "sensitivity/lilrb4a", "figures", "validation/models"]}
    required.update(f"{out}/{name}" for name in ["sensitivity/lilrb4a/summary.tsv", "sensitivity/lilrb4a/interpretation.txt",
        "sensitivity/lilrb4a/input-checksums.tsv", "sensitivity/lilrb4a/pair-de-status.tsv", "sensitivity/lilrb4a/locus-definition.tsv",
        "figures/interpretation.txt", "figures/input-checksums.tsv", "validation/models/models.tsv",
        "validation/models/input-checksums.tsv", "validation/models/diagnostics.json"])
    requested.add(f"{ref}/salmon_index")
    requested.update(f"{ref}/prepared/{name}" for name in ("gentrome.fa", "tx2gene.tsv", "decoys.txt", "reference.json"))
    for kind in ["cdna", "ncrna", "gtf"] + (["genome"] if cfg["reference"]["decoys"] else []):
        requested.update((f"{ref}/source/{kind}.gz", f"{ref}/source/{kind}.gz.source.json"))
    for sample in samples:
        sid = sample["sample_id"]
        for mate in (1, 2):
            requested.update((f"{raw}/{sid}_{mate}.fastq.gz", f"{raw}/{sid}_{mate}.fastq.gz.source.json"))
            requested.update(f"{out}/qc/raw/{sid}_{mate}_fastqc.{ext}" for ext in ("zip", "html"))
        requested.update(f"{out}/salmon/{sid}/{name}" for name in ("quant.sf", "aux_info/meta_info.json", "lib_format_counts.json"))
    for level in ("gene", "transcript"):
        for scope in ("qc", "day11/figures", "day33/figures"):
            for stem in ("pca", "sample-distances"):
                required.update(f"{out}/figures/{level}/{scope}/{stem}.{ext}" for ext in ("pdf", "png"))
        for stem in ("pca", "sample-distances", "variable-features-heatmap"):
            required.update(f"{out}/analysis/{level}/qc/{stem}.{ext}" for ext in ("pdf", "png"))
        if cfg["analysis"].get("interaction"):
            required.update(f"{out}/analysis/{level}/interaction/{name}" for name in ("results.tsv", "interpretation.txt", "deseq2.rds"))
        for day in (11, 33):
            prefix = f"{out}/analysis/{level}/day{day}"
            required.update(f"{prefix}/{name}" for name in ("results.tsv", "significant.tsv", "normalized-counts.tsv", "feature-filter.tsv", "deseq2.rds", "summary.tsv", "result-metadata.txt"))
            for stem in ("pca", "sample-distances", "variable-features-heatmap", "ma", "volcano"):
                required.update(f"{prefix}/figures/{stem}.{ext}" for ext in ("pdf", "png"))
    families = ["GO"] + (["KEGG"] if cfg.get("enrichment", {}).get("kegg_tsv") else [])
    if "KEGG" not in families:
        required.add(f"{out}/enrichment/kegg-status.txt")
    for day in (11, 33):
        required.add(f"{out}/comparison/day{day}-published-transcript-comparison.tsv")
        required.update(f"{out}/comparison/day{day}-effect-comparison.{ext}" for ext in ("pdf", "png"))
        for family in families:
            prefix = f"{out}/enrichment/day{day}/{family}"
            required.add(f"{prefix}/universe.tsv")
            for direction in ("up", "down"):
                required.add(f"{prefix}/{direction}.tsv")
                required.update(f"{prefix}/{direction}-enrichment.{ext}" for ext in ("pdf", "png"))
    return required, requested


def audit_project(project, output, rehash_raw=False, results="results", additional_outputs=()):
    root = Path(project).resolve()
    destination = Path(output).resolve()
    require(not destination.exists(), "Audit output already exists; choose a new directory")
    destination.mkdir(parents=True, exist_ok=False)
    errors, warnings, fingerprints = [], [], {}
    evidence = {"schema_version": 1, "started_utc": datetime.now(timezone.utc).isoformat(),
                "raw_verification_level": "streaming_md5_sha256_and_receipts" if rehash_raw else "receipt_manifest_and_file_size_only_no_raw_rehash",
                "raw_files": [], "quantifications": [], "fastqc": []}
    stages = []

    def attempt(label, function):
        try:
            return function()
        except (OSError, ValueError, TypeError, KeyError, IndexError, zipfile.BadZipFile) as exc:
            errors.append(label + ": " + str(exc))
            return None

    def fingerprint(path):
        key = str(path.relative_to(root))
        if key not in fingerprints:
            fingerprints[key] = {"bytes": path.stat().st_size, **digests(path)}
        return fingerprints[key]

    try:
        result_root = project_path(root, results)
        provenance = read_json(result_root / "provenance/run.json")
        cfg = provenance["effective_config"]
        require(project_path(root, cfg["output"]) == result_root, "Provenance output root differs from selected results")
        require(not cfg["quantification"]["trim"], "This complete-archive audit supports the executed untrimmed study only")
        sample_path = project_path(root, cfg["samples"])
        samples = read_samples(sample_path)
        require(all(int(row["read_pairs"]) > 0 for row in samples), "Missing/invalid archive read-pair counts")
        fingerprint(sample_path)
        reference = project_path(root, cfg["reference_dir"])
        mapping = read_tsv(reference / "prepared/tx2gene.tsv")
        expected_ids = {row["transcript_id"] for row in mapping}
        require(len(expected_ids) == len(mapping) and all(expected_ids), "Invalid prepared target identity")
        index = read_json(reference / "salmon_index/info.json")
        fingerprint(reference / "salmon_index/info.json")
        fingerprint(reference / "prepared/reference.json")
        previous_ids = None
        for sample in samples:
            sid = sample["sample_id"]
            for mate in (1, 2):
                raw = project_path(root, f"{cfg['raw']}/{sid}_{mate}.fastq.gz")
                record = attempt(f"Raw {sid}/{mate}", lambda: audit_raw(raw, sample, mate, rehash_raw))
                if record is not None:
                    evidence["raw_files"].append(record)
                qcpath = result_root / "qc/raw" / f"{sid}_{mate}_fastqc.zip"
                record = attempt(f"FastQC {sid}/{mate}", lambda: audit_fastqc(qcpath, f"{sid}_{mate}", int(sample["read_pairs"])))
                if record is not None:
                    evidence["fastqc"].append(record)
            folder = result_root / "salmon" / sid
            record = attempt("Quantification " + sid, lambda: audit_quant(folder / "quant.sf", read_json(folder / "aux_info/meta_info.json"), sample, index, expected_ids, previous_ids))
            if record is not None:
                quant, previous_ids = record
                fingerprints[str((folder / "quant.sf").relative_to(root))] = {k: quant[k] for k in ("md5", "sha256")}
                evidence["quantifications"].append(quant)
        require(len(evidence["raw_files"]) == 24 and len(evidence["quantifications"]) == 12 and len(evidence["fastqc"]) == 24,
                "Incomplete per-library audit coverage; inspect preceding errors")
        require(len({x["library_type"] for x in evidence["quantifications"]}) == 1, "Inconsistent library types")
        required, requested = declared_artifacts(cfg, samples)
        for name in additional_outputs:
            require(project_path(root, name).is_relative_to(result_root), "Additional output must be inside results")
            requested.add(name)
        evidence["artifacts"] = []
        for name in sorted(required | requested):
            path = project_path(root, name)
            if not path.exists() or (path.is_file() and path.stat().st_size == 0):
                errors.append("Missing/empty declared artifact: " + name)
            else:
                evidence["artifacts"].append({"path": name, "bytes": path.stat().st_size if path.is_file() else None})
        require((result_root / "launch/study.exit").read_text().strip() == "0", "Study driver did not record exit 0")
        require(read_json(result_root / "qc/quantification.json").get("passed") is True, "Workflow quantification gate did not pass")
        analyses = read_tsv(result_root / "analysis/summary.tsv")
        require(len(analyses) == 4 and {(x["level"], x["day"]) for x in analyses} == {(level, day) for level in ("gene", "transcript") for day in ("11", "33")}, "Analysis summary is not the four intended contrasts")
        require(all(int(x["n_control"]) == int(x["n_infected"]) == 3 for x in analyses), "Analysis replicate counts differ from design")
        enrichment = read_tsv(result_root / "enrichment/summary.tsv")
        families = {"GO", "KEGG"} if cfg.get("enrichment", {}).get("kegg_tsv") else {"GO"}
        require(len(enrichment) == 4 * len(families) and {(x["day"], x["family"], x["direction"]) for x in enrichment} == {(d, f, s) for d in ("11", "33") for f in families for s in ("up", "down")}, "Enrichment summary families are incomplete")
        comparison = read_tsv(result_root / "comparison/summary.tsv")
        require(len(comparison) == 4 and {(x["day"], x["selection"]) for x in comparison} == {(d, s) for d in ("11", "33") for s in ("significant", "nominal_sensitivity")}, "Comparison summary is incomplete")
        evidence["report_links"] = attempt("Report links", lambda: audit_links(result_root / "report/index.html", result_root))
        for directory in ("analysis", "enrichment", "comparison", "sensitivity/lilrb4a", "figures", "validation/models"):
            for row in read_tsv(result_root / directory / "input-checksums.tsv"):
                actual = attempt("Input checksum " + row["path"], lambda: fingerprint(project_path(root, row["path"])))
                if actual and actual["md5"] != row["md5"]:
                    errors.append("Executed analysis input changed: " + row["path"])
        for name in ("analysis/summary.tsv", "enrichment/summary.tsv", "comparison/summary.tsv", "provenance/run.json", "report/index.html", "qc/quantification.json"):
            fingerprint(result_root / name)
        evidence["launch_records"] = launch_records(root, result_root)
        stages, evidence["incomplete_markers"], stage_errors = stage_evidence(root, requested, evidence["launch_records"])
        errors.extend(stage_errors)
        if any(x["classification"] != "unresolved_requested_output" for x in evidence["incomplete_markers"]):
            warnings.append("Historical/superseded incomplete markers retained; none were deleted. A subsequent dry run should assess scheduler resume behavior.")
        require(bool(stages), "No stage execution metadata available")
        exit_path = result_root / "launch/study.exit"
        evidence["study_driver_exit"] = {"code": 0, "mtime_epoch": exit_path.stat().st_mtime,
                                        "scope": "Original driver completion; supplementary stages may run later"}
        if exit_path.stat().st_mtime < max(x["end_epoch"] or 0 for x in stages):
            warnings.append("Completed supplementary stages postdate the original study.exit marker; their completion is attested by the individual stage records, not that earlier exit file.")
        image_ids = {v["id"] for v in provenance["containers"]["images"].values()}
        for stage in stages:
            if stage["image_id"] not in image_ids:
                errors.append("Missing/unrecognized executed image: " + stage["rule"])
        evidence["provenance_recorded_end_commit"] = provenance.get("commit")
        evidence["revision_interpretation"] = "End commit is a provenance snapshot. Stage revisions are associated with actual launch records by time, not assumed equal to the end commit."
        if provenance.get("tracked_source_dirty"):
            errors.append("Workflow provenance reports dirty tracked source")
        try:
            evidence["audit_host_head"] = subprocess.check_output(["git", "-C", str(root), "rev-parse", "HEAD"], text=True, stderr=subprocess.DEVNULL).strip()
        except (OSError, subprocess.CalledProcessError):
            warnings.append("Git HEAD unavailable in this delivered copy; recorded execution revision retained")
    except (OSError, ValueError, TypeError, KeyError, IndexError) as exc:
        errors.append("Audit prerequisite or summary check: " + str(exc))
    evidence.update(status="PASS" if not errors else "FAIL", errors=errors, warnings=warnings,
                    fingerprints=fingerprints, completed_utc=datetime.now(timezone.utc).isoformat(),
                    limitations=["No scientific interpretation or independent reproduction claim", "No workflow rerun, browser rendering, container execution or external requests",
                                 "Raw content not rehashed" if not rehash_raw else "Raw hashes verified at audit time; subsequent mutations are not monitored"])
    columns = ["rule", "outputs", "start_epoch", "end_epoch", "incomplete", "image_id", "shell_command", "associated_launch", "associated_launch_revision", "revision_evidence", "metadata_records", "recorded_input_checksums"]
    with (destination / "stage-execution.tsv").open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(stages)
    evidence["stage_execution_sha256"] = digests(destination / "stage-execution.tsv")["sha256"]
    evidence["audit_script_sha256"] = digests(Path(__file__))["sha256"]
    with (destination / "audit.json").open("x") as handle:
        json.dump(evidence, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")
    for name in ("audit.json", "stage-execution.tsv"):
        (destination / name).chmod(0o444)
    return evidence


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", default=".")
    parser.add_argument("--results", default="results", help="Project-relative results directory containing provenance/run.json")
    parser.add_argument("--output", required=True, help="New audit directory; existing paths are refused")
    parser.add_argument("--rehash-raw", action="store_true")
    parser.add_argument("--additional-output", action="append", default=[],
                        help="Additional current Snakemake output, project-relative and inside results; repeat as needed")
    args = parser.parse_args()
    try:
        result = audit_project(args.project, args.output, args.rehash_raw, args.results, args.additional_output)
    except (OSError, ValueError) as exc:
        parser.exit(2, str(exc) + "\n")
    print(json.dumps({"status": result["status"], "errors": result["errors"], "output": str(Path(args.output).resolve())}))
    return 0 if result["status"] == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
