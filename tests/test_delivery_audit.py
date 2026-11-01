"""Small synthetic delivery fixtures; no network, containers or biological results."""
import base64
import csv
import gzip
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import zipfile

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("delivery_audit", ROOT / "scripts/audit_completed_run.py")
audit = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(audit)
IMAGE = "sha256:" + "a" * 64
INDEX = {value: str(i + 1) * (128 if "512" in value else 64)
         for i, value in enumerate(audit.HASH_FIELDS.values())}


def put_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value))


def put_tsv(path, rows, fields=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields or list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def quant_text(order=("tx1.1", "tx2.1"), counts=(40, 40), tpms=(500000, 500000)):
    text = "Name\tLength\tEffectiveLength\tTPM\tNumReads\n"
    return text + "".join(f"{name}\t1000\t500\t{tpm}\t{count}\n" for name, count, tpm in zip(order, counts, tpms))


class DeliveryAudit(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="toxoplasma-audit-test-")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def sample(self):
        content = gzip.compress(b"@synthetic\nACGT\n+\nIIII\n", mtime=0)
        path = self.root / "raw.fastq.gz"
        path.write_bytes(content)
        row = {"sample_id": "S1", "read_pairs": "100", "bytes_1": str(len(content)),
               "md5_1": hashlib.md5(content).hexdigest(), "fastq_1": "https://example.org/synthetic.gz"}
        receipt = {"bytes": len(content), "url": row["fastq_1"], "md5": row["md5_1"],
                   "source_md5": row["md5_1"], "sha256": hashlib.sha256(content).hexdigest()}
        put_json(Path(str(path) + ".source.json"), receipt)
        return path, row

    def meta(self):
        return {"num_processed": 100, "num_mapped": 80, "percent_mapped": 80, "quant_errors": [], "library_types": ["IU"], "salmon_version": "1.10.3",
                **{key: INDEX[value] for key, value in audit.HASH_FIELDS.items()}}

    def test_raw_default_is_truthfully_receipt_only_and_rehash_detects_same_size_mutation(self):
        path, row = self.sample()
        data = bytearray(path.read_bytes())
        data[-1] ^= 1
        path.write_bytes(data)
        self.assertFalse(audit.audit_raw(path, row, 1)["rehash_performed"])
        with self.assertRaisesRegex(ValueError, "MD5 mismatch"):
            audit.audit_raw(path, row, 1, True)

    def test_wrong_receipt_source_is_rejected(self):
        path, row = self.sample()
        receipt = audit.read_json(str(path) + ".source.json")
        receipt["url"] = "https://example.org/wrong.gz"
        put_json(Path(str(path) + ".source.json"), receipt)
        with self.assertRaisesRegex(ValueError, "URL mismatch"):
            audit.audit_raw(path, row, 1)

    def test_quant_order_counts_reference_and_nonfinite_checks(self):
        path = self.root / "quant.sf"
        path.write_text(quant_text())
        record, names = audit.audit_quant(path, self.meta(), {"sample_id": "S1", "read_pairs": "100"}, INDEX, {"tx1.1", "tx2.1"})
        self.assertEqual(record["assigned_fragments"], 80)
        self.assertEqual(record["md5"], hashlib.md5(path.read_bytes()).hexdigest())
        for text, message in [(quant_text(order=("tx2.1", "tx1.1")), "order differs"),
                              (quant_text(counts=(50, 50)), "mapped pairs"),
                              (quant_text(counts=("nan", 40)), "Invalid quantification"),
                              (quant_text(order=("tx1.1", "unknown")), "target IDs"),
                              (quant_text(tpms=(1, 2)), "TPM total")]:
            with self.subTest(message=message):
                path.write_text(text)
                with self.assertRaisesRegex(ValueError, message):
                    audit.audit_quant(path, self.meta(), {"sample_id": "S1", "read_pairs": "100"}, INDEX, set(names), names)
        path.write_text(quant_text())
        meta = self.meta()
        meta["index_decoy_seq_hash"] = "wrong"
        with self.assertRaisesRegex(ValueError, "Index signature"):
            audit.audit_quant(path, meta, {"sample_id": "S1", "read_pairs": "100"}, INDEX, set(names))
        for field, value, message in [("quant_errors", ["failed"], "quantification errors"),
                                      ("percent_mapped", 0.8, "mapping percentage"),
                                      ("salmon_version", "1.2.0", "Salmon version")]:
            meta = self.meta()
            meta[field] = value
            with self.subTest(field=field), self.assertRaisesRegex(ValueError, message):
                audit.audit_quant(path, meta, {"sample_id": "S1", "read_pairs": "100"}, INDEX, set(names))

    def test_fastqc_and_local_report_targets(self):
        path = self.root / "S1_1_fastqc.zip"
        with zipfile.ZipFile(path, "w") as archive:
            archive.writestr("S1_1_fastqc/summary.txt", "".join(f"PASS\t{x}\tS1_1.fastq.gz\n" for x in audit.CRITICAL_QC))
            archive.writestr("S1_1_fastqc/fastqc_data.txt", "Total Sequences\t100\n")
        self.assertEqual(audit.audit_fastqc(path, "S1_1", 100)["total_sequences"], 100)
        with self.assertRaisesRegex(ValueError, "count differs"):
            audit.audit_fastqc(path, "S1_1", 101)
        report = self.root / "index.html"
        report.write_text('<img src="figure%20one.png"><a href="https://example.org">Source</a>')
        with self.assertRaisesRegex(ValueError, "Missing/empty report target"):
            audit.audit_links(report, self.root)
        (self.root / "figure one.png").write_bytes(b"synthetic image placeholder for target existence only")
        self.assertEqual(len(audit.audit_links(report, self.root)["local_targets"]), 1)
        (self.root / "tables").mkdir()
        report.write_text('<a href="tables/">Tables</a>')
        self.assertEqual(audit.audit_links(report, self.root)["local_targets"][0]["kind"], "directory")
        report.write_text('<a href="../outside.tsv">Outside</a>')
        with self.assertRaisesRegex(ValueError, "escapes"):
            audit.audit_links(report, self.root)

    def metadata(self, output, start, end, incomplete=False):
        filename = base64.b64encode(output.encode()).decode()
        path = self.root / ".snakemake/metadata" / filename
        put_json(path, {"rule": "synthetic", "starttime": start, "endtime": end, "incomplete": incomplete,
                        "shellcmd": f"python wrapper --expected-id {IMAGE} -- synthetic", "input_checksums": {}})
        return filename

    def test_stage_revision_is_temporal_not_final_head_and_markers_are_classified(self):
        launches = [{"timestamp": 100, "path": "first.txt", "revision": "first"},
                    {"timestamp": 200, "path": "second.txt", "revision": "second"}]
        filename = self.metadata("results/current.tsv", 150, 160)
        (self.root / "results").mkdir()
        (self.root / "results/current.tsv").write_text("synthetic")
        marker = self.root / ".snakemake/incomplete" / filename
        put_json(marker, {"external_jobid": None})
        os.utime(marker, (140, 140))
        other = self.root / ".snakemake/incomplete" / base64.b64encode(b"results/old-attempt.tsv").decode()
        put_json(other, {})
        stages, markers, errors = audit.stage_evidence(self.root, {"results/current.tsv"}, launches)
        self.assertFalse(errors)
        self.assertEqual(stages[0]["associated_launch_revision"], "first")
        self.assertIn("not per-job Git attestation", stages[0]["revision_evidence"])
        self.assertEqual({x["classification"] for x in markers}, {"historical_or_unrelated_output", "superseded_by_later_complete_record"})
        os.utime(marker, (170, 170))
        self.assertIn("Unresolved incomplete", audit.stage_evidence(self.root, {"results/current.tsv"}, launches)[2][0])

    def fixture(self):
        cfg = {"samples": "metadata/samples.tsv", "raw": "data/fastq", "output": "results", "reference_dir": "resources/test",
               "reference": {"decoys": True}, "quantification": {"trim": False}, "analysis": {"interaction": False},
               "enrichment": {"kegg_tsv": ""}}
        samples = []
        payload = gzip.compress(b"@synthetic\nACGT\n+\nIIII\n", mtime=0)
        for day in (11, 33):
            for condition in ("control", "infected"):
                for replicate in (1, 2, 3):
                    n = len(samples) + 1
                    row = {"sample_id": f"S{n}", "run_accession": f"SRR{n}", "biosample": f"SAMN{n}", "condition": condition,
                           "day": str(day), "replicate": str(replicate), "layout": "PAIRED", "read_pairs": "100"}
                    for mate in (1, 2):
                        row.update({f"fastq_{mate}": f"https://ftp.sra.ebi.ac.uk/synthetic/SRR{n}_{mate}.fastq.gz",
                                    f"md5_{mate}": hashlib.md5(payload).hexdigest(), f"bytes_{mate}": str(len(payload))})
                    samples.append(row)
        required, requested = audit.declared_artifacts(cfg, samples)
        directory_outputs = {"results/analysis", "results/enrichment", "results/comparison", "resources/test/salmon_index",
                             "results/sensitivity/lilrb4a", "results/figures", "results/validation/models"}
        for name in sorted(required | requested):
            path = self.root / name
            if name in directory_outputs:
                path.mkdir(parents=True, exist_ok=True)
            else:
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text("Synthetic fixture: content existence only\n")
        put_tsv(self.root / cfg["samples"], samples)
        put_tsv(self.root / "resources/test/prepared/tx2gene.tsv", [{"transcript_id": "tx1.1", "gene_id": "gene1"}, {"transcript_id": "tx2.1", "gene_id": "gene2"}])
        put_json(self.root / "resources/test/prepared/reference.json", {"transcripts": 2, "decoy_sequences": 1})
        put_json(self.root / "resources/test/salmon_index/info.json", INDEX)
        put_json(self.root / "results/provenance/run.json", {"effective_config": cfg, "commit": "final-do-not-impute-to-stages", "tracked_source_dirty": False,
                                                            "containers": {"images": {"tools": {"id": IMAGE}}}})
        for row in samples:
            sid = row["sample_id"]
            for mate in (1, 2):
                path = self.root / f"data/fastq/{sid}_{mate}.fastq.gz"
                path.write_bytes(payload)
                put_json(Path(str(path) + ".source.json"), {"url": row[f"fastq_{mate}"], "bytes": len(payload), "md5": row[f"md5_{mate}"],
                                                          "source_md5": row[f"md5_{mate}"], "sha256": hashlib.sha256(payload).hexdigest()})
                with zipfile.ZipFile(self.root / f"results/qc/raw/{sid}_{mate}_fastqc.zip", "w") as archive:
                    name = f"{sid}_{mate}_fastqc"
                    archive.writestr(name + "/summary.txt", "".join(f"PASS\t{x}\tsynthetic\n" for x in audit.CRITICAL_QC))
                    archive.writestr(name + "/fastqc_data.txt", "Total Sequences\t100\n")
            (self.root / f"results/salmon/{sid}/quant.sf").write_text(quant_text())
            put_json(self.root / f"results/salmon/{sid}/aux_info/meta_info.json", self.meta())
        (self.root / "results/launch/study.exit").write_text("0\n")
        put_json(self.root / "results/qc/quantification.json", {"passed": True})
        put_tsv(self.root / "results/analysis/summary.tsv", [{"level": level, "day": day, "n_control": 3, "n_infected": 3} for level in ("gene", "transcript") for day in (11, 33)])
        put_tsv(self.root / "results/enrichment/summary.tsv", [{"day": day, "family": "GO", "direction": direction} for day in (11, 33) for direction in ("up", "down")])
        put_tsv(self.root / "results/comparison/summary.tsv", [{"day": day, "selection": selection} for day in (11, 33) for selection in ("significant", "nominal_sensitivity")])
        (self.root / "results/report/index.html").write_text('<img src="../analysis/gene/qc/pca.png"><a href="../provenance/run.json">Provenance</a>')
        for name in ("analysis", "enrichment", "comparison", "sensitivity/lilrb4a", "figures", "validation/models"):
            input_path = "metadata/samples.tsv"
            put_tsv(self.root / f"results/{name}/input-checksums.tsv", [{"path": input_path, "md5": hashlib.md5((self.root / input_path).read_bytes()).hexdigest()}])
        (self.root / "results/launch/20200913T122640Z.txt").write_text("a" * 40 + "\nCONFIG=config/config.yaml\nPROFILE=profiles/local\n./run.sh run\n")
        for output in requested:
            self.metadata(output, 1600000001, 1600000002)
        return cfg

    def test_complete_synthetic_audit_is_published_once_and_detects_missing_deliverable(self):
        self.fixture()
        result = audit.audit_project(self.root, self.root / "audit-one")
        self.assertEqual(result["status"], "PASS", result["errors"])
        self.assertEqual(len(result["raw_files"]), 24)
        self.assertEqual(len(result["quantifications"]), 12)
        self.assertIn("no_raw_rehash", result["raw_verification_level"])
        snapshot = (self.root / "audit-one/audit.json").read_bytes()
        with self.assertRaisesRegex(ValueError, "already exists"):
            audit.audit_project(self.root, self.root / "audit-one")
        self.assertEqual(snapshot, (self.root / "audit-one/audit.json").read_bytes())
        (self.root / "results/analysis/transcript/day33/results.tsv").unlink()  # This test's synthetic fixture only.
        failed = audit.audit_project(self.root, self.root / "audit-two", rehash_raw=True)
        self.assertEqual(failed["status"], "FAIL")
        self.assertTrue(any("transcript/day33/results.tsv" in e for e in failed["errors"]))
        self.assertEqual(failed["raw_verification_level"], "streaming_md5_sha256_and_receipts")

    def test_failing_audit_is_still_recorded(self):
        result = audit.audit_project(self.root, self.root / "failed-audit")
        self.assertEqual(result["status"], "FAIL")
        self.assertTrue((self.root / "failed-audit/audit.json").is_file())
        self.assertTrue((self.root / "failed-audit/stage-execution.tsv").is_file())


if __name__ == "__main__":
    unittest.main()
