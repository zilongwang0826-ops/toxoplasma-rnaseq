"""Restartable HTTPS downloads with archive MD5 and size validation."""
import argparse
import hashlib
import os
import subprocess
import time
from pathlib import Path
from urllib.parse import urlparse
from common import sha256, utc_now, write_json


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def transfer(url, part, max_attempts=4, expected_bytes=None):
    """Retry transient transport failures, resuming the current partial each time."""
    transient_codes = {5, 6, 7, 18, 28, 35, 52, 55, 56, 92}
    transient_http = {408, 429, 500, 502, 503, 504}
    if max_attempts < 1:
        raise ValueError("At least one transfer attempt is required")
    for attempt in range(max_attempts):
        command = ["curl", "--fail", "--location", "--retry", "0",
                   "--connect-timeout", "30", "--speed-limit", "1024", "--speed-time", "120",
                   "--proto", "=https", "--continue-at", "-", "--output", str(part),
                   "--write-out", "%{http_code}", url]
        result = subprocess.run(command, check=False, stdout=subprocess.PIPE, text=True)
        if expected_bytes and part.exists():
            size = part.stat().st_size
            if size > expected_bytes:
                raise ValueError("Partial download exceeds expected archive size; preserving it")
            if size == expected_bytes:
                # A socket can close after the final byte; fetch still verifies MD5.
                return attempt + 1
        if result.returncode == 0:
            return attempt + 1
        status = result.stdout.strip()
        retryable = result.returncode in transient_codes or (result.returncode == 22 and status.isdigit() and int(status) in transient_http)
        if not retryable or attempt + 1 == max_attempts:
            raise subprocess.CalledProcessError(result.returncode, command)
        # Public transfer diagnostics only; no sequence content is logged.
        print(f"Transient transfer failure (curl={result.returncode}, HTTP={status}); resuming attempt {attempt + 2}/{max_attempts}", flush=True)
        time.sleep(5 * (2 ** attempt))


def fetch(url, output, expected_md5=None, expected_bytes=None):
    if urlparse(url).scheme != "https":
        raise ValueError("Only HTTPS sources are accepted")
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    part = output.with_name(output.name + ".part")
    if output.exists():
        if not expected_md5:
            receipt = Path(str(output) + ".source.json")
            if not receipt.exists():
                raise ValueError("Existing reference lacks provenance; refusing to overwrite")
            import json
            prior = json.loads(receipt.read_text())
            if prior["url"] != url or prior["sha256"] != sha256(output):
                raise ValueError("Existing reference provenance mismatch")
            return
        if md5(output) != expected_md5 or (expected_bytes and output.stat().st_size != expected_bytes):
            raise ValueError("Existing input checksum/size mismatch; preserving it for investigation")
        receipt = Path(str(output) + ".source.json")
        if not receipt.exists():
            write_json(receipt, {"url":url, "validated_at":utc_now(), "bytes":output.stat().st_size,
                "md5":expected_md5, "source_md5":expected_md5, "sha256":sha256(output),
                "validation":"Existing file revalidated against archive MD5/size; retrieval date unknown"})
        return
    if expected_bytes and part.exists() and part.stat().st_size > expected_bytes:
        raise ValueError("Partial download exceeds expected archive size; preserving it")
    already_complete = bool(expected_bytes and part.exists() and part.stat().st_size == expected_bytes)
    attempts = 0 if already_complete else transfer(url, part, expected_bytes=expected_bytes)
    digest = md5(part)
    if expected_md5 and digest != expected_md5:
        raise ValueError("MD5 mismatch; retained .part for investigation")
    if expected_bytes and part.stat().st_size != expected_bytes:
        raise ValueError("Byte-count mismatch; retained .part for investigation")
    receipt = {"url": url, "retrieved_at": utc_now(), "bytes": part.stat().st_size,
               "md5": digest, "source_md5": expected_md5, "sha256": sha256(part),
               "transfer_attempts_this_invocation": attempts,
               "validation": "source MD5 and size" if expected_md5 else "local SHA256; upstream SHA256 not supplied"}
    # Publish the receipt first: a crash leaves an incomplete file, never a trusted unrecorded input.
    write_json(str(output) + ".source.json", receipt)
    os.replace(part, output)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--url", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--md5")
    p.add_argument("--bytes", type=int)
    a = p.parse_args()
    fetch(a.url, a.output, a.md5, a.bytes)
