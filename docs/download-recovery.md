# Interrupted HTTPS acquisition and recovery

The first full-cohort run stopped while acquiring T332 read 2 at 06:29:50 UTC on 28 October 2026. Curl reported exit 56 (`SSL_ERROR_SYSCALL`) after an interrupted TLS read. The partial file was retained at 1,076,799,745 bytes of the expected 3,339,038,092 bytes. This was a transport interruption; no MD5 failure or biological-data corruption was demonstrated, and the incomplete file was never published as a validated input.

The downloader now makes at most four transfer attempts per invocation, with 5/10/20-second backoff and a new `--continue-at -` operation for each attempt. Explicitly enumerated transient network failures and HTTP 408/429/500/502/503/504 can retry. Permanent HTTP failures, unsupported range requests, local write errors and checksum mismatches stop for investigation. TLS verification is not disabled. Curl's internal retries are disabled to keep the invocation's transfer budget explicit.

A complete-sized partial file, including one completed immediately before a socket error, is verified directly, without requesting an invalid end-of-file range. Oversized or checksum-mismatching partial files are preserved and not published. Final files and receipts become outputs only after archive MD5/size verification. Interrupted bytes remain available for resume; no automatic deletion or full restart is used.

Acquisition rules treat downloader/helper source timestamps as `ancient` inputs so a transport implementation repair does not force reacquisition of existing protected files. URL, expected MD5/bytes and container identity remain rule parameters; changes to source identity still require reviewed reconciliation. The exact executing script remains recorded by Git and source provenance. This exception applies only to acquisition helper timestamps, not analysis inputs or statistical outputs.

Nine targeted tests cover connection-reset continuation, transient HTTP retries, permanent failures, exhausted retry budgets, complete partial recovery, corrupt partials oversized partial retention, final-byte transport resets, and certificate/range/local-write failures. The real partial's final MD5 must still pass before it is accepted; synthetic retry tests do not establish the integrity of downloaded study data.

Recovery execution: the container tests passed all 29 Python cases, and the remote dry-run scheduled only 61 remaining jobs (21 FASTQ downloads); the completed index, three validated mate files and C111 quantification were reused. The restarted transfer explicitly resumed from byte 1,076,799,745 and the partial grew beyond that offset. The final archive checksum remains a required gate, not an assumption from successful range resumption.

## Verified recovery outcomes

The interrupted acquisitions subsequently completed with exact archive MD5 and size validation. All 24 final FASTQ files were independently rehashed during the final delivery audit.

| File | Expected/final bytes | Matched archive MD5 | Completed UTC |
| --- | ---: | --- | --- |
| T332 mate 2 | 3,339,038,092 | `cd8e7011652af8acb5ecb870cf678440` | 2026-10-28 07:14:49 |
| T331 mate 1 | 2,598,313,464 | `62d634dfe758350ae86135829ea9a6b1` | 2026-10-28 11:34:58 |
| T112 mate 2 | 3,029,940,866 | `d5043bd577c4439df34e44a12debe74e` | 2026-10-28 15:00:15 |

No whole-file restart, checksum override or disabling of TLS verification was used. These transport recoveries changed no biological processing parameters.
