# Observed resource validation

## Initial full-genome-decoy index attempt

The real mouse reference index used Salmon 1.10.3, two threads, k=31, sparse indexing, retained duplicate transcript targets and 66 primary-assembly decoys. The target/reference configuration was unchanged during this attempt.

The Docker limit was `18000m`, meaning 18,874,368,000 bytes (17.578 GiB). At 00:51:49 UTC on 28 October 2026, Docker recorded an `oom` event for the index container, followed by exit 137. The failure occurred during presence-vector construction after the minimal perfect hash had been built. Snakemake elapsed time was 2:42:31. The exact peak allocation was not captured; periodic memory samples are not peak measurements.

All original sequencing files, downloaded reference sources and the prepared reference were retained. Snakemake's former default failure cleanup removed the incomplete index directory before the monitoring check. Logs and Docker event evidence were preserved. The local profile now uses `keep-incomplete: true` so future failed outputs remain available for investigation; preserve or archive these outputs deliberately before any retry.

## Bounded retry

The approved retry limit is `22888m`: 23,999,807,488 bytes (22.352 GiB), strictly below 24 decimal GB. The configuration retains the legacy `index_mem_mb` name; the wrapper passes that value to Docker's `m` suffix, whose units are MiB. The local scheduling budget is 24,000 MB and one job. Two threads and every scientific indexing parameter remain unchanged.

The VM reported 32,454,103,040 bytes of available memory after the failed container exited. The higher cap therefore leaves headroom for the host. It is an evidence-based retry, not proof that the build will fit or complete. No further increase beyond the approved bound is authorized automatically. Reuse the existing protected input files; do not download them again or remove decoys to conceal a resource failure.

## Upstream memory consideration

The exact [pufferfish dependency for Salmon 1.10.3](https://github.com/COMBINE-lab/salmon/blob/v1.10.3/scripts/fetchPufferfish.sh) retains a `std::unordered_set<uint64_t>` of sampled positions during [sparse presence-vector construction](https://github.com/COMBINE-lab/pufferfish/blob/salmon-v1.10.3/src/PufferfishIndexer.cpp#L799-L1036). The inspected code uses this set for a diagnostic size message; its nodes, buckets and allocator overhead add to the fixed sparse arrays. The observed reference has 268,781,739 sampled positions, so sparse index file size alone is not an adequate peak-RAM estimate. No source patch or alternative indexing method has been applied. If the bounded retry fails, a separately reviewed software fix or larger execution allocation must be considered; success is not assumed.

## Completed index and quantification load adjustment

The bounded retry completed at 05:11:27 UTC on 28 October 2026. Docker recorded index-container exit 0, and Snakemake marked the index rule complete at 05:11:28. The sparse index occupies approximately 9.9 GiB and contains 2,275,675,295 k-mers, 268,781,739 sampled k-mers and 66 decoys. This proves index construction completed; successful loading and biological quantification are separate checks.

C111 quantification then exceeded its initial `8000m` cap (8,388,608,000 bytes) while loading sparse-index vectors, before read mapping began. Docker recorded an `oom` event and exit 137. The index's component sizes make the initial cap insufficient; this was a resource estimate failure, not evidence of damaged input reads.

Quantification now uses `16000m` (16,777,216,000 bytes; 15.625 GiB), two threads and one concurrent job. This remains below the approved single-job resource envelope and leaves substantial host memory headroom. No scientific quantification parameter, reference or software version changed. The completed index and verified reads are reused; only the failed quantification is retried. Its partial directory and logs are retained under a named failed attempt. Actual mapping and peak usage remain to be measured.

The C111 retry successfully loaded the completed index and had processed at least 9.5 million fragments by the 05:29 UTC check. Reported container memory was 11.76 GiB under the 15.625 GiB cap. Quantification completion and final mapping/QC metrics were still pending.

## Complete-cohort outcome

All 12 libraries subsequently completed with the same index, two threads and 16,000 MiB quantification cap. Host assignment ranged 73.46–88.44%; all processed-pair totals matched their archive counts. No additional memory increase or reference change was needed. These successful executions resolve the load/completion uncertainty above; the exact peak memory remains unmeasured because periodic observations do not establish a peak. See `validation-status.md` for final evidence.
