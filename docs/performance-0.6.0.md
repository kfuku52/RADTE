# 0.6.0 validation and performance measurements

Local comparison on 2026-08-31: macOS arm64, R 4.4.3, ape 5.8.1. Baseline is `fb1ebc8` (0.5.1); the updated source is 0.6.0. Measurements ran sequentially, without other RADTE validation processes. They are observations on one machine, not CI timing thresholds.

Each scaling workload used one warm-up and five measured repetitions. Inputs were left/caterpillar trees with 250, 500, 1,000, 2,000, and 4,000 tips. Calibration intervals were strictly feasible so output equality was meaningful. The benchmark also measures clade matching, conflict scanning, and gene/species mapping.

| Workload | Tips | Before median | After median |
| --- | ---: | ---: | ---: |
| Constraint validation | 1,000 | 0.046 s | 0.001 s |
| Constraint validation | 2,000 | 0.096 s | 0.003 s |
| Constraint validation | 4,000 | 0.202 s | 0.005 s |
| MCMCTree writer | 1,000 | 0.086 s | 0.007 s |
| MCMCTree writer | 2,000 | C stack error | 0.015 s |
| MCMCTree writer | 4,000 | C stack error | 0.030 s |

All 23 workloads that completed in both versions had identical serialized output checksums. Both previously failing writer sizes now complete. This does not imply that a particular PAML executable supports trees of this size; PAML's compiled limits are separate from RADTE's writer.

At 1,000 tips the writer allocated 32,521,936 bytes before and 684,072 bytes after. The new two-pass constraint validator uses additional linear-size arrays (801,296 allocated bytes at 4,000 tips versus 160,128 previously). Whole-benchmark peak RSS was 279,986,176 bytes before and 282,836,992 after. Allocated bytes are cumulative R allocations from a separate `Rprofmem()` pass; they are not peak memory.

`make check` took 20.17 s before and 20.55 s after, with peak RSS 501,121,024 and 438,353,920 bytes respectively. The new suite includes additional regression properties, so this comparison does **not** establish faster total local checks. It demonstrates approximately unchanged wall time with broader checks and less peak memory in this run. The CI cache repair targets the much larger repeated dependency-install cost and must be verified through cache save/restore logs on GitHub.

Raw measurements: [before](../benchmark/results/0.6.0/before.tsv), [after](../benchmark/results/0.6.0/after.tsv). Reproduce with:

```sh
Rscript benchmark/benchmark_scaling.R --root=/path/to/baseline --output=before.tsv
Rscript benchmark/benchmark_scaling.R --output=after.tsv
```

Use `/usr/bin/time -l` on macOS or `/usr/bin/time -v` on Linux for whole-process peak RSS. Regression tests separately cover infeasible and overlapping bounds, finite-timeout RNG inheritance, equivalent reconciliation formats, output node identity, protected inputs, external-tool failures, output rollback, and the installed package API.
