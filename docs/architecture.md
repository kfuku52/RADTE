# RADTE architecture

RADTE has two source representations with a one-way generation boundary:

```text
R/ version + cli + tree_navigation + input_parsers + calibration
  + chronos_backend + output + mcmctree_backend + generax
  + preparation + dating + run_output + main
                         |
                         v  tools/build_radte.R
                  executable radte.r
```

`R/load.R` defines the authoritative module order. Tests source these modules directly for speed and isolation. The generated `radte.r` is still exercised as a real subprocess so direct-download and Bioconda users receive the same behavior.

## Responsibilities

- `R/cli.R`: declarative option schema, coercion, validation, and generated help.
- `R/tree_navigation.R`: node lookup, canonical clade matching, and linear ancestor-context propagation.
- `R/input_parsers.R`: tree, NOTUNG, species parser/map, and age-bound validation.
- `R/calibration.R`: immutable calibration feasibility, dated-tree validation, node/age transfer, and fitted-age metadata.
- `R/chronos_backend.R`: retry budgets, numerical safeguards, seed inheritance, and feasible deterministic construction.
- `R/mcmctree_backend.R`: mirror constraints, staging, control generation, external execution, and posterior validation.
- `R/output.R`: atomic artifact writes, run transactions and prefix locks, plots, SHA-256 calculation, and manifest serialization.
- `R/generax.R`: NHX normalization and import.
- `R/preparation.R`: species/gene contexts, common format-independent age lookup, and calibration selection.
- `R/dating.R`: backend and fallback strategy selection.
- `R/run_output.R`: fitted-age tables, plots, effective configuration manifest, and complete-run publication.
- `R/main.R`: a short entrypoint connecting those independent phases.

## Runtime isolation and reproducibility

User-facing artifacts are written under `--outdir` with `--prefix`. All artifacts are staged before publication. A prefix lock excludes concurrent writers, previous artifacts are backed up, and the manifest is published last with rollback on errors. MCMCTree receives a unique work directory unless the caller explicitly supplies one.

`--seed` initializes chronos retry seeds and the MCMCTree control file. A successful run records effective options including defaults, backend parameters, RNG kind, calibration policy, tool versions, executable identity, timestamps, run identity, and input/output SHA-256 hashes in `<prefix>_run_manifest.tsv`.

## Performance invariants

Clades are represented by interned child-subtree IDs rather than repeatedly concatenated descendant-tip strings. Ancestor constraint minima are propagated once in preorder rather than recomputing every ancestor chain for every calibrated node. MCMCTree Newick generation uses an explicit stack and token buffer. Calibration feasibility uses two tree passes; original bounds are never rewritten. Correctness tests include deep trees, and `benchmark/benchmark_scaling.R` reports warmed repeated timings, allocations, and output equivalence checksums.

## CI layers

A base-R preflight checks generated artifacts/version consistency and selects lanes before dependency installation. Code pull requests run quality and current-R full/PAML tests. Code pushes retain Linux R 4.3–4.6 compatibility; relevant runtime/platform changes, releases, schedules, and manual runs include macOS and Windows. Documentation-only changes stop after preflight. Superseded PR/development runs are cancelled, while release pushes have distinct concurrency groups. Weekly runs install fresh dependencies without cache restoration. Cache save/restore uses the actual library exported by `setup-r`, with dependency-only fingerprints and seven-day coverage artifact retention.
