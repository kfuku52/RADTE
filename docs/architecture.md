# RADTE architecture

RADTE has two source representations with a one-way generation boundary:

```text
R/version.R + R/cli.R + R/tree_navigation.R + R/input_parsers.R
  + R/chronos_backend.R + R/output.R + R/mcmctree_backend.R
  + R/generax.R + R/main.R
                         |
                         v  tools/build_radte.R
                  executable radte.r
```

`R/load.R` defines the authoritative module order. Tests source these modules directly for speed and isolation. The generated `radte.r` is still exercised as a real subprocess so direct-download and Bioconda users receive the same behavior.

## Responsibilities

- `R/cli.R`: declarative option schema, coercion, validation, and generated help.
- `R/tree_navigation.R`: node lookup, canonical clade matching, and linear ancestor-context propagation.
- `R/input_parsers.R`: tree, NOTUNG, species parser/map, and age-bound validation.
- `R/chronos_backend.R`: constraint preparation, retry budgets, numerical safeguards, and deterministic fallback.
- `R/mcmctree_backend.R`: mirror constraints, staging, control generation, external execution, and posterior validation.
- `R/output.R`: atomic artifact writes, plots, SHA-256 calculation, and manifest serialization.
- `R/generax.R`: NHX normalization and import.
- `R/main.R`: orchestration only: mode selection, inputs, backend invocation, outputs, and manifest.

## Runtime isolation and reproducibility

User-facing artifacts are written under `--outdir` with `--prefix`. Each completed artifact replaces its destination through a same-directory temporary file, with rollback protection when an old file exists. MCMCTree receives a unique work directory unless the caller explicitly supplies one.

`--seed` initializes chronos retry seeds and the MCMCTree control file. A successful run records options, effective backend parameters, tool versions, executable identity, timestamps, and input SHA-256 hashes in `<prefix>_run_manifest.tsv`.

## Performance invariants

Clades are represented by interned child-subtree IDs rather than repeatedly concatenated descendant-tip strings. Ancestor constraint minima are propagated once in preorder rather than recomputing every ancestor chain for every calibrated node. The scaling test and `benchmark/benchmark_scaling.R` guard these paths.

## CI layers

Pull requests run a cached quality lane plus one current-R full/PAML lane. Pushes and the weekly schedule run the full Linux R 4.3–4.6, macOS, Windows, and PAML compatibility matrix. Generated bundle/help checks and version consistency fail before the expensive lanes can publish a release.
