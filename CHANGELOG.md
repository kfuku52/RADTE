# Changelog

## 0.5.1

- Install the libcurl development headers required by the quality job's R development dependencies.

## 0.5.0

- Split the canonical implementation into responsibility-focused `R/` modules while retaining a generated executable `radte.r` distribution.
- Added a typed CLI schema, `--help`, `--version`, `--outdir`, `--prefix`, `--seed`, and `--mcmctree_timeout_sec`.
- Added atomic outputs, per-run manifests with input SHA-256 hashes, unique MCMCTree work directories, executable/banner capture, and deterministic backend seeds.
- Replaced quadratic clade matching and ancestor-constraint scans with interned subtree IDs and preorder propagation.
- Accelerated validation tests by running application validation in-process while retaining generated-CLI integration coverage.
- Added performance/runtime-contract tests, a benchmark, dependency metadata and lockfile, lint/coverage tooling, cached layered CI, and version/generated-artifact gates.
- Separated immutable example inputs from reference output and removed obsolete debug artifacts.
