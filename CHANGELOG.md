# Changelog

## 0.6.0

- Preserve original chronos age bounds, reject infeasible constraints, and validate dated ages, topology, tip sets, and positive branches. Retries no longer widen or soften hard constraints.
- Fix finite-timeout seed reproducibility while restoring the caller RNG state; record effective defaults and RNG kind.
- Normalize NOTUNG and GeneRax age mapping, reject unknown species nodes, preserve quoted/numeric-looking IDs, and correct the bundled NOTUNG species-node aliases.
- Preserve gene-tree node numbering in allS and MCMCTree outputs and reject discordant allS topologies.
- Protect MCMCTree inputs before workdir cleanup, use copied alignment snapshots, support relative executables, validate posterior structure, and generate deep-tree Newick iteratively.
- Publish complete output generations with prefix locks, rollback, manifest-last ordering, run IDs, and artifact hashes. Add fitted ages and constraint status to calibration tables.
- Resolve configuration from one schema and separate input, calibration, inference, and output phases.
- Fix final-attempt dependency installation, cache the actual R library in CI, gate documentation/platform lanes, and preserve release runs while cancelling superseded development runs.
- Restore numerical/retry coverage, load fast-test helpers once, add an installed-package smoke test and regression properties, and replace timing thresholds with warmed repeated benchmarks.

See [migration and run contracts](docs/run-contracts.md) for intentional behavior changes.

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
