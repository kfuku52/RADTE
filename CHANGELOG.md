# Changelog

## 0.6.3

- Make the bundled usage examples runnable from the repository root with separate output directories, and label external-tool templates and historical reference outputs explicitly.
- Align calibration, retry-order, padding, timeout, and output descriptions with the current implementation; fix the mirror diagram link and replace the old ape documentation link.
- Add runtime installation and verification guidance and an explicit local-library workflow for restoring the R 4.4.3 development dependencies; include the missing Bioconductor version marker in the lockfile.
- Distinguish unset option defaults from required conditions in generated CLI help and documentation without changing argument validation or dating behavior.

## 0.6.2

- Keep documentation-only CI lightweight when the required patch version bump also updates metadata and the generated CLI. Dependency, runtime, executable-mode, and release changes still receive full validation.

## 0.6.1

- Preserve the executable link name when invoking PAML wrappers that select a program by their invocation name, including the Debian/Ubuntu package.
- Capture CLI validation errors portably in the infeasible-constraint regression test, including on Windows.
- Give the finite-timeout numerical reproducibility test scheduling headroom without changing production timeouts or dedicated timeout assertions.

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
