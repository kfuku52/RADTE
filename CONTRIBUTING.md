# Contributing to RADTE

RADTE keeps maintainable modules and a portable one-file CLI in the same repository. Edit files under `R/`; do not edit generated sections of `radte.r`.

## Local setup

Install normal and development dependencies:

```sh
RADTE_INSTALL_DEV=true Rscript test/install_dependencies.R
```

For the exact R 4.4 development environment, install `renv` and run `renv::restore()`. Compatibility work should still be tested through the same R 4.3–4.6 matrix used by CI; the lockfile is not a substitute for that matrix.

## Change loop

```sh
make build
make cli-docs
make check
make test
```

Run `make paml` when MCMCTree behavior changes and `make benchmark` when tree traversal or calibration propagation changes. `make coverage` writes local reports that match the CI artifact.

Before pushing changes to `master`, increase `R/version.R` and the matching `DESCRIPTION` version; individual unpushed commits need no bump. `Rscript tools/check_version.R --against=origin/master` performs the same comparison as pull-request CI.

## Test boundaries

- Unit and validation tests source the canonical modules directly.
- Runtime-contract and integration tests execute the generated `radte.r` bundle.
- `tests/smoke.R` exercises the installed `RADTE::radte_main()` API during `R CMD check`.
- External PAML tests are opt-in locally and mandatory in the current-R integration lane.
- Example fixtures are immutable inputs; tests write only to temporary output directories.

Keep dependency declarations honest: add a direct dependency only when it is used in source, configuration, tests, or documentation. The compatibility installer contains one demonstrated R 4.3 workaround for `treeio`/`tidytree`; do not turn it into a general version-pinning policy.

## Coverage and performance

Fast tests use a single `test_dir()` selection, so helpers and fixtures are loaded once. Coverage includes `test-tree-manipulation.R`; mocks bind at the dependency environment, not an unrelated global binding. CLI subprocess integration remains a separate boundary from in-process coverage.

Benchmarks warm up each workload and report the median/minimum of five runs, allocated bytes when R memory profiling is available, and output checksums. Compare a clean baseline with `Rscript benchmark/benchmark_scaling.R --root=/path/to/baseline --output=before.tsv`, then run against the changed tree. Use `/usr/bin/time -l` on macOS or `/usr/bin/time -v` on Linux to record whole-process peak RSS separately. CI retains deep-tree correctness checks without a single-machine time threshold.

CI caches the `R_LIBS_USER` established by `setup-r`. Cache keys use dependency declarations and installer content, not the project version. Weekly runs bypass caches to exercise fresh compatible dependencies. Documentation-only changes run preflight checks; runtime changes, releases, schedules, and manual runs retain the relevant R/OS/PAML coverage.
