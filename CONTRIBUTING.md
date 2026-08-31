# Contributing to RADTE

RADTE keeps maintainable modules and a portable one-file CLI in the same repository. Edit files under `R/`; do not edit generated sections of `radte.r`.

## Local setup

For runtime prerequisites and platform-specific commands, see [installation](docs/installation.md). Install normal and development dependencies from the repository root:

```sh
Rscript -e "Sys.setenv(RADTE_INSTALL_DEV='true'); source('test/install_dependencies.R')"
```

### Recorded dependency environment

`renv.lock` records R **4.4.3**, Bioconductor **3.20**, and specific R package versions. Install R 4.4.3 separately before reproducing this environment; renv does not install R, PAML, compilers, or system libraries. Compatibility work should still use the R 4.3–4.6 CI matrix; the lockfile is a separate reproducibility reference.

This repository does not activate renv automatically. In a POSIX shell, select a local library explicitly before restoring or running checks:

```sh
mkdir -p .r-lib
export R_LIBS_USER="$PWD/.r-lib"
Rscript -e "stopifnot(getRversion() == '4.4.3'); install.packages('renv', repos='https://cloud.r-project.org')"
Rscript -e "renv::restore(project='.', library=Sys.getenv('R_LIBS_USER'), prompt=FALSE)"
Rscript -e "stopifnot(normalizePath(.libPaths()[1]) == normalizePath(Sys.getenv('R_LIBS_USER'))); renv::status(project='.', library=.libPaths())"
make check
```

In PowerShell, replace the first two lines with `New-Item -ItemType Directory -Force .r-lib` and `$env:R_LIBS_USER = (Resolve-Path .r-lib).Path`; use the same R commands afterwards. Windows development also needs build tools and `make` (for example those supplied by the matching [Rtools](https://cran.r-project.org/bin/windows/Rtools/)).

Keep `R_LIBS_USER` set to this absolute path in every terminal that should use the recorded packages; a new terminal needs the setting again. The check above runs in a fresh R process to verify the selected library and all libraries R searches. Restore can reuse matching package versions already installed in R's system library. `.r-lib/` is ignored by Git. This explicit-library workflow does not create an `.Rprofile` or switch other sessions to renv. A bare `renv::restore()` in an unactivated project can instead target the current user library. If you prefer renv auto-activation, follow its [activation instructions](https://rstudio.github.io/renv/reference/activate.html) before restoring, and keep that setup separate from unrelated changes.

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
