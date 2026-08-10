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

Every change merged to `master` must increase `R/version.R` and the matching `DESCRIPTION` version. `Rscript tools/check_version.R --against=origin/master` performs the same comparison as pull-request CI.

## Test boundaries

- Unit and validation tests source the canonical modules directly.
- Runtime-contract and integration tests execute the generated `radte.r` bundle.
- External PAML tests are opt-in locally and mandatory in the current-R integration lane.
- Example fixtures are immutable inputs; tests write only to temporary output directories.

Keep dependency declarations honest: add a direct dependency only when it is used in source, configuration, tests, or documentation. The compatibility installer contains one demonstrated R 4.3 workaround for `treeio`/`tidytree`; do not turn it into a general version-pinning policy.
