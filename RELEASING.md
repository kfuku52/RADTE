# Releasing RADTE

Keep `radte_version` in `R/version.R` on the next semantic version for every
change merged to `master`, and keep the `DESCRIPTION` version identical. Run
`make build`, `make cli-docs`, and `Rscript tools/check_version.R` before the
change is merged. Never edit the generated `radte.r` version directly.

The `RADTE ci` workflow validates each push. After it succeeds, the release
workflow checks `R/version.R` from the exact tested commit. Pull-request CI
also verifies that the version increased relative to its base branch. Only a successful
`master` push from this repository can trigger a release; pull requests, fork
runs, and manual workflow dispatches cannot publish a tag:

- Versions whose patch component is nonzero (for example, `0.3.4`) remain
  available from `master`, but do not receive a Git tag or GitHub Release.
- Major and minor versions whose patch component is zero (for example,
  `0.4.0` or `1.0.0`) receive an annotated `v<version>` tag and a GitHub
  Release automatically.

Bioconda discovers tagged upstream releases. Consequently, its `radte` recipe
is updated for major and minor releases only; patch-only versions are
intentionally not autobumped.

Do not create release tags manually unless recovering the automated workflow.
If recovery is necessary, point the annotated tag at the commit that passed
`RADTE ci` and preserve the existing `v<version>` tag format.
