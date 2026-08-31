# Installation and verification

These instructions describe the source version in this checkout. Commands assume the repository root and a POSIX shell unless marked otherwise. [Bioconda releases](#bioconda-packaged-release) may contain an older CLI.

## Source version

Install [R](https://www.r-project.org/) 4.3 or later and make `Rscript` available in your terminal. RADTE's CI tests R 4.3–4.6. Obtain the source, then install the dependency needed for NOTUNG input with the chronos backend:

```sh
git clone https://github.com/kfuku52/RADTE
cd RADTE
Rscript -e 'install.packages("ape", repos="https://cloud.r-project.org")'
./radte.r --version
```

If R asks for a personal package library, use a writable location. Do not install into a read-only system library. `--version` and `--help` intentionally work without loading ape or treeio, so they verify CLI startup, not the full installation.

### GeneRax input

`--generax_nhx` additionally needs treeio. Use the repository's compatibility installer:

```sh
Rscript test/install_dependencies.R
```

This installs ape, treeio, and the installer/test support packages BiocManager and testthat. It selects the Bioconductor release for the running R version and includes the documented R 4.3 treeio/tidytree compatibility workaround. It does not install the external NOTUNG, GeneRax, or PAML programs. NOTUNG and GeneRax are needed to create reconciliations for your own data; the bundled examples already contain their outputs.

### MCMCTree input

Install [PAML](https://github.com/abacus-gene/paml) separately for `--dating_backend=mcmctree`. Either put its `mcmctree` executable in `PATH` or pass `--mcmctree_bin=/absolute/path/to/mcmctree` (use `mcmctree.exe` on Windows). To check a `PATH` installation without starting an analysis:

```sh
Rscript -e 'stopifnot(nzchar(Sys.which("mcmctree"))); cat(Sys.which("mcmctree"), "\n")'
```

An alignment whose taxon names exactly match the gene-tree tips is also required. The repository does not include alignments for the illustrated gene families. See the [MCMCTree templates](usage.md#example-4-radte-with-the-mcmctree-backend), and assess MCMC convergence separately from RADTE's structural validation.

### Windows invocation

Use `Rscript radte.r` instead of `./radte.r`. For example, this single-line command works from the repository root once ape is installed:

```powershell
Rscript radte.r --species_tree=data/example_notung_01/species_tree.nwk --gene_tree=data/example_notung_01/gene_tree.nwk.reconciled --notung_parsable=data/example_notung_01/gene_tree.nwk.reconciled.parsable.txt --max_age=1000 --chronos_model=discrete --chronos_lambda=1 --outdir=results/install-check --prefix=example --seed=1
```

The multi-line examples elsewhere use POSIX `\` continuations; do not copy those continuations into PowerShell. R package source installation may require the matching [Rtools](https://cran.r-project.org/bin/windows/Rtools/). The [chronos timeout limitations](usage.md#--chronos_attempt_timeout_sec) on Windows also apply.

## Verify a source installation

The [README Quick start](../README.md#quick-start) exercises input parsing, chronos, tables, plots, and manifest publication. It uses a bundled NOTUNG reconciliation, so Java and NOTUNG are not needed for that check. Confirm that it completes and creates `results/example_gene_tree_output.nwk` and `results/example_run_manifest.tsv`.

Run [example 2](usage.md#example-2-radte-after-generax) to verify treeio/NHX input and [example 3](usage.md#example-3-transfer-species-tree-node-age-bounds) to verify interval calibrations. All examples write to separate output directories and leave the input fixtures unchanged.

## Bioconda packaged release

Follow the [Bioconda setup guidance](https://bioconda.github.io/#usage). With conda available, create an environment with explicit channel order and strict priority:

```sh
conda create -n radte --strict-channel-priority -c conda-forge -c bioconda radte
conda activate radte
conda list radte
```

Availability depends on the platform and on the dependencies of the packaged release. Use the source instructions above on platforms without a suitable package.

The packaged executable is `radte`, whereas this repository's script is `./radte.r` (or `Rscript radte.r`). Check the installed version with `conda list radte` and read the README/documentation at the corresponding [upstream tag](https://github.com/kfuku52/RADTE/tags). Older releases may not support options such as `--outdir`, `--prefix`, `--seed`, `--help`, or `--version`. The latest source examples are not automatically compatible with an older package; use the source installation when you need features absent from that release. Installing the package also does not provide a repository checkout at your current directory.

## Development

See [CONTRIBUTING.md](../CONTRIBUTING.md) for development dependencies, the R 4.4.3 recorded library, and checks. Runtime dependencies, reproducibility locks, and external system tools are managed separately.
