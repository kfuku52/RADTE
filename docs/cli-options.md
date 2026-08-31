# RADTE command-line reference

This file is generated from `R/cli.R`. Run `make cli-docs` after changing the option schema.

```text
RADTE 0.6.3

Usage:
  radte --species_tree=FILE (--notung_parsable=FILE --gene_tree=FILE | --generax_nhx=FILE) [options]

Meta options:
  --help, -h       Show this help and exit.
  --version, -V    Show the RADTE version and exit.

Options use --key=value syntax:
  --allow_constraint_drop=VALUE
      Allow RS -> S -> R fallback within each control profile; retained bounds are never relaxed. [default: true; backend: chronos].
  --chronos_attempt_timeout_sec=VALUE
      Per-attempt budget; enforced by child termination on Unix, R interrupt checks on Windows. 0/inf/none/off disables it. [default: 60; backend: chronos].
  --chronos_dual_iter_max=VALUE
      Initial chronos dual.iter.max. [default: 20; backend: chronos].
  --chronos_eval_max=VALUE
      Initial chronos eval.max. [default: 250; backend: chronos].
  --chronos_high_control_fallback=VALUE
      Enable the high-cost chronos control fallback. [default: true; backend: chronos].
  --chronos_iter_max=VALUE
      Initial chronos iter.max. [default: 250; backend: chronos].
  --chronos_lambda=VALUE
      Chronos smoothing parameter. [default: not set; required: chronos; backend: chronos].
  --chronos_model=VALUE
      Chronos model. [default: not set; required: chronos; backend: chronos]. Choices: discrete, relaxed, correlated.
  --chronos_total_timeout_sec=VALUE
      Total chronos retry budget; Windows attempts may overrun it. 0/inf/none/off disables it. [default: 300; backend: chronos].
  --dating_backend=VALUE
      Dating engine. [default: chronos; backend: all]. Choices: chronos, mcmctree.
  --gene_tree=VALUE
      Rooted Newick gene tree. [default: not set; required: notung; backend: all].
  --generax_nhx=VALUE
      GeneRax NHX reconciliation input; choose this or --notung_parsable. [default: not set; backend: all].
  --max_age=VALUE
      Upper age for duplication nodes above the species-tree root. [default: not set; required: always; backend: all].
  --mcmctree_bin=VALUE
      MCMCTree executable name or path. [default: mcmctree; backend: mcmctree].
  --mcmctree_burnin=VALUE
      MCMCTree burn-in samples. [default: 2000; backend: mcmctree].
  --mcmctree_clock=VALUE
      MCMCTree clock model. [default: 2; backend: mcmctree].
  --mcmctree_model=VALUE
      MCMCTree substitution model. [default: 0; backend: mcmctree].
  --mcmctree_ncatG=VALUE
      MCMCTree number of gamma categories. [default: 5; backend: mcmctree].
  --mcmctree_nsample=VALUE
      MCMCTree number of posterior samples. [default: 20000; backend: mcmctree].
  --mcmctree_sampfreq=VALUE
      MCMCTree sampling frequency. [default: 10; backend: mcmctree].
  --mcmctree_seqfile=VALUE
      MCMCTree alignment file. [default: not set; required: mcmctree; backend: mcmctree].
  --mcmctree_seqtype=VALUE
      MCMCTree sequence type. [default: 0; backend: mcmctree].
  --mcmctree_timeout_sec=VALUE
      Wall-time limit for the external MCMCTree process. [default: inf; backend: mcmctree].
  --mcmctree_usedata=VALUE
      MCMCTree usedata value; RADTE currently supports only 1. [default: 1; backend: mcmctree].
  --mcmctree_workdir=VALUE
      Explicit MCMCTree staging directory; default is a unique directory under outdir. [default: not set; backend: mcmctree].
  --notung_parsable=VALUE
      NOTUNG parsable reconciliation input; choose this or --generax_nhx. [default: not set; backend: all].
  --outdir=VALUE
      Directory for RADTE output artifacts. [default: .; backend: all].
  --pad_short_edge=VALUE
      Floor for input gene branches; also minimum chronos dated branches within original bounds. No MCMCTree output minimum. [default: not set; backend: all].
  --prefix=VALUE
      Output filename prefix. [default: radte; backend: all].
  --seed=VALUE
      Base random seed for chronos retries and MCMCTree. [default: 1; backend: all].
  --species_map_tsv=VALUE
      Species mapping TSV; required with --species_parser=map. [default: not set; backend: all].
  --species_node_bounds_tsv=VALUE
      Species-tree node age bounds TSV. [default: not set; backend: all].
  --species_parser=VALUE
      Species-label parser. [default: legacy; backend: all]. Choices: legacy, taxonomic, regex, map.
  --species_regex=VALUE
      Species extraction regex; required with --species_parser=regex. [default: not set; backend: all].
  --species_tree=VALUE
      Rooted, dated, fully bifurcating species tree. [default: not set; required: always; backend: all].
```
