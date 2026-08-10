# RADTE command-line reference

This file is generated from `R/cli.R`. Run `make cli-docs` after changing the option schema.

```text
RADTE 0.5.0

Usage:
  radte --species_tree=FILE (--notung_parsable=FILE --gene_tree=FILE | --generax_nhx=FILE) [options]

Meta options:
  --help, -h       Show this help and exit.
  --version, -V    Show the RADTE version and exit.

Options use --key=value syntax:
  --allow_constraint_drop=VALUE
      Allow S/R constraint-drop fallback stages. [default: true; backend: chronos].
  --chronos_attempt_timeout_sec=VALUE
      Wall-time limit per chronos attempt; 0/inf/none/off disables it. [default: 60; backend: chronos].
  --chronos_dual_iter_max=VALUE
      Initial chronos dual.iter.max. [default: 20; backend: chronos].
  --chronos_eval_max=VALUE
      Initial chronos eval.max. [default: 250; backend: chronos].
  --chronos_high_control_fallback=VALUE
      Enable the high-cost chronos control fallback. [default: true; backend: chronos].
  --chronos_iter_max=VALUE
      Initial chronos iter.max. [default: 250; backend: chronos].
  --chronos_lambda=VALUE
      Chronos smoothing parameter (required for chronos). [default: required/conditional; backend: chronos].
  --chronos_model=VALUE
      Chronos model (required for chronos). [default: required/conditional; backend: chronos]. Choices: discrete, relaxed, correlated.
  --chronos_total_timeout_sec=VALUE
      Total wall-time budget for chronos retries; 0/inf/none/off disables it. [default: 300; backend: chronos].
  --dating_backend=VALUE
      Dating engine. [default: chronos; backend: all]. Choices: chronos, mcmctree.
  --gene_tree=VALUE
      Rooted Newick gene tree (required in NOTUNG mode). [default: required/conditional; backend: all].
  --generax_nhx=VALUE
      GeneRax NHX reconciliation input. [default: required/conditional; backend: all].
  --max_age=VALUE
      Upper age for duplication nodes above the species-tree root. [default: required/conditional; backend: all].
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
      MCMCTree alignment file (required for mcmctree). [default: required/conditional; backend: mcmctree].
  --mcmctree_seqtype=VALUE
      MCMCTree sequence type. [default: 0; backend: mcmctree].
  --mcmctree_timeout_sec=VALUE
      Wall-time limit for the external MCMCTree process. [default: inf; backend: mcmctree].
  --mcmctree_usedata=VALUE
      MCMCTree usedata value; RADTE currently supports only 1. [default: 1; backend: mcmctree].
  --mcmctree_workdir=VALUE
      Explicit MCMCTree staging directory; default is a unique directory under outdir. [default: required/conditional; backend: mcmctree].
  --notung_parsable=VALUE
      NOTUNG parsable reconciliation input. [default: required/conditional; backend: all].
  --outdir=VALUE
      Directory for RADTE output artifacts. [default: .; backend: all].
  --pad_short_edge=VALUE
      Minimum dated branch length. [default: required/conditional; backend: all].
  --prefix=VALUE
      Output filename prefix. [default: radte; backend: all].
  --seed=VALUE
      Base random seed for chronos retries and MCMCTree. [default: 1; backend: all].
  --species_map_tsv=VALUE
      Species mapping TSV for --species_parser=map. [default: required/conditional; backend: all].
  --species_node_bounds_tsv=VALUE
      Species-tree node age bounds TSV. [default: required/conditional; backend: all].
  --species_parser=VALUE
      Species-label parser. [default: legacy; backend: all]. Choices: legacy, taxonomic, regex, map.
  --species_regex=VALUE
      Species extraction regex for --species_parser=regex. [default: required/conditional; backend: all].
  --species_tree=VALUE
      Rooted, dated, fully bifurcating species tree. [default: required/conditional; backend: all].
```
