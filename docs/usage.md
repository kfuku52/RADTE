# Usage and examples

Commands assume the repository root unless an example changes directories. See [run contracts and migration notes](run-contracts.md) for calibration, reproducibility, and output guarantees.

## Options
Arguments use `--key=value` syntax. Unknown or duplicated option names are rejected before input processing.
Run `./radte.r --help` for the generated help, or see the complete [command-line reference](cli-options.md). `./radte.r --version` prints the version without loading R packages or input files.

#### `--outdir`, `--prefix`
All artifacts can be directed to one directory with `--outdir` (default: the current directory) and renamed with `--prefix` (default: `radte`). For example, `--outdir=results --prefix=family42` writes `results/family42_gene_tree_output.nwk` and the corresponding tables, plots, and manifest. RADTE creates the output directory when needed, stages the complete run, and publishes its manifest last. A failed writer leaves the previous completed run intact.

#### `--seed`
Base random seed for reproducible `chronos` retries and MCMCTree runs. The default is `1`. The requested and ultimately used seeds are recorded in the run manifest.

#### `--species_tree`
Species tree with estimated divergence time.
By default, leaves (species) should be labeled as `GENUS_SPECIES` (e.g., Homo_sapiens).
If `--species-parser=taxonomic` is used, taxonomically qualified labels such as `Dictyostelium_cf_discoideum` are also accepted.
The tree is expected to be ultrametric and branch lengths should represent evolutionary time (e.g., million years).
Internal nodes including the root node must be uniquely labeled and the same file should be consistently used for **NOTUNG/GeneRax** and **RADTE**.
Don't know how to label internal nodes? Try this R one-liner.
```
R -q -e "library(ape); t=read.tree('species_tree_noLabel.nwk'); \
t[['node.label']]=paste0('s',1:Nnode(t)); \
write.tree(t, 'species_tree.nwk')"
```
#### `--gene_tree`
Rooted newick tree. By default, leaves (genes) should be labeled as `GENUS_SPECIES_GENEID` (e.g., Homo_sapiens_ENSG00000102144).
If `--species-parser` is switched to `taxonomic`, `regex`, or `map`, gene tips are interpreted with that parser instead.
The tree is expected to be non-ultrametric and branch lengths should represent substitutions per site.
Use the tree that **NOTUNG** produces because its internal nodes are correctly labeled in accordance with the **NOTUNG parsable file**, another input for this program.
#### `--notung_parsable`
An output file from **NOTUNG** (tested with version 2.9) can be used to acquire the species–gene relationships in phylogeny reconciliation. See **Examples** for details.
#### `--generax_nhx`
Instead of the **NOTUNG** output, the NHX tree from **GeneRax** can also be used as an input. In GeneRax mode, `--gene_tree` is ignored and `--notung_parsable` must not be specified. See **Examples** for details.
The NHX species annotation tag `S` is required for all nodes and must match species-tree node labels.
The duplication tag `D` is optional (missing values are treated as non-duplication). Accepted duplication values are `Y`, `YES`, `TRUE`, `T`, `1`; accepted non-duplication values are `N`, `NO`, `FALSE`, `F`, `0`.
#### `--species-parser`
Optional species-label parser. Default is `legacy`.
Supported values are:
* `legacy`: current `GENUS_SPECIES` / `GENUS_SPECIES_GENEID` behavior.
* `taxonomic`: accepts taxonomically qualified species labels such as `Dictyostelium_cf_discoideum_gene123` or `Arabidopsis_thaliana_subsp_lyrata_gene456`.
* `regex`: extracts species labels from gene tips using `--species-regex`.
* `map`: resolves species labels from `--species-map-tsv`.
#### `--species-regex`
Required when `--species-parser=regex`.
RADTE uses the first capture group if the regex contains captures; otherwise it uses the full match as the extracted species label.
#### `--species-map-tsv`
Required when `--species-parser=map`.
This should be a tab-delimited file with `label` and `species` columns.
An optional `taxonomy_query` column can be used to override scientific-name conversion.
#### `--species_node_bounds_tsv`
Optional tab-delimited file for species-tree node age constraints.
This file should contain a `node` column plus either `age` or the pair `age_min` / `age_max`.
The node names must match the labeled internal/root nodes in `--species_tree`.
The point age implied by the species-tree branch lengths must fall within the supplied interval.
RADTE transfers these bounds to gene-tree speciation nodes. If the same species-tree node corresponds to multiple gene-tree speciation nodes, `mcmctree` enforces a shared age parameter, while `chronos` only reuses the same interval bounds.
Accepted examples:
```
node	age
n1	10
root	30
```
or
```
node	age_min	age_max
n1	8	12
root	28	32
```
Use `age` if you want an exact calibration for that species-tree node.
Use `age_min` / `age_max` if you want a confidence interval to be propagated to the corresponding gene-tree speciation nodes.
#### `--max_age`
If duplication nodes are deeper than the root node of the species tree, this value will be used as an upper limit of the root node.
#### `--chronos_lambda`
Passed to `chronos` for divergence time estimation. See `chronos` in the [**ape** documentation](https://www.rdocumentation.org/packages/ape/versions/5.2/topics/chronos).
#### `--chronos_model`
Passed to `chronos` for divergence time estimation. Supported values are `discrete`, `relaxed`, and `correlated`.
If an unsupported value is given (e.g., typo like `difscrete`), RADTE now exits with a clear error and suggestion.
See `chronos` in the [**ape** documentation](https://www.rdocumentation.org/packages/ape/versions/5.2/topics/chronos).
#### `--pad_short_edge`
Prohibit dated branches shorter than this value. If detected, the branch length is readjusted by transferring a small portion of branch length from the parent branch.
#### `--allow_constraint_drop`
Default: `true`. After RS retries fail, chronos may try S-only and then R-only constraints. Original ranges are never changed. With `false`, RADTE retains every RS constraint and may use a deterministic feasible tree if numerical optimization fails; the manifest identifies that fallback. Infeasible input is always an error. MCMCTree ignores this option.

#### `--chronos_attempt_timeout_sec`
Per-attempt timeout (seconds) for each `chronos` call.
Use a non-negative number, or `inf`/`none`/`off` to disable per-attempt timeout.
Default is `60` seconds and applies to both the fast and high-cost profiles.
Every attempt is also bounded by the remaining `--chronos_total_timeout_sec` budget.
#### `--chronos_total_timeout_sec`
Total timeout budget (seconds) across all `chronos` retries (RS + retry strategies + S/R if enabled).
Use a non-negative number, or `inf`/`none`/`off` to disable total budgeting.
Default is `300` seconds. When `--allow_constraint_drop=false`, exhausting the budget proceeds to the deterministic no-drop fallback.
#### `--chronos_iter_max`, `--chronos_eval_max`, `--chronos_dual_iter_max`
Controls the initial fast `ape::chronos` optimization profile. Defaults are `250`, `250`, and `20`.
RADTE tries all calibration/model/seed strategies with this bounded profile first, avoiding long waits on numerically unproductive attempts.
#### `--chronos_high_control_fallback`
Boolean; default is `true`.
If every fast-profile strategy fails, RADTE retries with the historical high-cost limits (`iter.max=100000`, `eval.max=100000`, `dual.iter.max=200`) while respecting the same total timeout budget.
Set this to `false` to disable the high-cost fallback.
#### `--dating_backend`
Dating engine. Supported values are `chronos` (default) and `mcmctree`.
`chronos` uses the current `ape::chronos` workflow and supports `--allow_constraint_drop` plus the `--chronos_*` options.
When speciation-node intervals are supplied through `--species_node_bounds_tsv`, `chronos` uses those `age.min` / `age.max` bounds but cannot force separate internal nodes to share one estimated age unless the bounds are exact.
`mcmctree` runs the external **PAML** program **MCMCTree** on the reconciled gene tree using the transferred root/speciation calibrations.
For repeated speciation events caused by ancestral duplications, RADTE uses **MCMCTree** mirror labels (`#1`, `#2`, ...) so that corresponding gene-tree speciation nodes can share the same age parameter.
Both backends preserve the original calibration table and reject temporally infeasible input. MCMCTree also rejects incomplete or inconsistent mirror groups. Its emitted calibrations are PAML soft priors (exact ages use narrow intervals); see [run contracts](run-contracts.md).
After MCMCTree finishes, RADTE independently recalculates every mirrored node age from the posterior mean time tree and stops if members of a shared group differ beyond numerical tolerance.
The current RADTE integration supports `usedata=1` only and requires a MCMCTree-compatible alignment file.
**BEAST is not yet integrated.**
#### `--mcmctree_seqfile`
Required when `--dating_backend=mcmctree`.
This should be an alignment file readable by **MCMCTree** whose taxon names exactly match the gene-tree tip labels.
#### `--mcmctree_bin`
Optional path to the `mcmctree` executable. Default is `mcmctree` in `PATH`.
#### `--mcmctree_workdir`
Optional staging directory for the **MCMCTree** run. RADTE writes the generated tree/control files there and captures `out.txt`, `mcmc.txt`, `FigTree.tre`, and the stdout/stderr logs. By default each invocation receives a unique `<prefix>_mcmctree_run_<timestamp>_<pid>` directory under `--outdir`, so concurrent jobs do not share scratch files.
#### `--mcmctree_timeout_sec`
Wall-time limit for the external MCMCTree process. The default is unlimited. A timed-out process fails with an explicit diagnostic while preserving its isolated work directory for inspection.
#### `--mcmctree_usedata`
Passed to **MCMCTree** as `usedata`.
The current RADTE integration supports only `1`.
#### `--mcmctree_seqtype`
Passed to **MCMCTree** as `seqtype`. Default is `0`.
#### `--mcmctree_clock`
Passed to **MCMCTree** as `clock`. Default is `2`.
#### `--mcmctree_model`
Passed to **MCMCTree** as `model`. Default is `0`.
#### `--mcmctree_burnin`, `--mcmctree_sampfreq`, `--mcmctree_nsample`, `--mcmctree_ncatG`
Passed to **MCMCTree** as `burnin`, `sampfreq`, `nsample`, and `ncatG`.

## Example 1: RADTE after NOTUNG
For input data, see `data/example_notung_01`.
```
# Run NOTUNG in the reconciliation mode
# Don't forget to specify --parsable
java -jar -Xmx2g Notung-2.9.jar \
-s species_tree.nwk \
-g gene_tree.nwk \
--reconcile \
--infertransfers "false" \
--treeoutput newick \
--parsable \
--speciestag prefix \
--maxtrees 1 \
--nolosses \
--outputdir .

# Run RADTE
./radte.r \
--species_tree=species_tree.nwk \
--gene_tree=gene_tree.nwk.reconciled \
--notung_parsable=gene_tree.nwk.reconciled.parsable.txt \
--max_age=1000 \
--chronos_lambda=1 \
--chronos_model=discrete \
--pad_short_edge=0.001
```
#### species_tree.nwk
![](../img/notung_radte_species_tree.svg)

#### gene_tree.nwk.reconciled
![](../img/notung_radte_gene_tree_input.svg)

#### radte_gene_tree_output.nwk
![](../img/notung_radte_gene_tree_output.svg)

## Example 2: RADTE after GeneRax
For input data, see `data/example_generax_01`.
For your own data, please run GeneRax and obtain a `nhx` file for the gene tree.
In Generax, `--rec-model UndatedDTL` may not be compatible with RADTE, so please use `--rec-model UndatedDL`.
```
./radte.r \
--species_tree=species_tree.nwk \
--generax_nhx=gene_tree.nhx \
--max_age=1000 \
--chronos_lambda=1 \
--chronos_model=discrete \
--pad_short_edge=0.001
```

#### species_tree.nwk
![](../img/generax_radte_species_tree.svg)

#### gene_tree.nhx
![](../img/generax_radte_gene_tree_input.svg)

#### radte_gene_tree_output.nwk
![](../img/generax_radte_gene_tree_output.svg)

## Example 3: transfer species-tree node age CI to gene-tree speciation nodes
The following file is included as `data/example_generax_01/species_node_bounds.tsv`.
```
cat > species_node_bounds.tsv <<'EOF'
node	age_min	age_max
s2	8	12
s1	110	125
EOF
```

Run RADTE with the default `chronos` backend.
```
./radte.r \
--species_tree=species_tree.nwk \
--generax_nhx=gene_tree.nhx \
--species_node_bounds_tsv=species_node_bounds.tsv \
--max_age=1000 \
--chronos_lambda=1 \
--chronos_model=discrete \
--pad_short_edge=0.001
```

This transfers the interval for `s2` and `s1` to the corresponding gene-tree speciation nodes.
With `chronos`, repeated speciation nodes created by ancestral duplications receive the same interval bounds, but they are not forced to have exactly identical estimated ages unless the bounds are exact.

Check the transferred constraints in:
* `radte_species_tree.tsv`: branch-length point ages and the effective `age_min` / `age_max`
* `radte_gene_tree.tsv`: transferred `lower_age` / `upper_age`, `constraint_sp_node`, and `shared_speciation_group`

## Example 4: RADTE with the MCMCTree backend
`MCMCTree` requires an alignment file whose taxon labels match the gene-tree tips.
If you use PHYLIP sequential format, separate each taxon name from the sequence by at least two spaces because **MCMCTree** is strict about this.
```
./radte.r \
--species_tree=species_tree.nwk \
--gene_tree=gene_tree.nwk.reconciled \
--notung_parsable=gene_tree.nwk.reconciled.parsable.txt \
--max_age=1000 \
--dating_backend=mcmctree \
--mcmctree_seqfile=gene_alignment.phy \
--mcmctree_clock=2 \
--mcmctree_model=0 \
--mcmctree_burnin=2000 \
--mcmctree_sampfreq=10 \
--mcmctree_nsample=20000
```

To combine `MCMCTree` with species-node CI:
```
./radte.r \
--species_tree=species_tree.nwk \
--generax_nhx=gene_tree.nhx \
--species_node_bounds_tsv=species_node_bounds.tsv \
--max_age=1000 \
--dating_backend=mcmctree \
--mcmctree_seqfile=gene_alignment.phy
```

When the same labeled species-tree node is mapped to multiple gene-tree speciation nodes, RADTE writes **MCMCTree** mirror labels (`#1`, `#2`, ...) so that those nodes share one age parameter.
Only one member (the `driver`) receives the calibration prior in the generated MCMCTree tree; the remaining `mirror` members reuse it through their common mirror label.

![How RADTE mirrors repeated speciation nodes in MCMCTree](img/mcmctree_speciation_mirror.svg)

The ancestral duplication in this example produces two gene-tree nodes for the same `sp_ab` species-tree event. RADTE emits the calibration prior once at the `driver` node and assigns the same mirror label to both nodes, making their estimated ages identical.

## Output files
See `data/example_notung_01/expected` and `data/example_generax_01/expected` for illustrative reference output. Tests generate fresh artifacts in temporary directories instead of modifying these fixtures.

Every filename below starts with the selected `--prefix`; `radte` is shown as the default.

#### radte_gene_tree_output.nwk
This is the main output file of RADTE. Branch lengths represent the estimated evolutionary time.
Node ages represent the estimated divergence time.
The unit of the branch length is the same as that in the input species tree.

#### radte_*.pdf
RADTE generates pdf files for input and output trees in which nodes are colored (see above examples). Red and blue respectively indicate unconstrained and constrained nodes.
While the divergence time of blue nodes is transferred from the species tree, that of red nodes is estimated.
When the root node is blue, it means the divergence time is either transferred from the species tree or bounded by `--max_age`.

#### radte_calibration_all.tsv
This table contains all identified calibration nodes where the divergence time may be transferred from the species tree to the gene tree.

#### radte_calibration_used.tsv
This table is a subset of `radte_calibration_all.tsv` and contain only calibration nodes that are used to transfer the divergence time.
With the MCMCTree backend, `shared_speciation_group` identifies cross-braced nodes, `mirror_role` distinguishes the `driver` from its `mirror` nodes, and `prior_emitted` records whether a calibration prior was written at that node in the generated MCMCTree tree.
With the chronos backend, RADTE preserves the original hard bounds and validates estimated ages against every retained constraint. Temporally infeasible input is rejected before inference.
If `--allow_constraint_drop=true` (default), a part of calibration points may still be dropped only when all no-drop attempts fail.

Both calibration tables also record `estimated_age`, `within_original_bounds`, `constraint_status`, and `bound_policy`; see [run contracts](run-contracts.md). Node numbers refer to the output gene tree as well as the input gene tree.

#### radte_shared_speciation_ages.tsv
This file is written by the MCMCTree backend. It reports each shared species-tree event, its gene-tree member nodes, the posterior mean/minimum/maximum ages, the maximum age difference among members, and the numerical tolerance used for the final mirror consistency check.

#### radte_gene_tree.tsv
This table summarizes gene tree nodes.
In the column `event`, `S` and `D` respectively denote `speciation node` or `duplication node` inferred by Notung or GeneRax.
The root node is indicated as `S(R)` or `D(R)`.
`lower_sp_node` and `upper_sp_node` together indicate which node/branch of the species tree the gene tree node is mapped.
`constraint_sp_node` identifies speciation constraints that correspond to a single labeled species-tree node, and `shared_speciation_group` marks repeated speciation nodes that are linked to the same species-tree event.

#### radte_species_tree.tsv
This table summarizes species tree nodes.
The table records `age_min` and `age_max` alongside the branch-length point age. Without `--species_node_bounds_tsv`, both bounds equal the point age.

#### radte_calibrated_nodes.txt
This file records what types of gene tree nodes are constrained in the divergence time estimation.
RADTE first attempts to constrain all available calibration points transferred from the species tree (**R**, root node; **S**, speciation node) for the divergence time estimation by `chronos` from the **ape** package.
If the estimation succeeded, the content of this file should be **RS**, because both **R** and **S** nodes were used.
If you supply `--species_node_bounds_tsv`, RADTE may report **RS** even when the input gene tree has no duplication nodes, because `chronos` needs to run to satisfy interval constraints.
If the first estimation failed, RADTE retries while preserving **RS** through edge scaling, multi-seed restarts, and alternative model/lambda parameters. Bounds are never widened or converted to soft bounds.
If all **RS** retries fail and `--allow_constraint_drop=true`, RADTE repeats the same exhaustive retry pipeline at **S**, then at **R** (order: **RS** -> **S** -> **R**).
This differs from the method described in Fukushima and Pollock (2020), where duplication nodes (**D**) may be constrained with the upper and lower limits.

#### radte_mcmctree_*
When `--dating_backend=mcmctree` is used, RADTE also copies the generated **MCMCTree** artifacts into the output directory with the prefix `radte_mcmctree_` (for example `radte_mcmctree_out.txt`, `radte_mcmctree_mcmc.txt`, `radte_mcmctree_FigTree.tre`, and the generated control/tree files).

#### radte_run_manifest.tsv
An audit record written after a successful run. It includes the RADTE, R, `ape`, and `treeio` versions; start/end timestamps; effective options including defaults; dating backend and effective parameters; requested/effective seeds; output location; MCMCTree executable/banner/work directory when applicable; and SHA-256 hashes of input and output files, run ID, completion status, RNG kind, and calibration policy.
