# Run contracts and migration to 0.6

## Calibration semantics

Chronos retains each selected calibration's original `age.min` and `age.max`. Overlapping ancestor/descendant ranges are valid when a positive-length dated tree exists; the ancestor's lower bound need not exceed the descendant's lower bound. Feasibility uses postorder lower-bound and preorder upper-bound propagation without modifying the input table. Infeasible ranges fail before optimization, even when constraint dropping is enabled.

Retries may change optimizer controls, starting seeds, input edge scaling, model, or lambda. They never widen exact ages, lower bounds, or turn hard constraints into soft bounds. `--allow_constraint_drop=true` permits RS → S → R selection after failures. With `false`, only RS is retained; a deterministic feasible tree may be returned when optimization fails. Its manifest has `chronos_model_used=deterministic-fallback` and no effective random seed. That tree is a constraint-respecting construction, not a maximum-likelihood or posterior estimate.

Chronos results must have finite positive edges, the same topology and tips as the gene input, ultrametric tips, and ages inside retained bounds. Age/ultrametric comparison uses `1e-7 * max(1, tree_height)` for floating-point tolerance. Feasibility requires positive edges; the default numerical margin is `16 * .Machine$double.eps * max(1, maximum_calibration_age)`. `--pad_short_edge` additionally projects short edges within the same hard bounds and fails if that minimum is infeasible. It preserves existing inferred ages where feasible, records `branch_padding_applied`, and no longer changes species-tree calibration ages.

MCMCTree retains its PAML prior semantics: `B{lower, upper}` and root `>lower<upper` are **soft** priors with default tail probabilities, not hard intervals. This follows PAML's [ProcessNodeAnnotation implementation](https://github.com/abacus-gene/paml/blob/master/src/treesub.c). Exact input ages are represented by a narrow prior interval of ±`max(1e-8, abs(age) * 1e-6)` because PAML requires an interval. A posterior mean outside a nominal bound is reported, not silently clipped or rejected as a hard-bound violation. Output structure/topology and shared speciation equality are validated; structural age tolerance is `1e-5 * max(1, tree_height)` to accommodate printed posterior rounding. None of these checks establishes MCMC convergence.

## Input and node identity

- NOTUNG and GeneRax resolve every lower/upper species node against the same dated species-node table. Unknown labels are errors. A missing upper node is allowed only for a duplication whose lower node is the species-tree root; its upper age is `--max_age`.
- A node's outer Newick single quotes are removed as syntax, while IDs such as `001` retain their zeros. TSV labels are read as strings; age columns are converted separately.
- allS requires one gene per represented species and the same topology as the pruned species tree. Node ages are transferred onto the original gene topology and numbering. MCMCTree output is likewise mapped back onto gene node IDs.
- The bundled NOTUNG metadata had aliases from a different species-tree export. Its numeric aliases are now the supplied tree's `n` labels, and `n54`'s obsolete upper alias is the parent `Root`. Historical `expected/` files are illustrative only; regenerate new output separately.

## Calibration tables

Existing filenames and node IDs remain available. `age.min`/`age.max` now always describe the original input range. Both `_calibration_all.tsv` and `_calibration_used.tsv` include:

| Column | Meaning |
| --- | --- |
| `estimated_age` | Age at this node in the output gene tree |
| `within_original_bounds` | Whether the estimated age falls within the original bounds, using the backend's numerical tolerance |
| `constraint_status` | `retained`, `dropped` (an omitted R/S constraint), or `not-selected` (a non-root D interval) |
| `bound_policy` | `hard` for chronos or `PAML-soft-prior` for MCMCTree |
| `soft.bounds` | `FALSE` for chronos; `TRUE` for MCMCTree prior semantics |

Mirror roles and `prior_emitted` keep their existing meanings. A mirror can be retained through a shared driver without emitting a second prior.

## Reproducibility

The same seed and the same R/RNG/ape environment give the same chronos random stream in both finite-timeout Unix children and inline execution. The caller's RNG state is restored. Timeout-driven strategy selection can still vary with machine load; record the effective model, lambda, seed, and timeout settings when comparing runs. No cross-version or cross-platform bitwise guarantee is made.

Manifests record effective options including defaults, `rng_kind`, actual backend parameters, original input hashes, output hashes, `run_id`, and `run_status=complete`. No seed is reported as effective for allS or deterministic construction.

## Output and input safety

All user-facing artifacts are built in a unique staging directory on the output filesystem. A prefix lock prevents concurrent writers. If any writer, including PDF generation, fails, the preceding complete generation remains unchanged. Publication backs up the previous generation, publishes the manifest last, and rolls back on errors. Artifacts from an earlier backend that are absent from the new generation are removed on successful publication.

The flat multi-file layout cannot provide a single atomic snapshot to unlocked readers or survive every filesystem/power failure. Readers requiring consistency should wait for the prefix lock to disappear and verify the `output_sha256.*` entries against the manifest. A killed process can leave a `.PREFIX.radte-lock` directory; inspect its owner file and ensure no run is active before removing it. Incomplete rollback retains its backup directory and reports the recovery path.

Inputs must not alias output artifacts, including through symlinks. MCMCTree checks protected inputs before cleanup, validates and stages all new control/tree/alignment files before replacing old workdir artifacts, normalizes its executable before changing directories, and receives an alignment copy. Default workdirs are unique. Explicit workdirs are diagnostic scratch space and may contain a failed new run; user-facing outputs are published only on success.
