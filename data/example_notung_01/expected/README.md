# NOTUNG reference output

These are historical, illustrative outputs from before RADTE 0.6.0; they are not numerical golden fixtures. Earlier releases silently replaced unresolved or infeasible ages. They must not be used to assess the new bound-preserving algorithm.

The input reconciliation now uses the species tree's actual labels (`n47`, `n42`, etc.). The obsolete upper alias `n170` for duplication `n54` is replaced by `Root`, the parent of `n47` in the supplied dated tree. Fresh results should be generated with `--outdir` in a separate directory. Tests never overwrite these inputs or historical files.
