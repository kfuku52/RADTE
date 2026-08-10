# Example data

`example_notung_01` and `example_generax_01` contain immutable input fixtures. Their `expected/` directories contain illustrative reference artifacts from a successful RADTE run; tests do not treat floating-point/PDF bytes as golden snapshots and never write into these directories.

Fresh integration outputs are generated in temporary directories and validated structurally (tree topology, ultrametricity, calibration bounds, tables, and manifests). Historical debug-only tables and PDFs were removed because they were unreferenced and could be confused with supported outputs.
