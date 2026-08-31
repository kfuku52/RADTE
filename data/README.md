# Example data

`example_notung_01` and `example_generax_01` contain immutable input fixtures. Their `expected/` directories contain historical illustrations from before RADTE 0.6.0; the exact generating versions were not recorded. They retain older table schemas and calibration values that must not be used to assess the current bound-preserving algorithm. Tests do not treat these files as numerical/PDF golden snapshots and never write into these directories.

Run the [current examples](../docs/usage.md#example-1-radte-after-notung) to generate current tables, plots, and manifests in separate `results/` directories. See each `expected/README.md` for the historical limitations.

Fresh integration outputs are generated in temporary directories and validated structurally (tree topology, ultrametricity, calibration bounds, tables, and manifests). Historical debug-only tables and PDFs were removed because they were unreferenced and could be confused with supported outputs.
