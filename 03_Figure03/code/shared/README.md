# Shared Figure 3 computation

[simulation_velocity_metrics.ipynb](simulation_velocity_metrics.ipynb) is the recovered upstream notebook formerly named `Figure3abcdf_cell_level_metric_generation.ipynb`. Its code cells are unchanged. It computes reusable metrics from the same simulation H5ADs; panel-specific plotting is in the separate a-h directories.

Cell numbers below are one-based positions in the notebook, not execution counts.

| Cells | Computation | Consumers |
|---|---|---|
| 1-13 | Imports, method lists and method-specific preparation | a, b, c, d, f |
| 14-15 | High-dimensional velocity metrics | a, b, c, d, f |
| 16-17 | PCA30 velocity metrics | Auxiliary endpoint |
| 18-19 | PCA2 velocity metrics | b, c, f |
| 20-21 | Transition-matrix metrics | f |
| 22-23 | PCA2-30 dimensional sweep | d |

The notebook requires external inputs and original environments; it is an upstream reference, not an independently executed reproduction of the current release. Saved CSVs define the released values. [figure3_common.R](figure3_common.R) provides shared R readers, method labels, colors and export helpers.

[Panel code guide](../README.md)
