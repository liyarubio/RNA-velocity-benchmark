# dyngen simulation pipeline

The original `Simulation_by_dyngen` directory is retained here with its internal structure and source bytes unchanged.

1. [Download, generate and process datasets](1_simulation/): R generation scripts and the AnnData processing notebook.
2. [Run velocity methods](2_run_velocity/): method-specific implementations.
3. [Evaluate saved outputs](3_eval_metrics/): ground-truth and comparison notebooks.
4. [Dimension comparisons](4_corplot_dim/): R plotting scripts.

Inspect each script's input paths and environment before execution. The method environments and external data are not supplied by `requirements-demo.txt`, which supports saved-data checks only. No simulation or model rerun was performed during this reorganization.

[Simulation overview](../README.md)
