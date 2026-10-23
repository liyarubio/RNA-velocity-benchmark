# Figure 3f: Agreement between evaluation metrics

[Plotting script](figure3f.R) · [Current reference data](../../data/figure3_detail/Figure3f_correlation.csv)

Run from the repository root:

```bash
Rscript --vanilla 03_Figure03/code/f/figure3f.R
```

This entrypoint renders only panel f. Set `FIGURE3_OUT` to override `demo_output/figure3`.

Upstream: [shared simulation metric notebook](../shared/simulation_velocity_metrics.ipynb). Uses velocity and transition-matrix cosines together with [global, in-cluster and directional metrics](global_incluster_direction_metrics.ipynb).

Historical metric generators require external H5ADs and original method environments. They were not rerun for this release. See the [code guide](../README.md) for dependencies and validation limits.
