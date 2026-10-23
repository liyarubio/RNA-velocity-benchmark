# Figure 3d: PCA dimension sensitivity (2-30)

[Plotting script](figure3d.R) · [Current reference data](../../data/figure3_detail/Figure3d_plot_data.csv)

Run from the repository root:

```bash
Rscript --vanilla 03_Figure03/code/d/figure3d.R
```

This entrypoint renders only panel d. Set `FIGURE3_OUT` to override `demo_output/figure3`.

Upstream: [shared simulation metric notebook](../shared/simulation_velocity_metrics.ipynb). Uses the PCA2-30 dimension sweep and high-dimensional velocity cosines.

Historical metric generators require external H5ADs and original method environments. They were not rerun for this release. See the [code guide](../README.md) for dependencies and validation limits.
