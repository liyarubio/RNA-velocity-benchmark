# Figure 3b: Method-level high- versus low-dimensional comparison

[Plotting script](figure3b.R) · [Current reference data](../../data/figure3_detail/Figure3bc_paired_data.csv)

Run from the repository root:

```bash
Rscript --vanilla 03_Figure03/code/b/figure3b.R
```

This entrypoint renders only panel b. Set `FIGURE3_OUT` to override `demo_output/figure3`.

Upstream: [panel-specific metric notebook](generate_metrics.ipynb). Uses high-dimensional and PCA2 velocity cosines. The plot summarizes each method; deterministic Dynamo and VeloAE are excluded as in the recovered implementation.

Historical metric generators require external H5ADs and original method environments. They were not rerun for this release. See the [code guide](../README.md) for dependencies and validation limits.
