# Figure 3c: Paired high- versus low-dimensional distributions

[Plotting script](figure3c.R) · [Current reference data](../../data/figure3_detail/Figure3bc_paired_data.csv)

Run from the repository root:

```bash
Rscript --vanilla 03_Figure03/code/c/figure3c.R
```

This entrypoint renders only panel c. Set `FIGURE3_OUT` to override `demo_output/figure3`.

Upstream: [shared simulation metric notebook](../shared/simulation_velocity_metrics.ipynb). Uses high-dimensional and PCA2 velocity cosines. The plot retains dataset/method pairs, including deterministic Dynamo, as in the recovered implementation. It also writes the shared pairing table `Figure3bc_paired_data.csv`.

Historical metric generators require external H5ADs and original method environments. They were not rerun for this release. See the [code guide](../README.md) for dependencies and validation limits.
