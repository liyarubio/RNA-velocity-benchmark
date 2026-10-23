# Figure 3g: HuBMAP performance

[Plotting script](figure3g.R) · [Current reference data](../../data/figure3_detail/Figure3g_hubmap_global_performance.csv)

Run from the repository root:

```bash
Rscript --vanilla 03_Figure03/code/g/figure3g.R
```

This entrypoint renders only panel g. Set `FIGURE3_OUT` to override `demo_output/figure3`.

- [Historical metric generation: Figure3g_01a_hubmap_global_eval_cell2fate.py](Figure3g_01a_hubmap_global_eval_cell2fate.py)
- [Historical metric generation: Figure3g_01b_hubmap_global_eval_topicvelo.py](Figure3g_01b_hubmap_global_eval_topicvelo.py)

Historical metric generators require external H5ADs and original method environments. They were not rerun for this release. See the [code guide](../README.md) for dependencies and validation limits.
