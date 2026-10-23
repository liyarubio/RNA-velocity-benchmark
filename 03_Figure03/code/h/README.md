# Figure 3h: Ten conventional datasets

[Plotting script](figure3h.R) · [Current reference data](../../data/figure3_detail/Figure3h_dataset10_global_performance.csv)

Run from the repository root:

```bash
Rscript --vanilla 03_Figure03/code/h/figure3h.R
```

This entrypoint renders only panel h. Set `FIGURE3_OUT` to override `demo_output/figure3`.

- [Historical metric generation: Figure3h_01_dataset10_compute_global_metrics.py](Figure3h_01_dataset10_compute_global_metrics.py)

Historical metric generators require external H5ADs and original method environments. They were not rerun for this release. See the [code guide](../README.md) for dependencies and validation limits.

## Formal panel: descending mean order

Run `python 03_Figure03/code/h/figure3h_formal.py --out demo_output/figure3h_mean` from the repository root. Install [the shared rendering dependencies](../../../09_Core/visualization/mean_top3/requirements.txt).

Methods run from highest to lowest mean GDC across their available values in the ten datasets. Missing values remain missing, and ties retain the previous display order. The per-dataset **Best performance** envelope remains at the far right and is excluded from the ranking. The command writes a PDF/SVG/PNG panel, `method_order.csv` with the mean and valid dataset count, and validation JSON. The R alternative applies the same ordering. No model evaluation is rerun.

The **Mean** row below the heatmap shows the 19 method means as colored cells, using the same RdYlBu_r palette and fixed GDC range [-1, 1] as the heatmap, without numeric annotations. It excludes missing dataset values. The best-performance envelope has no method-mean cell. See the [Mean-row validation](../../../12_Provenance/figure3h_mean_color_validation.json).
