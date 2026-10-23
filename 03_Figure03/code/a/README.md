# Figure 3a: High-dimensional simulation cosines

[Plotting script](figure3a.R) · [Current reference data](../../data/figure3_detail/Figure3a_high_dim_cosine_similarity.csv)

Run from the repository root:

```bash
Rscript --vanilla 03_Figure03/code/a/figure3a.R
```

This entrypoint renders only panel a. Set `FIGURE3_OUT` to override `demo_output/figure3`.

Upstream: [panel-specific metric notebook](generate_metrics.ipynb). Uses the high-dimensional velocity cosine block.

Historical metric generators require external H5ADs and original method environments. They were not rerun for this release. See the [code guide](../README.md) for dependencies and validation limits.

[Optional saved-data heatmap diagnostic](figure3a_saved_data.py) is the single-panel extraction from the former Python batch script. The R file above remains the detailed renderer.

## Current formal presentation

Run `python 03_Figure03/code/a/figure3a_formal.py --out demo_output/mean_top3` from the repository root. [Formal entrypoint](figure3a_formal.py) removes the Median display column and aligns Mean ranks 1–3. Install [renderer dependencies](../../../09_Core/visualization/mean_top3/requirements.txt). Raw source values and the original normalization are retained.
