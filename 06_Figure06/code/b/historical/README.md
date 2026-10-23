# Historical read-depth calculation, split by output

The old Figure5c multi-plot script belongs to the read-depth analysis now shown in Figure 6b. It has been split into:

- [compute_metrics.py](compute_metrics.py): original H5AD alignment and pairwise cosine computation only; external paths and historical method coverage retained.
- [plot_reference_heatmap.py](plot_reference_heatmap.py): one methods-by-depth heatmap.
- [plot_pairwise_heatmaps.py](plot_pairwise_heatmaps.py): one faceted pairwise heatmap figure.
- [plot_reference_lines.py](plot_reference_lines.py): one reference-consistency line plot.

Plot scripts accept `--input velocity_cosine_sim.csv --out OUTPUT_DIR` and require seaborn in addition to the demo dependencies. Their original visual encodings and method ordering are preserved; they are historical diagnostic views, not the current released renderer. Use [figure6b.py](../figure6b.py) for the current saved-data plot. The upstream H5AD computation was not executed.
