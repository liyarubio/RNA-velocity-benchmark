# Figure 3 code, by panel

The earlier `historical/figure3` folder contained only supplementary g/h metric scripts. It was not a complete panel-code index. The saved-data renderers below were recovered from `Manuscript_2607/Figure3_sources`; all nine unique numerical input CSVs are byte-identical to the current Manuscript_1021 inputs. No older data or formal figures were imported.

| Panel | Content | Plotting entrypoint | Current data / derived reference |
|---|---|---|---|
| a | High-dimensional simulation cosines | [figure3a.R](panels/figure3a.R) | [Figure3a_high_dim_cosine_similarity.csv](../data/figure3_detail/Figure3a_high_dim_cosine_similarity.csv) |
| b, c | High- versus low-dimensional comparisons | [figure3bc.R](panels/figure3bc.R) | [Figure3bc_paired_data.csv](../data/figure3_detail/Figure3bc_paired_data.csv) |
| d | PCA dimension sensitivity (2-30) | [figure3d.R](panels/figure3d.R) | [Figure3d_plot_data.csv](../data/figure3_detail/Figure3d_plot_data.csv) |
| e | Metric-definition schematic | Manual artwork; no numerical analysis code | [Retained formal figure](../figures/Figure03.pdf) |
| f | Agreement between evaluation metrics | [figure3f.R](panels/figure3f.R) | [Figure3f_correlation.csv](../data/figure3_detail/Figure3f_correlation.csv) |
| g | HuBMAP performance | [figure3g.R](panels/figure3g.R) | [Figure3g_hubmap_global_performance.csv](../data/figure3_detail/Figure3g_hubmap_global_performance.csv) |
| h | Ten conventional datasets | [figure3h.R](panels/figure3h.R) | [Figure3h_dataset10_global_performance.csv](../data/figure3_detail/Figure3h_dataset10_global_performance.csv) |

## Saved-data redraws

Install the R packages listed in [requirements-R.txt](panels/requirements-R.txt) in your R environment, then run from the repository root:

```bash
FIGURE3_OUT=/tmp/figure3-panels Rscript --vanilla 03_Figure03/code/panels/render_all.R
```

Individual scripts can be run the same way. Without `FIGURE3_OUT`, outputs go to `demo_output/figure3`. Panel e remains in the retained formal figure; `render_all.R` covers the seven numerical panels a-d and f-h. Redraws are separate from the published layout, which remains unchanged.

The shared helpers are [figure3_common.R](panels/figure3_common.R) and [figure3a_plot_core.R](panels/figure3a_plot_core.R). The latter contains the complete bubble-matrix plotting implementation, rather than a placeholder.

Validation distinguishes source recovery from execution: R syntax and current-data identities were checked, and d/f numerical derivations were recalculated independently. Full R plotting execution was not completed because the available R environments lack required packages. The recovery report records this limitation; these checks do not certify identical formal rendering or a model rerun.

## Independently check saved panel data

```bash
python 03_Figure03/code/validate_panel_data.py --out /tmp/figure3-panel-validation.json
```

This uses the existing Python demo dependencies and recalculates all 435 panel-d dimension summaries, 798 panel-f aligned units and their pairwise-complete Pearson correlations. Missing values remain missing.

## Upstream metric generation

- [Cell-level velocity and transition metrics for a-d/f](historical/figure3/Figure3abcdf_cell_level_metric_generation.ipynb)
- [Global, in-cluster and directional metrics for f](historical/figure3/Figure3f_global_incluster_direction_generation.ipynb)
- [HuBMAP cell2fate metrics for g](historical/figure3/Figure3g_01a_hubmap_global_eval_cell2fate.py)
- [HuBMAP TopicVelo metrics for g](historical/figure3/Figure3g_01b_hubmap_global_eval_topicvelo.py)
- [Dataset10 metrics for h](historical/figure3/Figure3h_01_dataset10_compute_global_metrics.py)
- [Shared dyngen simulation, model execution and evaluation](../../09_Core/simulation/dyngen/README.md)
- [Earlier a/b/c visualization implementations](legacy/Visualization/Vis_summary/)

The recovered notebooks preserve every code cell; historical outputs and execution counts were cleared. They require external H5ADs, original paths and method environments, and have not been executed in this release. Their method sets may differ from current panel coverage; the current CSVs and workbook define the released values.

[Recovery hashes and scope](../../12_Provenance/figure3_code_recovery.json) · [Figure overview](../README.md)
