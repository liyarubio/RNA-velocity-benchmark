# Figure03 — Evaluation on simulated and conventional datasets

Current release: **Manuscript_1021**, including figure presentation updates through **23 October 2026**. Historical numerical source ID: **Figure03**.

740 pairs in the main figure and 698 in Extended Data. Per-dataset maxima are a post hoc envelope.

![Figure03](<figures/Figure03.png>)

## Panel code a-h

[Open the complete panel-by-panel code guide](code/README.md). Each panel has its own directory under `code/a` through `code/h`; each numerical panel has its own drawing entrypoint, and a/b/c/d/f have separate metric notebooks; e is retained manual artwork. The earlier g/h-only historical folder was a subset of the available code. Recovery uses the current 1021 CSVs, with hashes and execution limits recorded in the guide.

## Mean-only summary update

The simulation summary columns now display Mean only, with ranks 1–3 aligned inside the bars. Numerical inputs and normalization are unchanged. Companion caption text is retained for author revision. See the panel code guide for independent formal rendering commands.

## Panel h method order

Methods are ordered left to right by descending available-case mean GDC across the ten classic datasets. The bottom Mean row shows each method mean using the same heatmap colors and fixed GDC scale, without numeric annotations. Missing values are retained; the best-performance envelope stays at the right. See the [independent formal entrypoint](code/h/figure3h_formal.py).

## Current files

| Role | File |
|---|---|
| results | [figures/Figure03.pdf](<figures/Figure03.pdf>) |
| results | [figures/Figure03.png](<figures/Figure03.png>) |
| results | [figures/Figure03.svg](<figures/Figure03.svg>) |
| results | [figures/legends/Figure03_legend.pdf](<figures/legends/Figure03_legend.pdf>) |
| data | [data/source_data/Figure03.xlsx](<data/source_data/Figure03.xlsx>) |
| data | [data/figure3_detail/Figure3a_high_dim_cosine_similarity.csv](<data/figure3_detail/Figure3a_high_dim_cosine_similarity.csv>) |
| data | [data/figure3_detail/Figure3bc_low_dim_cosine_similarity.csv](<data/figure3_detail/Figure3bc_low_dim_cosine_similarity.csv>) |
| data | [data/figure3_detail/Figure3bc_paired_data.csv](<data/figure3_detail/Figure3bc_paired_data.csv>) |
| data | [data/figure3_detail/Figure3d_pca2_to_30_cosine_similarity.csv](<data/figure3_detail/Figure3d_pca2_to_30_cosine_similarity.csv>) |
| data | [data/figure3_detail/Figure3d_plot_data.csv](<data/figure3_detail/Figure3d_plot_data.csv>) |
| data | [data/figure3_detail/Figure3f_aligned_metrics.csv](<data/figure3_detail/Figure3f_aligned_metrics.csv>) |
| data | [data/figure3_detail/Figure3f_correlation.csv](<data/figure3_detail/Figure3f_correlation.csv>) |
| data | [data/figure3_detail/Figure3f_direction_metric.csv](<data/figure3_detail/Figure3f_direction_metric.csv>) |
| data | [data/figure3_detail/Figure3f_global_metric.csv](<data/figure3_detail/Figure3f_global_metric.csv>) |
| data | [data/figure3_detail/Figure3f_incluster_metric.csv](<data/figure3_detail/Figure3f_incluster_metric.csv>) |
| data | [data/figure3_detail/Figure3f_transition_matrix_cosine_similarity.csv](<data/figure3_detail/Figure3f_transition_matrix_cosine_similarity.csv>) |
| data | [data/figure3_detail/Figure3g_hubmap_global_performance.csv](<data/figure3_detail/Figure3g_hubmap_global_performance.csv>) |
| data | [data/figure3_detail/Figure3g_long_data.csv](<data/figure3_detail/Figure3g_long_data.csv>) |
| data | [data/figure3_detail/Figure3g_method_medians.csv](<data/figure3_detail/Figure3g_method_medians.csv>) |
| data | [data/figure3_detail/Figure3h_data10_global.csv](<data/figure3_detail/Figure3h_data10_global.csv>) |
| data | [data/figure3_detail/Figure3h_dataset10_global_performance.csv](<data/figure3_detail/Figure3h_dataset10_global_performance.csv>) |
| data | [data/figure3_detail/data_manifest.csv](<data/figure3_detail/data_manifest.csv>) |
| code | [code/g/Figure3g_01a_hubmap_global_eval_cell2fate.py](<code/g/Figure3g_01a_hubmap_global_eval_cell2fate.py>) |
| code | [code/g/Figure3g_01b_hubmap_global_eval_topicvelo.py](<code/g/Figure3g_01b_hubmap_global_eval_topicvelo.py>) |
| code | [code/h/Figure3h_01_dataset10_compute_global_metrics.py](<code/h/Figure3h_01_dataset10_compute_global_metrics.py>) |
| code | [../09_Core/model_execution/Run_velocity/config.py](<../09_Core/model_execution/Run_velocity/config.py>) |
| code | [../09_Core/model_execution/Run_velocity/run_bivi.py](<../09_Core/model_execution/Run_velocity/run_bivi.py>) |
| code | [../09_Core/model_execution/Run_velocity/run_cell2fate.py](<../09_Core/model_execution/Run_velocity/run_cell2fate.py>) |
| code | [../09_Core/model_execution/Run_velocity/run_celldancer.py](<../09_Core/model_execution/Run_velocity/run_celldancer.py>) |
| code | [../09_Core/model_execution/Run_velocity/run_deepkinet.py](<../09_Core/model_execution/Run_velocity/run_deepkinet.py>) |
| code | [../09_Core/model_execution/Run_velocity/run_deepvelo_gcn.py](<../09_Core/model_execution/Run_velocity/run_deepvelo_gcn.py>) |
| code | [../09_Core/model_execution/Run_velocity/run_deepvelo_vae.py](<../09_Core/model_execution/Run_velocity/run_deepvelo_vae.py>) |
| code | [../09_Core/model_execution/Run_velocity/run_dynamo.py](<../09_Core/model_execution/Run_velocity/run_dynamo.py>) |
| code | [../09_Core/model_execution/Run_velocity/run_latentvelo.py](<../09_Core/model_execution/Run_velocity/run_latentvelo.py>) |
| code | [../09_Core/model_execution/Run_velocity/run_multivelo.py](<../09_Core/model_execution/Run_velocity/run_multivelo.py>) |
| code | [../09_Core/model_execution/Run_velocity/run_phylovelo.py](<../09_Core/model_execution/Run_velocity/run_phylovelo.py>) |
| code | [../09_Core/model_execution/Run_velocity/run_scvelo.py](<../09_Core/model_execution/Run_velocity/run_scvelo.py>) |
| code | [../09_Core/model_execution/Run_velocity/run_sdevelo.py](<../09_Core/model_execution/Run_velocity/run_sdevelo.py>) |
| code | [../09_Core/model_execution/Run_velocity/run_stt.py](<../09_Core/model_execution/Run_velocity/run_stt.py>) |
| code | [../09_Core/model_execution/Run_velocity/run_tfvelo.py](<../09_Core/model_execution/Run_velocity/run_tfvelo.py>) |
| code | [../09_Core/model_execution/Run_velocity/run_topicvelo.py](<../09_Core/model_execution/Run_velocity/run_topicvelo.py>) |
| code | [../09_Core/model_execution/Run_velocity/run_unitvelo.py](<../09_Core/model_execution/Run_velocity/run_unitvelo.py>) |
| code | [../09_Core/model_execution/Run_velocity/run_veloae.py](<../09_Core/model_execution/Run_velocity/run_veloae.py>) |
| code | [../09_Core/model_execution/Run_velocity/run_velocyto.py](<../09_Core/model_execution/Run_velocity/run_velocyto.py>) |
| code | [../09_Core/model_execution/Run_velocity/run_velovae.py](<../09_Core/model_execution/Run_velocity/run_velovae.py>) |
| code | [../09_Core/model_execution/Run_velocity/run_velovi.py](<../09_Core/model_execution/Run_velocity/run_velovi.py>) |
| code | [../09_Core/model_execution/Run_velogif.py](<../09_Core/model_execution/Run_velogif.py>) |
| code | [../09_Core/model_execution/run_velogif_batch.py](<../09_Core/model_execution/run_velogif_batch.py>) |
| code | [code/a/figure3a.R](<code/a/figure3a.R>) |
| code | [code/a/figure3a_plot_core.R](<code/a/figure3a_plot_core.R>) |
| code | [code/b/figure3b.R](<code/b/figure3b.R>) |
| code | [code/c/figure3c.R](<code/c/figure3c.R>) |
| code | [code/d/figure3d.R](<code/d/figure3d.R>) |
| code | [code/f/figure3f.R](<code/f/figure3f.R>) |
| code | [code/g/figure3g.R](<code/g/figure3g.R>) |
| code | [code/h/figure3h.R](<code/h/figure3h.R>) |
| code | [code/shared/figure3_common.R](<code/shared/figure3_common.R>) |
| code | [code/render_all.R](<code/render_all.R>) |
| code | [code/validate_panel_data.py](<code/validate_panel_data.py>) |

## Code and validation

[Portable checks and workflow boundaries](<../09_Core/README.md>). Saved-data checks do not rerun models. Historical scripts require their original environments and external inputs.

[Retained figure-specific code](<code>)

[Repository overview](<../README.md>) · [Release notes](<../12_Provenance/RELEASE_1021.md>)

## Panel-specific metric notebooks

- [a: generate_metrics.ipynb](code/a/generate_metrics.ipynb)
- [b: generate_metrics.ipynb](code/b/generate_metrics.ipynb)
- [c: generate_metrics.ipynb](code/c/generate_metrics.ipynb)
- [d: generate_metrics.ipynb](code/d/generate_metrics.ipynb)
- [f: generate_metrics.ipynb](code/f/generate_metrics.ipynb)
