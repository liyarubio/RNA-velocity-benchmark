# Figure03 — Evaluation on simulated and conventional datasets

Current figure: **Figure03**. Historical source ID: **Figure03**.

![Figure03](<figures/Figure03.png>)

## Evidence and reproduction

740 pairs in the main figure and 698 in Extended Data. Per-dataset maxima are a post hoc envelope.

Data status: `available`. Analysis code: `historical_code_available`.

Portable saved-data checks: [commands and dependencies](<../09_Core/README.md>). Historical scripts retain their original environment and external-input requirements. Layout scripts are provenance, not standalone numerical generators.

## Files

| Role | File | Purpose |
|---|---|---|
| Figures | [Figure03.pdf](<figures/Figure03.pdf>) | Current publication figure or preview |
| Figures | [Figure03.png](<figures/Figure03.png>) | Current publication figure or preview |
| Figures | [Figure03_legend.pdf](<figures/Figure03_legend.pdf>) | Current publication figure or preview |
| Data | [Figure03.xlsx](<data/source_data/Figure03.xlsx>) | Source workbook, saved results or input metadata |
| Data | [Figure3a_high_dim_cosine_similarity.csv](<data/figure3_detail/Figure3a_high_dim_cosine_similarity.csv>) | Source workbook, saved results or input metadata |
| Data | [Figure3bc_low_dim_cosine_similarity.csv](<data/figure3_detail/Figure3bc_low_dim_cosine_similarity.csv>) | Source workbook, saved results or input metadata |
| Data | [Figure3bc_paired_data.csv](<data/figure3_detail/Figure3bc_paired_data.csv>) | Source workbook, saved results or input metadata |
| Data | [Figure3d_pca2_to_30_cosine_similarity.csv](<data/figure3_detail/Figure3d_pca2_to_30_cosine_similarity.csv>) | Source workbook, saved results or input metadata |
| Data | [Figure3d_plot_data.csv](<data/figure3_detail/Figure3d_plot_data.csv>) | Source workbook, saved results or input metadata |
| Data | [Figure3f_aligned_metrics.csv](<data/figure3_detail/Figure3f_aligned_metrics.csv>) | Source workbook, saved results or input metadata |
| Data | [Figure3f_correlation.csv](<data/figure3_detail/Figure3f_correlation.csv>) | Source workbook, saved results or input metadata |
| Data | [Figure3f_direction_metric.csv](<data/figure3_detail/Figure3f_direction_metric.csv>) | Source workbook, saved results or input metadata |
| Data | [Figure3f_global_metric.csv](<data/figure3_detail/Figure3f_global_metric.csv>) | Source workbook, saved results or input metadata |
| Data | [Figure3f_incluster_metric.csv](<data/figure3_detail/Figure3f_incluster_metric.csv>) | Source workbook, saved results or input metadata |
| Data | [Figure3f_transition_matrix_cosine_similarity.csv](<data/figure3_detail/Figure3f_transition_matrix_cosine_similarity.csv>) | Source workbook, saved results or input metadata |
| Data | [Figure3g_hubmap_global_performance.csv](<data/figure3_detail/Figure3g_hubmap_global_performance.csv>) | Source workbook, saved results or input metadata |
| Data | [Figure3g_long_data.csv](<data/figure3_detail/Figure3g_long_data.csv>) | Source workbook, saved results or input metadata |
| Data | [Figure3g_method_medians.csv](<data/figure3_detail/Figure3g_method_medians.csv>) | Source workbook, saved results or input metadata |
| Data | [Figure3h_data10_global.csv](<data/figure3_detail/Figure3h_data10_global.csv>) | Source workbook, saved results or input metadata |
| Data | [Figure3h_dataset10_global_performance.csv](<data/figure3_detail/Figure3h_dataset10_global_performance.csv>) | Source workbook, saved results or input metadata |
| Data | [data_manifest.csv](<data/figure3_detail/data_manifest.csv>) | Source workbook, saved results or input metadata |
| Code | [reproduce.py](<../09_Core/validation/reproduce.py>) | Numerical consistency checks on saved data; no model fitting |
| Code | [plot_saved_panels.py](<../09_Core/visualization/plot_saved_panels.py>) | Diagnostic plots of saved data; not replacements for publication layouts |
| Code | [Figure3g_01a_hubmap_global_eval_cell2fate.py](<code/historical/figure3/Figure3g_01a_hubmap_global_eval_cell2fate.py>) | Historical analysis or plotting implementation; original path, environment and external-input dependencies remain |
| Code | [Figure3g_01b_hubmap_global_eval_topicvelo.py](<code/historical/figure3/Figure3g_01b_hubmap_global_eval_topicvelo.py>) | Historical analysis or plotting implementation; original path, environment and external-input dependencies remain |
| Code | [Figure3h_01_dataset10_compute_global_metrics.py](<code/historical/figure3/Figure3h_01_dataset10_compute_global_metrics.py>) | Historical analysis or plotting implementation; original path, environment and external-input dependencies remain |
| Code | [config.py](<../09_Core/model_execution/Run_velocity/config.py>) | Historical analysis or plotting implementation; original path, environment and external-input dependencies remain |
| Code | [run_bivi.py](<../09_Core/model_execution/Run_velocity/run_bivi.py>) | Historical analysis or plotting implementation; original path, environment and external-input dependencies remain |
| Code | [run_cell2fate.py](<../09_Core/model_execution/Run_velocity/run_cell2fate.py>) | Historical analysis or plotting implementation; original path, environment and external-input dependencies remain |
| Code | [run_celldancer.py](<../09_Core/model_execution/Run_velocity/run_celldancer.py>) | Historical analysis or plotting implementation; original path, environment and external-input dependencies remain |
| Code | [run_deepkinet.py](<../09_Core/model_execution/Run_velocity/run_deepkinet.py>) | Historical analysis or plotting implementation; original path, environment and external-input dependencies remain |
| Code | [run_deepvelo_gcn.py](<../09_Core/model_execution/Run_velocity/run_deepvelo_gcn.py>) | Historical analysis or plotting implementation; original path, environment and external-input dependencies remain |
| Code | [run_deepvelo_vae.py](<../09_Core/model_execution/Run_velocity/run_deepvelo_vae.py>) | Historical analysis or plotting implementation; original path, environment and external-input dependencies remain |
| Code | [run_dynamo.py](<../09_Core/model_execution/Run_velocity/run_dynamo.py>) | Historical analysis or plotting implementation; original path, environment and external-input dependencies remain |
| Code | [run_latentvelo.py](<../09_Core/model_execution/Run_velocity/run_latentvelo.py>) | Historical analysis or plotting implementation; original path, environment and external-input dependencies remain |
| Code | [run_multivelo.py](<../09_Core/model_execution/Run_velocity/run_multivelo.py>) | Historical analysis or plotting implementation; original path, environment and external-input dependencies remain |
| Code | [run_phylovelo.py](<../09_Core/model_execution/Run_velocity/run_phylovelo.py>) | Historical analysis or plotting implementation; original path, environment and external-input dependencies remain |
| Code | [run_scvelo.py](<../09_Core/model_execution/Run_velocity/run_scvelo.py>) | Historical analysis or plotting implementation; original path, environment and external-input dependencies remain |
| Code | [run_sdevelo.py](<../09_Core/model_execution/Run_velocity/run_sdevelo.py>) | Historical analysis or plotting implementation; original path, environment and external-input dependencies remain |
| Code | [run_stt.py](<../09_Core/model_execution/Run_velocity/run_stt.py>) | Historical analysis or plotting implementation; original path, environment and external-input dependencies remain |
| Code | [run_tfvelo.py](<../09_Core/model_execution/Run_velocity/run_tfvelo.py>) | Historical analysis or plotting implementation; original path, environment and external-input dependencies remain |
| Code | [run_topicvelo.py](<../09_Core/model_execution/Run_velocity/run_topicvelo.py>) | Historical analysis or plotting implementation; original path, environment and external-input dependencies remain |
| Code | [run_unitvelo.py](<../09_Core/model_execution/Run_velocity/run_unitvelo.py>) | Historical analysis or plotting implementation; original path, environment and external-input dependencies remain |
| Code | [run_veloae.py](<../09_Core/model_execution/Run_velocity/run_veloae.py>) | Historical analysis or plotting implementation; original path, environment and external-input dependencies remain |
| Code | [run_velocyto.py](<../09_Core/model_execution/Run_velocity/run_velocyto.py>) | Historical analysis or plotting implementation; original path, environment and external-input dependencies remain |
| Code | [run_velovae.py](<../09_Core/model_execution/Run_velocity/run_velovae.py>) | Historical analysis or plotting implementation; original path, environment and external-input dependencies remain |
| Code | [run_velovi.py](<../09_Core/model_execution/Run_velocity/run_velovi.py>) | Historical analysis or plotting implementation; original path, environment and external-input dependencies remain |
| Code | [Run_velogif.py](<../09_Core/model_execution/Run_velogif.py>) | Historical analysis or plotting implementation; original path, environment and external-input dependencies remain |
| Code | [run_velogif_batch.py](<../09_Core/model_execution/run_velogif_batch.py>) | Historical analysis or plotting implementation; original path, environment and external-input dependencies remain |

## Additional retained files

- [code/legacy/HuBMAP/1_Annotate_and_process/1_Annotate.R](<code/legacy/HuBMAP/1_Annotate_and_process/1_Annotate.R>)
- [code/legacy/HuBMAP/1_Annotate_and_process/2_Process.ipynb](<code/legacy/HuBMAP/1_Annotate_and_process/2_Process.ipynb>)
- [code/legacy/HuBMAP/2_Run_RNA_velocity/Run_DeepVeloVAE.py](<code/legacy/HuBMAP/2_Run_RNA_velocity/Run_DeepVeloVAE.py>)
- [code/legacy/HuBMAP/2_Run_RNA_velocity/Run_DeepVelo_GCN.py](<code/legacy/HuBMAP/2_Run_RNA_velocity/Run_DeepVelo_GCN.py>)
- [code/legacy/HuBMAP/2_Run_RNA_velocity/Run_LatentVelo.py](<code/legacy/HuBMAP/2_Run_RNA_velocity/Run_LatentVelo.py>)
- [code/legacy/HuBMAP/2_Run_RNA_velocity/Run_STT.py](<code/legacy/HuBMAP/2_Run_RNA_velocity/Run_STT.py>)
- [code/legacy/HuBMAP/2_Run_RNA_velocity/Run_TFvelo.py](<code/legacy/HuBMAP/2_Run_RNA_velocity/Run_TFvelo.py>)
- [code/legacy/HuBMAP/2_Run_RNA_velocity/Run_UniTVelo.py](<code/legacy/HuBMAP/2_Run_RNA_velocity/Run_UniTVelo.py>)
- [code/legacy/HuBMAP/2_Run_RNA_velocity/Run_cellDancer.py](<code/legacy/HuBMAP/2_Run_RNA_velocity/Run_cellDancer.py>)
- [code/legacy/HuBMAP/2_Run_RNA_velocity/Run_scVelo_dynamic.py](<code/legacy/HuBMAP/2_Run_RNA_velocity/Run_scVelo_dynamic.py>)
- [code/legacy/HuBMAP/2_Run_RNA_velocity/Run_veloAE.py](<code/legacy/HuBMAP/2_Run_RNA_velocity/Run_veloAE.py>)
- [code/legacy/HuBMAP/2_Run_RNA_velocity/Run_veloVAE.ipynb](<code/legacy/HuBMAP/2_Run_RNA_velocity/Run_veloVAE.ipynb>)
- [code/legacy/HuBMAP/2_Run_RNA_velocity/Run_veloVAE.py](<code/legacy/HuBMAP/2_Run_RNA_velocity/Run_veloVAE.py>)
- [code/legacy/HuBMAP/2_Run_RNA_velocity/Run_veloVI.py](<code/legacy/HuBMAP/2_Run_RNA_velocity/Run_veloVI.py>)
- [code/legacy/HuBMAP/2_Run_RNA_velocity/Run_velocyto.py](<code/legacy/HuBMAP/2_Run_RNA_velocity/Run_velocyto.py>)
- [code/legacy/HuBMAP/2_Run_RNA_velocity/run_Dynamo.py](<code/legacy/HuBMAP/2_Run_RNA_velocity/run_Dynamo.py>)
- [code/legacy/HuBMAP/2_Run_RNA_velocity/run_Dynamo_stochastic.py](<code/legacy/HuBMAP/2_Run_RNA_velocity/run_Dynamo_stochastic.py>)
- [code/legacy/HuBMAP/3_eval.py](<code/legacy/HuBMAP/3_eval.py>)
- [code/legacy/Real_Datasets10/1_data_process.py](<code/legacy/Real_Datasets10/1_data_process.py>)
- [code/legacy/Real_Datasets10/2_Run_velocity/run_DeepVelo_GCN.py](<code/legacy/Real_Datasets10/2_Run_velocity/run_DeepVelo_GCN.py>)
- [code/legacy/Real_Datasets10/2_Run_velocity/run_DeepVelo_VAE.py](<code/legacy/Real_Datasets10/2_Run_velocity/run_DeepVelo_VAE.py>)
- [code/legacy/Real_Datasets10/2_Run_velocity/run_Dynamo.py](<code/legacy/Real_Datasets10/2_Run_velocity/run_Dynamo.py>)
- [code/legacy/Real_Datasets10/2_Run_velocity/run_LatentVelo.py](<code/legacy/Real_Datasets10/2_Run_velocity/run_LatentVelo.py>)
- [code/legacy/Real_Datasets10/2_Run_velocity/run_STT.py](<code/legacy/Real_Datasets10/2_Run_velocity/run_STT.py>)
- [code/legacy/Real_Datasets10/2_Run_velocity/run_TFvelo.py](<code/legacy/Real_Datasets10/2_Run_velocity/run_TFvelo.py>)
- [code/legacy/Real_Datasets10/2_Run_velocity/run_UnitVelo.py](<code/legacy/Real_Datasets10/2_Run_velocity/run_UnitVelo.py>)
- [code/legacy/Real_Datasets10/2_Run_velocity/run_cellDancer.py](<code/legacy/Real_Datasets10/2_Run_velocity/run_cellDancer.py>)
- [code/legacy/Real_Datasets10/2_Run_velocity/run_scVelo.py](<code/legacy/Real_Datasets10/2_Run_velocity/run_scVelo.py>)
- [code/legacy/Real_Datasets10/2_Run_velocity/run_veloAE.py](<code/legacy/Real_Datasets10/2_Run_velocity/run_veloAE.py>)
- [code/legacy/Real_Datasets10/2_Run_velocity/run_veloVAE.py](<code/legacy/Real_Datasets10/2_Run_velocity/run_veloVAE.py>)
- [code/legacy/Real_Datasets10/2_Run_velocity/run_veloVI.py](<code/legacy/Real_Datasets10/2_Run_velocity/run_veloVI.py>)
- [code/legacy/Real_Datasets10/2_Run_velocity/run_velocyto.py](<code/legacy/Real_Datasets10/2_Run_velocity/run_velocyto.py>)
- [code/legacy/Real_Datasets10/3_eval.py](<code/legacy/Real_Datasets10/3_eval.py>)
- [code/legacy/Visualization/Vis_summary/1_Summary_dyngen_methods.R](<code/legacy/Visualization/Vis_summary/1_Summary_dyngen_methods.R>)
- [code/legacy/Visualization/Vis_summary/2_scatter_and_distribution_dyngen.R](<code/legacy/Visualization/Vis_summary/2_scatter_and_distribution_dyngen.R>)
- [code/legacy/Visualization/Vis_summary/4_HuBMAP.R](<code/legacy/Visualization/Vis_summary/4_HuBMAP.R>)
- [code/legacy/Visualization/Vis_summary/5_summary_data10.R](<code/legacy/Visualization/Vis_summary/5_summary_data10.R>)
- [code/legacy/Visualization/Vis_velocity_stream/Vis_datasets10.ipynb](<code/legacy/Visualization/Vis_velocity_stream/Vis_datasets10.ipynb>)

## Related figures

[ED01](<../08_Extended_Data/ExtendedDataFig01/README.md>) · [ED02](<../08_Extended_Data/ExtendedDataFig02/README.md>) · [ED03](<../08_Extended_Data/ExtendedDataFig03/README.md>) · [ED04](<../08_Extended_Data/ExtendedDataFig04/README.md>) · [ED05](<../08_Extended_Data/ExtendedDataFig05/README.md>) · [ED06](<../08_Extended_Data/ExtendedDataFig06/README.md>)

[Repository overview](<../README.md>) · [Current / historical numbering](<../12_Provenance/figure_index.csv>)
