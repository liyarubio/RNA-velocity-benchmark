# Figure05 — Spatial evaluation, runtime and read-depth sensitivity

Current figure: **Figure05**. Historical source ID: **Figure05**.

![Figure05](<figures/Figure05.png>)

## Evidence and reproduction

ST_small contains 608 spots. The TopicVelo source value was corrected to 0.4838869706599916, matching the final figure.

Data status: `available`. Analysis code: `historical_code_available`.

Portable saved-data checks: [commands and dependencies](<../09_Core/README.md>). Historical scripts retain their original environment and external-input requirements. Layout scripts are provenance, not standalone numerical generators.

## Files

| Role | File | Purpose |
|---|---|---|
| Figures | [Figure05.pdf](<figures/Figure05.pdf>) | Current publication figure or preview |
| Figures | [Figure05.png](<figures/Figure05.png>) | Current publication figure or preview |
| Figures | [Figure05_legend.pdf](<figures/Figure05_legend.pdf>) | Current publication figure or preview |
| Data | [Figure05.xlsx](<data/source_data/Figure05.xlsx>) | Source workbook, saved results or input metadata |
| Data | [GDC_display_order.csv](<../10_Shared/data/spatial_completion/GDC_display_order.csv>) | Source workbook, saved results or input metadata |
| Data | [ST_raw_status.json](<../10_Shared/data/spatial_completion/ST_raw_status.json>) | Source workbook, saved results or input metadata |
| Data | [ST_small_status.json](<../10_Shared/data/spatial_completion/ST_small_status.json>) | Source workbook, saved results or input metadata |
| Data | [TFvelo_evaluation_summary.csv](<../10_Shared/data/spatial_completion/TFvelo_evaluation_summary.csv>) | Source workbook, saved results or input metadata |
| Data | [completion_audit.json](<../10_Shared/data/spatial_completion/completion_audit.json>) | Source workbook, saved results or input metadata |
| Data | [evaluation_probe.json](<../10_Shared/data/spatial_completion/evaluation_probe.json>) | Source workbook, saved results or input metadata |
| Data | [spatial_cells.csv](<../10_Shared/data/spatial_completion/spatial_cells.csv>) | Source workbook, saved results or input metadata |
| Data | [spatial_directions.csv](<../10_Shared/data/spatial_completion/spatial_directions.csv>) | Source workbook, saved results or input metadata |
| Data | [ST_raw.json](<../10_Shared/data/spatial_completion/topicvelo_retry/ST_raw.json>) | Source workbook, saved results or input metadata |
| Data | [ST_raw_directions.csv](<../10_Shared/data/spatial_completion/topicvelo_retry/ST_raw_directions.csv>) | Source workbook, saved results or input metadata |
| Data | [ST_raw_effective_transition_function.py](<../10_Shared/data/spatial_completion/topicvelo_retry/ST_raw_effective_transition_function.py>) | Source workbook, saved results or input metadata |
| Data | [ST_raw_native_transition.npz](<../10_Shared/data/spatial_completion/topicvelo_retry/ST_raw_native_transition.npz>) | Source workbook, saved results or input metadata |
| Data | [ST_raw_status.json](<../10_Shared/data/spatial_completion/topicvelo_retry/ST_raw_status.json>) | Source workbook, saved results or input metadata |
| Data | [ST_small.json](<../10_Shared/data/spatial_completion/topicvelo_retry/ST_small.json>) | Source workbook, saved results or input metadata |
| Data | [ST_small_directions.csv](<../10_Shared/data/spatial_completion/topicvelo_retry/ST_small_directions.csv>) | Source workbook, saved results or input metadata |
| Data | [ST_small_effective_transition_function.py](<../10_Shared/data/spatial_completion/topicvelo_retry/ST_small_effective_transition_function.py>) | Source workbook, saved results or input metadata |
| Data | [ST_small_native_transition.npz](<../10_Shared/data/spatial_completion/topicvelo_retry/ST_small_native_transition.npz>) | Source workbook, saved results or input metadata |
| Data | [ST_small_status.json](<../10_Shared/data/spatial_completion/topicvelo_retry/ST_small_status.json>) | Source workbook, saved results or input metadata |
| Code | [plot_saved_panels.py](<../09_Core/visualization/plot_saved_panels.py>) | Diagnostic plots of saved data; not replacements for publication layouts |
| Code | [complete_topicvelo.py](<../09_Core/shared/spatial_completion/complete_topicvelo.py>) | Historical analysis or plotting implementation; original path, environment and external-input dependencies remain |
| Code | [evaluate_saved_outputs.py](<../09_Core/shared/spatial_completion/evaluate_saved_outputs.py>) | Historical analysis or plotting implementation; original path, environment and external-input dependencies remain |
| Code | [evaluate_topicvelo.py](<../09_Core/shared/spatial_topicvelo_retry/evaluate_topicvelo.py>) | Historical analysis or plotting implementation; original path, environment and external-input dependencies remain |
| Code | [retry_topicvelo.py](<../09_Core/shared/spatial_topicvelo_retry/retry_topicvelo.py>) | Historical analysis or plotting implementation; original path, environment and external-input dependencies remain |
| Code | [browser_check.py](<../09_Core/shared/stability_completion/browser_check.py>) | Historical analysis or plotting implementation; original path, environment and external-input dependencies remain |
| Code | [check_sparse_assembly.py](<../09_Core/shared/stability_completion/check_sparse_assembly.py>) | Historical analysis or plotting implementation; original path, environment and external-input dependencies remain |
| Code | [complete_topicvelo.py](<../09_Core/shared/stability_completion/complete_topicvelo.py>) | Historical analysis or plotting implementation; original path, environment and external-input dependencies remain |
| Code | [evaluate_stability.py](<../09_Core/shared/stability_completion/evaluate_stability.py>) | Historical analysis or plotting implementation; original path, environment and external-input dependencies remain |
| Code | [prepare_publication.py](<../09_Core/shared/stability_completion/prepare_publication.py>) | Historical analysis or plotting implementation; original path, environment and external-input dependencies remain |
| Code | [render_original_style.py](<../09_Core/shared/stability_completion/render_original_style.py>) | Historical analysis or plotting implementation; original path, environment and external-input dependencies remain |
| Code | [update_html.py](<../09_Core/shared/stability_completion/update_html.py>) | Historical analysis or plotting implementation; original path, environment and external-input dependencies remain |
| Code | [validate_delivery.py](<../09_Core/shared/stability_completion/validate_delivery.py>) | Historical analysis or plotting implementation; original path, environment and external-input dependencies remain |
| Code | [Figure5c_velocity_matrix_similarity_original.py](<code/historical/Figure5_sources/code/panel_c/Figure5c_velocity_matrix_similarity_original.py>) | Historical analysis or plotting implementation; original path, environment and external-input dependencies remain |

## Additional retained files

- [code/legacy/Scalability/1_sampling.py](<code/legacy/Scalability/1_sampling.py>)
- [code/legacy/Scalability/2_Run_velocity/run_Latentvelo.py](<code/legacy/Scalability/2_Run_velocity/run_Latentvelo.py>)
- [code/legacy/Scalability/2_Run_velocity/run_STT.py](<code/legacy/Scalability/2_Run_velocity/run_STT.py>)
- [code/legacy/Scalability/2_Run_velocity/run_UnitVelo.py](<code/legacy/Scalability/2_Run_velocity/run_UnitVelo.py>)
- [code/legacy/Scalability/2_Run_velocity/run_celldancer.py](<code/legacy/Scalability/2_Run_velocity/run_celldancer.py>)
- [code/legacy/Scalability/2_Run_velocity/run_deepvelo_GCN.py](<code/legacy/Scalability/2_Run_velocity/run_deepvelo_GCN.py>)
- [code/legacy/Scalability/2_Run_velocity/run_deepvelo_VAE.py](<code/legacy/Scalability/2_Run_velocity/run_deepvelo_VAE.py>)
- [code/legacy/Scalability/2_Run_velocity/run_dynamo.py](<code/legacy/Scalability/2_Run_velocity/run_dynamo.py>)
- [code/legacy/Scalability/2_Run_velocity/run_scvelo.py](<code/legacy/Scalability/2_Run_velocity/run_scvelo.py>)
- [code/legacy/Scalability/2_Run_velocity/run_veloAE.py](<code/legacy/Scalability/2_Run_velocity/run_veloAE.py>)
- [code/legacy/Scalability/2_Run_velocity/run_veloVAE.py](<code/legacy/Scalability/2_Run_velocity/run_veloVAE.py>)
- [code/legacy/Scalability/2_Run_velocity/run_veloVI.py](<code/legacy/Scalability/2_Run_velocity/run_veloVI.py>)
- [code/legacy/Scalability/2_Run_velocity/run_velocyto.py](<code/legacy/Scalability/2_Run_velocity/run_velocyto.py>)
- [code/legacy/Scalability/3_eval.py](<code/legacy/Scalability/3_eval.py>)
- [code/legacy/Visualization/Vis_summary/7_time.R](<code/legacy/Visualization/Vis_summary/7_time.R>)

## Related figures

[ED09](<../08_Extended_Data/ExtendedDataFig09/README.md>) · [ED10](<../08_Extended_Data/ExtendedDataFig10/README.md>) · [ED11](<../08_Extended_Data/ExtendedDataFig11/README.md>)

[Repository overview](<../README.md>) · [Current / historical numbering](<../12_Provenance/figure_index.csv>)
