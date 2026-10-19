# Figure05 — 空间、运行时间及 reads

Current figure: **Figure05**. Historical source ID: **Figure05**.

![Figure05](<figures/Figure05.png>)

## Evidence and reproduction

ST_small 608 位点；TopicVelo 源表修正至 0.4838869706599916，与最终图一致。

Data status: `available`. Analysis code: `historical_code_available`.

Portable saved-data checks: [commands and dependencies](<../09_Core/README.md>). Historical scripts retain their original environment and external-input requirements. Layout scripts are provenance, not standalone numerical generators.

## Files

| Role | File | Purpose |
|---|---|---|
| Figures | [Figure05.pdf](<figures/Figure05.pdf>) | 当前正式图件或预览 |
| Figures | [Figure05.png](<figures/Figure05.png>) | 当前正式图件或预览 |
| Figures | [Figure05_legend.pdf](<figures/Figure05_legend.pdf>) | 当前正式图件或预览 |
| Data | [Figure05.xlsx](<data/source_data/Figure05.xlsx>) | 源表、保存结果或输入元数据 |
| Data | [GDC_display_order.csv](<../10_Shared/data/spatial_completion/GDC_display_order.csv>) | 源表、保存结果或输入元数据 |
| Data | [ST_raw_status.json](<../10_Shared/data/spatial_completion/ST_raw_status.json>) | 源表、保存结果或输入元数据 |
| Data | [ST_small_status.json](<../10_Shared/data/spatial_completion/ST_small_status.json>) | 源表、保存结果或输入元数据 |
| Data | [TFvelo_evaluation_summary.csv](<../10_Shared/data/spatial_completion/TFvelo_evaluation_summary.csv>) | 源表、保存结果或输入元数据 |
| Data | [completion_audit.json](<../10_Shared/data/spatial_completion/completion_audit.json>) | 源表、保存结果或输入元数据 |
| Data | [evaluation_probe.json](<../10_Shared/data/spatial_completion/evaluation_probe.json>) | 源表、保存结果或输入元数据 |
| Data | [spatial_cells.csv](<../10_Shared/data/spatial_completion/spatial_cells.csv>) | 源表、保存结果或输入元数据 |
| Data | [spatial_directions.csv](<../10_Shared/data/spatial_completion/spatial_directions.csv>) | 源表、保存结果或输入元数据 |
| Data | [ST_raw.json](<../10_Shared/data/spatial_completion/topicvelo_retry/ST_raw.json>) | 源表、保存结果或输入元数据 |
| Data | [ST_raw_directions.csv](<../10_Shared/data/spatial_completion/topicvelo_retry/ST_raw_directions.csv>) | 源表、保存结果或输入元数据 |
| Data | [ST_raw_effective_transition_function.py](<../10_Shared/data/spatial_completion/topicvelo_retry/ST_raw_effective_transition_function.py>) | 源表、保存结果或输入元数据 |
| Data | [ST_raw_native_transition.npz](<../10_Shared/data/spatial_completion/topicvelo_retry/ST_raw_native_transition.npz>) | 源表、保存结果或输入元数据 |
| Data | [ST_raw_status.json](<../10_Shared/data/spatial_completion/topicvelo_retry/ST_raw_status.json>) | 源表、保存结果或输入元数据 |
| Data | [ST_small.json](<../10_Shared/data/spatial_completion/topicvelo_retry/ST_small.json>) | 源表、保存结果或输入元数据 |
| Data | [ST_small_directions.csv](<../10_Shared/data/spatial_completion/topicvelo_retry/ST_small_directions.csv>) | 源表、保存结果或输入元数据 |
| Data | [ST_small_effective_transition_function.py](<../10_Shared/data/spatial_completion/topicvelo_retry/ST_small_effective_transition_function.py>) | 源表、保存结果或输入元数据 |
| Data | [ST_small_native_transition.npz](<../10_Shared/data/spatial_completion/topicvelo_retry/ST_small_native_transition.npz>) | 源表、保存结果或输入元数据 |
| Data | [ST_small_status.json](<../10_Shared/data/spatial_completion/topicvelo_retry/ST_small_status.json>) | 源表、保存结果或输入元数据 |
| Code | [plot_saved_panels.py](<../09_Core/visualization/plot_saved_panels.py>) | 保存数据的诊断图；不替代正式版布局 |
| Code | [complete_topicvelo.py](<../09_Core/shared/spatial_completion/complete_topicvelo.py>) | 历史计算或绘图实现；保留原路径、环境及外部输入依赖 |
| Code | [evaluate_saved_outputs.py](<../09_Core/shared/spatial_completion/evaluate_saved_outputs.py>) | 历史计算或绘图实现；保留原路径、环境及外部输入依赖 |
| Code | [evaluate_topicvelo.py](<../09_Core/shared/spatial_topicvelo_retry/evaluate_topicvelo.py>) | 历史计算或绘图实现；保留原路径、环境及外部输入依赖 |
| Code | [retry_topicvelo.py](<../09_Core/shared/spatial_topicvelo_retry/retry_topicvelo.py>) | 历史计算或绘图实现；保留原路径、环境及外部输入依赖 |
| Code | [browser_check.py](<../09_Core/shared/stability_completion/browser_check.py>) | 历史计算或绘图实现；保留原路径、环境及外部输入依赖 |
| Code | [check_sparse_assembly.py](<../09_Core/shared/stability_completion/check_sparse_assembly.py>) | 历史计算或绘图实现；保留原路径、环境及外部输入依赖 |
| Code | [complete_topicvelo.py](<../09_Core/shared/stability_completion/complete_topicvelo.py>) | 历史计算或绘图实现；保留原路径、环境及外部输入依赖 |
| Code | [evaluate_stability.py](<../09_Core/shared/stability_completion/evaluate_stability.py>) | 历史计算或绘图实现；保留原路径、环境及外部输入依赖 |
| Code | [prepare_publication.py](<../09_Core/shared/stability_completion/prepare_publication.py>) | 历史计算或绘图实现；保留原路径、环境及外部输入依赖 |
| Code | [render_original_style.py](<../09_Core/shared/stability_completion/render_original_style.py>) | 历史计算或绘图实现；保留原路径、环境及外部输入依赖 |
| Code | [update_html.py](<../09_Core/shared/stability_completion/update_html.py>) | 历史计算或绘图实现；保留原路径、环境及外部输入依赖 |
| Code | [validate_delivery.py](<../09_Core/shared/stability_completion/validate_delivery.py>) | 历史计算或绘图实现；保留原路径、环境及外部输入依赖 |
| Code | [Figure5c_velocity_matrix_similarity_original.py](<code/historical/Figure5_sources/code/panel_c/Figure5c_velocity_matrix_similarity_original.py>) | 历史计算或绘图实现；保留原路径、环境及外部输入依赖 |

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
