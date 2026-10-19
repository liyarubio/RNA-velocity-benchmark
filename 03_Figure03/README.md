# Figure03 — 模拟与常规数据评价

Current figure: **Figure03**. Historical source ID: **Figure03**.

![Figure03](<figures/Figure03.png>)

## Evidence and reproduction

740 对主图配对；拓展图 698 对。逐数据集最大值为事后包络。

Data status: `available`. Analysis code: `historical_code_available`.

Portable saved-data checks: [commands and dependencies](<../09_Core/README.md>). Historical scripts retain their original environment and external-input requirements. Layout scripts are provenance, not standalone numerical generators.

## Files

| Role | File | Purpose |
|---|---|---|
| Figures | [Figure03.pdf](<figures/Figure03.pdf>) | 当前正式图件或预览 |
| Figures | [Figure03.png](<figures/Figure03.png>) | 当前正式图件或预览 |
| Figures | [Figure03_legend.pdf](<figures/Figure03_legend.pdf>) | 当前正式图件或预览 |
| Data | [Figure03.xlsx](<data/source_data/Figure03.xlsx>) | 源表、保存结果或输入元数据 |
| Data | [Figure3a_high_dim_cosine_similarity.csv](<data/figure3_detail/Figure3a_high_dim_cosine_similarity.csv>) | 源表、保存结果或输入元数据 |
| Data | [Figure3bc_low_dim_cosine_similarity.csv](<data/figure3_detail/Figure3bc_low_dim_cosine_similarity.csv>) | 源表、保存结果或输入元数据 |
| Data | [Figure3bc_paired_data.csv](<data/figure3_detail/Figure3bc_paired_data.csv>) | 源表、保存结果或输入元数据 |
| Data | [Figure3d_pca2_to_30_cosine_similarity.csv](<data/figure3_detail/Figure3d_pca2_to_30_cosine_similarity.csv>) | 源表、保存结果或输入元数据 |
| Data | [Figure3d_plot_data.csv](<data/figure3_detail/Figure3d_plot_data.csv>) | 源表、保存结果或输入元数据 |
| Data | [Figure3f_aligned_metrics.csv](<data/figure3_detail/Figure3f_aligned_metrics.csv>) | 源表、保存结果或输入元数据 |
| Data | [Figure3f_correlation.csv](<data/figure3_detail/Figure3f_correlation.csv>) | 源表、保存结果或输入元数据 |
| Data | [Figure3f_direction_metric.csv](<data/figure3_detail/Figure3f_direction_metric.csv>) | 源表、保存结果或输入元数据 |
| Data | [Figure3f_global_metric.csv](<data/figure3_detail/Figure3f_global_metric.csv>) | 源表、保存结果或输入元数据 |
| Data | [Figure3f_incluster_metric.csv](<data/figure3_detail/Figure3f_incluster_metric.csv>) | 源表、保存结果或输入元数据 |
| Data | [Figure3f_transition_matrix_cosine_similarity.csv](<data/figure3_detail/Figure3f_transition_matrix_cosine_similarity.csv>) | 源表、保存结果或输入元数据 |
| Data | [Figure3g_hubmap_global_performance.csv](<data/figure3_detail/Figure3g_hubmap_global_performance.csv>) | 源表、保存结果或输入元数据 |
| Data | [Figure3g_long_data.csv](<data/figure3_detail/Figure3g_long_data.csv>) | 源表、保存结果或输入元数据 |
| Data | [Figure3g_method_medians.csv](<data/figure3_detail/Figure3g_method_medians.csv>) | 源表、保存结果或输入元数据 |
| Data | [Figure3h_data10_global.csv](<data/figure3_detail/Figure3h_data10_global.csv>) | 源表、保存结果或输入元数据 |
| Data | [Figure3h_dataset10_global_performance.csv](<data/figure3_detail/Figure3h_dataset10_global_performance.csv>) | 源表、保存结果或输入元数据 |
| Data | [data_manifest.csv](<data/figure3_detail/data_manifest.csv>) | 源表、保存结果或输入元数据 |
| Code | [reproduce.py](<../09_Core/validation/reproduce.py>) | 保存数据的数值一致性检查；不运行模型 |
| Code | [plot_saved_panels.py](<../09_Core/visualization/plot_saved_panels.py>) | 保存数据的诊断图；不替代正式版布局 |
| Code | [Figure3g_01a_hubmap_global_eval_cell2fate.py](<code/historical/figure3/Figure3g_01a_hubmap_global_eval_cell2fate.py>) | 历史计算或绘图实现；保留原路径、环境及外部输入依赖 |
| Code | [Figure3g_01b_hubmap_global_eval_topicvelo.py](<code/historical/figure3/Figure3g_01b_hubmap_global_eval_topicvelo.py>) | 历史计算或绘图实现；保留原路径、环境及外部输入依赖 |
| Code | [Figure3h_01_dataset10_compute_global_metrics.py](<code/historical/figure3/Figure3h_01_dataset10_compute_global_metrics.py>) | 历史计算或绘图实现；保留原路径、环境及外部输入依赖 |
| Code | [config.py](<../09_Core/model_execution/Run_velocity/config.py>) | 历史计算或绘图实现；保留原路径、环境及外部输入依赖 |
| Code | [run_bivi.py](<../09_Core/model_execution/Run_velocity/run_bivi.py>) | 历史计算或绘图实现；保留原路径、环境及外部输入依赖 |
| Code | [run_cell2fate.py](<../09_Core/model_execution/Run_velocity/run_cell2fate.py>) | 历史计算或绘图实现；保留原路径、环境及外部输入依赖 |
| Code | [run_celldancer.py](<../09_Core/model_execution/Run_velocity/run_celldancer.py>) | 历史计算或绘图实现；保留原路径、环境及外部输入依赖 |
| Code | [run_deepkinet.py](<../09_Core/model_execution/Run_velocity/run_deepkinet.py>) | 历史计算或绘图实现；保留原路径、环境及外部输入依赖 |
| Code | [run_deepvelo_gcn.py](<../09_Core/model_execution/Run_velocity/run_deepvelo_gcn.py>) | 历史计算或绘图实现；保留原路径、环境及外部输入依赖 |
| Code | [run_deepvelo_vae.py](<../09_Core/model_execution/Run_velocity/run_deepvelo_vae.py>) | 历史计算或绘图实现；保留原路径、环境及外部输入依赖 |
| Code | [run_dynamo.py](<../09_Core/model_execution/Run_velocity/run_dynamo.py>) | 历史计算或绘图实现；保留原路径、环境及外部输入依赖 |
| Code | [run_latentvelo.py](<../09_Core/model_execution/Run_velocity/run_latentvelo.py>) | 历史计算或绘图实现；保留原路径、环境及外部输入依赖 |
| Code | [run_multivelo.py](<../09_Core/model_execution/Run_velocity/run_multivelo.py>) | 历史计算或绘图实现；保留原路径、环境及外部输入依赖 |
| Code | [run_phylovelo.py](<../09_Core/model_execution/Run_velocity/run_phylovelo.py>) | 历史计算或绘图实现；保留原路径、环境及外部输入依赖 |
| Code | [run_scvelo.py](<../09_Core/model_execution/Run_velocity/run_scvelo.py>) | 历史计算或绘图实现；保留原路径、环境及外部输入依赖 |
| Code | [run_sdevelo.py](<../09_Core/model_execution/Run_velocity/run_sdevelo.py>) | 历史计算或绘图实现；保留原路径、环境及外部输入依赖 |
| Code | [run_stt.py](<../09_Core/model_execution/Run_velocity/run_stt.py>) | 历史计算或绘图实现；保留原路径、环境及外部输入依赖 |
| Code | [run_tfvelo.py](<../09_Core/model_execution/Run_velocity/run_tfvelo.py>) | 历史计算或绘图实现；保留原路径、环境及外部输入依赖 |
| Code | [run_topicvelo.py](<../09_Core/model_execution/Run_velocity/run_topicvelo.py>) | 历史计算或绘图实现；保留原路径、环境及外部输入依赖 |
| Code | [run_unitvelo.py](<../09_Core/model_execution/Run_velocity/run_unitvelo.py>) | 历史计算或绘图实现；保留原路径、环境及外部输入依赖 |
| Code | [run_veloae.py](<../09_Core/model_execution/Run_velocity/run_veloae.py>) | 历史计算或绘图实现；保留原路径、环境及外部输入依赖 |
| Code | [run_velocyto.py](<../09_Core/model_execution/Run_velocity/run_velocyto.py>) | 历史计算或绘图实现；保留原路径、环境及外部输入依赖 |
| Code | [run_velovae.py](<../09_Core/model_execution/Run_velocity/run_velovae.py>) | 历史计算或绘图实现；保留原路径、环境及外部输入依赖 |
| Code | [run_velovi.py](<../09_Core/model_execution/Run_velocity/run_velovi.py>) | 历史计算或绘图实现；保留原路径、环境及外部输入依赖 |
| Code | [Run_velogif.py](<../09_Core/model_execution/Run_velogif.py>) | 历史计算或绘图实现；保留原路径、环境及外部输入依赖 |
| Code | [run_velogif_batch.py](<../09_Core/model_execution/run_velogif_batch.py>) | 历史计算或绘图实现；保留原路径、环境及外部输入依赖 |

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
