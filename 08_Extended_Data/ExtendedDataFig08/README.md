# ED08 — scLT 谱系参照

Current figure: **ED08**. Historical source ID: **ED14**.

![ED08](<figures/Extended_Data_Fig08.png>)

## Evidence and reproduction

1,756 参照细胞、19 配置；五个克隆示例的原成员与随机选择记录不完整。

Data status: `available`. Analysis code: `historical_code_available`.

Portable saved-data checks: [commands and dependencies](<../../09_Core/README.md>). Historical scripts retain their original environment and external-input requirements. Layout scripts are provenance, not standalone numerical generators.

## Files

| Role | File | Purpose |
|---|---|---|
| Figures | [Extended_Data_Fig08.pdf](<figures/Extended_Data_Fig08.pdf>) | 当前正式图件或预览 |
| Figures | [Extended_Data_Fig08.png](<figures/Extended_Data_Fig08.png>) | 当前正式图件或预览 |
| Figures | [Extended_Data_Fig08_legend.pdf](<figures/Extended_Data_Fig08_legend.pdf>) | 当前正式图件或预览 |
| Data | [Extended_Data_Fig14.xlsx](<data/source_data/Extended_Data_Fig14.xlsx>) | 源表、保存结果或输入元数据 |
| Data | [scLT_projected_directions.csv](<data/reference_restoration/scLT_projected_directions.csv>) | 源表、保存结果或输入元数据 |
| Code | [reproduce.py](<../../09_Core/validation/reproduce.py>) | 保存数据的数值一致性检查；不运行模型 |
| Code | [plot_saved_panels.py](<../../09_Core/visualization/plot_saved_panels.py>) | 保存数据的诊断图；不替代正式版布局 |
| Code | [restore_reference_figures.py](<../../09_Core/shared/plots/restore_reference_figures.py>) | 历史计算或绘图实现；保留原路径、环境及外部输入依赖 |
| Code | [evaluate_reference.py](<../../09_Core/reference_evaluation/scLT_scNT_velogif/evaluate_reference.py>) | 历史计算或绘图实现；保留原路径、环境及外部输入依赖 |
| Code | [finalize.py](<../../09_Core/reference_evaluation/scLT_scNT_velogif/finalize.py>) | 历史计算或绘图实现；保留原路径、环境及外部输入依赖 |
| Code | [prepare.py](<../../09_Core/reference_evaluation/scLT_scNT_velogif/prepare.py>) | 历史计算或绘图实现；保留原路径、环境及外部输入依赖 |
| Code | [run_all.py](<../../09_Core/reference_evaluation/scLT_scNT_velogif/run_all.py>) | 历史计算或绘图实现；保留原路径、环境及外部输入依赖 |
| Code | [test_reference.py](<../../09_Core/reference_evaluation/scLT_scNT_velogif/test_reference.py>) | 历史计算或绘图实现；保留原路径、环境及外部输入依赖 |
| Code | [restore_tfvelo_projection.py](<../../09_Core/shared/restore_tfvelo_projection.py>) | 历史计算或绘图实现；保留原路径、环境及外部输入依赖 |
| Code | [build_figure_order_assets.py](<../../09_Core/layout/build_figure_order_assets.py>) | 当前图件排版、编号或图注处理；不生成新数值 |
| Code | [figure_order_1019.py](<../../09_Core/layout/figure_order_1019.py>) | 当前图件排版、编号或图注处理；不生成新数值 |
| Code | [figure2_legend_style.py](<../../09_Core/layout/figure2_legend_style.py>) | 当前图件排版、编号或图注处理；不生成新数值 |
| Data | [scLT_TFvelo_restored_projection.csv](<data/reference_restoration/scLT_TFvelo_restored_projection.csv>) | 源表、保存结果或输入元数据 |
| Data | [scLT_method_summary.csv](<data/reference_restoration/scLT_method_summary.csv>) | 源表、保存结果或输入元数据 |
| Data | [scLT_per_cell_cosine.csv](<data/reference_restoration/scLT_per_cell_cosine.csv>) | 源表、保存结果或输入元数据 |
| Data | [scLT_reference_cells.csv](<data/reference_restoration/scLT_reference_cells.csv>) | 源表、保存结果或输入元数据 |
| Data | [scLT_saved_expected_targets.csv](<data/reference_restoration/scLT_saved_expected_targets.csv>) | 源表、保存结果或输入元数据 |
| Data | [scLT_timepoints.csv](<data/reference_restoration/scLT_timepoints.csv>) | 源表、保存结果或输入元数据 |

## Additional retained files

- [code/legacy/Visualization/Vis_velocity_stream/Vis_scLT.ipynb](<code/legacy/Visualization/Vis_velocity_stream/Vis_scLT.ipynb>)
- [code/legacy/scLT/1_process/1_run_velocyto.sh](<code/legacy/scLT/1_process/1_run_velocyto.sh>)
- [code/legacy/scLT/1_process/2_add_velocyto.R](<code/legacy/scLT/1_process/2_add_velocyto.R>)
- [code/legacy/scLT/1_process/3_build_cell_pairwise.py](<code/legacy/scLT/1_process/3_build_cell_pairwise.py>)
- [code/legacy/scLT/1_process/4_process_and_bulid_clone_based_velocity.ipynb](<code/legacy/scLT/1_process/4_process_and_bulid_clone_based_velocity.ipynb>)
- [code/legacy/scLT/2_Run_velocity/10_velocyto.ipynb](<code/legacy/scLT/2_Run_velocity/10_velocyto.ipynb>)
- [code/legacy/scLT/2_Run_velocity/11_veloAE.ipynb](<code/legacy/scLT/2_Run_velocity/11_veloAE.ipynb>)
- [code/legacy/scLT/2_Run_velocity/12_veloVAE.ipynb](<code/legacy/scLT/2_Run_velocity/12_veloVAE.ipynb>)
- [code/legacy/scLT/2_Run_velocity/13_STT.ipynb](<code/legacy/scLT/2_Run_velocity/13_STT.ipynb>)
- [code/legacy/scLT/2_Run_velocity/1_scvelo_dynamic.ipynb](<code/legacy/scLT/2_Run_velocity/1_scvelo_dynamic.ipynb>)
- [code/legacy/scLT/2_Run_velocity/1_scvelo_stochastic.ipynb](<code/legacy/scLT/2_Run_velocity/1_scvelo_stochastic.ipynb>)
- [code/legacy/scLT/2_Run_velocity/2_TFvelo.py](<code/legacy/scLT/2_Run_velocity/2_TFvelo.py>)
- [code/legacy/scLT/2_Run_velocity/3_celldancer.ipynb](<code/legacy/scLT/2_Run_velocity/3_celldancer.ipynb>)
- [code/legacy/scLT/2_Run_velocity/4_unitvelo.ipynb](<code/legacy/scLT/2_Run_velocity/4_unitvelo.ipynb>)
- [code/legacy/scLT/2_Run_velocity/5_deepveloSA.ipynb](<code/legacy/scLT/2_Run_velocity/5_deepveloSA.ipynb>)
- [code/legacy/scLT/2_Run_velocity/6_deepveloGB.ipynb](<code/legacy/scLT/2_Run_velocity/6_deepveloGB.ipynb>)
- [code/legacy/scLT/2_Run_velocity/7_veloVI.ipynb](<code/legacy/scLT/2_Run_velocity/7_veloVI.ipynb>)
- [code/legacy/scLT/2_Run_velocity/8_latentvelo.ipynb](<code/legacy/scLT/2_Run_velocity/8_latentvelo.ipynb>)
- [code/legacy/scLT/2_Run_velocity/9_dynamo_deterministic.ipynb](<code/legacy/scLT/2_Run_velocity/9_dynamo_deterministic.ipynb>)
- [code/legacy/scLT/2_Run_velocity/9_dynamo_stochastic.ipynb](<code/legacy/scLT/2_Run_velocity/9_dynamo_stochastic.ipynb>)
- [code/legacy/scLT/3_eval_with_clone.ipynb](<code/legacy/scLT/3_eval_with_clone.ipynb>)

## Related figures

[Figure04](<../../04_Figure04/README.md>)

[Repository overview](<../../README.md>) · [Current / historical numbering](<../../12_Provenance/figure_index.csv>)
