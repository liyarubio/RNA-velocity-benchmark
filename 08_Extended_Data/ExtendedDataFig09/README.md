# ED09 — ST_raw 空间评价

Current figure: **ED09**. Historical source ID: **ED15**.

![ED09](<figures/Extended_Data_Fig09.png>)

## Evidence and reproduction

1,292 位点、五条关系；与主图 5a 的 ST_small 区分。

Data status: `available`. Analysis code: `historical_code_available`.

Portable saved-data checks: [commands and dependencies](<../../09_Core/README.md>). Historical scripts retain their original environment and external-input requirements. Layout scripts are provenance, not standalone numerical generators.

## Files

| Role | File | Purpose |
|---|---|---|
| Figures | [Extended_Data_Fig09.pdf](<figures/Extended_Data_Fig09.pdf>) | 当前正式图件或预览 |
| Figures | [Extended_Data_Fig09.png](<figures/Extended_Data_Fig09.png>) | 当前正式图件或预览 |
| Figures | [Extended_Data_Fig09_legend.pdf](<figures/Extended_Data_Fig09_legend.pdf>) | 当前正式图件或预览 |
| Data | [Extended_Data_Fig15.xlsx](<data/source_data/Extended_Data_Fig15.xlsx>) | 源表、保存结果或输入元数据 |
| Data | [GDC_display_order.csv](<../../10_Shared/data/spatial_completion/GDC_display_order.csv>) | 源表、保存结果或输入元数据 |
| Data | [ST_raw_status.json](<../../10_Shared/data/spatial_completion/ST_raw_status.json>) | 源表、保存结果或输入元数据 |
| Data | [ST_small_status.json](<../../10_Shared/data/spatial_completion/ST_small_status.json>) | 源表、保存结果或输入元数据 |
| Data | [TFvelo_evaluation_summary.csv](<../../10_Shared/data/spatial_completion/TFvelo_evaluation_summary.csv>) | 源表、保存结果或输入元数据 |
| Data | [completion_audit.json](<../../10_Shared/data/spatial_completion/completion_audit.json>) | 源表、保存结果或输入元数据 |
| Data | [evaluation_probe.json](<../../10_Shared/data/spatial_completion/evaluation_probe.json>) | 源表、保存结果或输入元数据 |
| Data | [spatial_cells.csv](<../../10_Shared/data/spatial_completion/spatial_cells.csv>) | 源表、保存结果或输入元数据 |
| Data | [spatial_directions.csv](<../../10_Shared/data/spatial_completion/spatial_directions.csv>) | 源表、保存结果或输入元数据 |
| Data | [ST_raw.json](<../../10_Shared/data/spatial_completion/topicvelo_retry/ST_raw.json>) | 源表、保存结果或输入元数据 |
| Data | [ST_raw_directions.csv](<../../10_Shared/data/spatial_completion/topicvelo_retry/ST_raw_directions.csv>) | 源表、保存结果或输入元数据 |
| Data | [ST_raw_effective_transition_function.py](<../../10_Shared/data/spatial_completion/topicvelo_retry/ST_raw_effective_transition_function.py>) | 源表、保存结果或输入元数据 |
| Data | [ST_raw_native_transition.npz](<../../10_Shared/data/spatial_completion/topicvelo_retry/ST_raw_native_transition.npz>) | 源表、保存结果或输入元数据 |
| Data | [ST_raw_status.json](<../../10_Shared/data/spatial_completion/topicvelo_retry/ST_raw_status.json>) | 源表、保存结果或输入元数据 |
| Data | [ST_small.json](<../../10_Shared/data/spatial_completion/topicvelo_retry/ST_small.json>) | 源表、保存结果或输入元数据 |
| Data | [ST_small_directions.csv](<../../10_Shared/data/spatial_completion/topicvelo_retry/ST_small_directions.csv>) | 源表、保存结果或输入元数据 |
| Data | [ST_small_effective_transition_function.py](<../../10_Shared/data/spatial_completion/topicvelo_retry/ST_small_effective_transition_function.py>) | 源表、保存结果或输入元数据 |
| Data | [ST_small_native_transition.npz](<../../10_Shared/data/spatial_completion/topicvelo_retry/ST_small_native_transition.npz>) | 源表、保存结果或输入元数据 |
| Data | [ST_small_status.json](<../../10_Shared/data/spatial_completion/topicvelo_retry/ST_small_status.json>) | 源表、保存结果或输入元数据 |
| Code | [plot_saved_panels.py](<../../09_Core/visualization/plot_saved_panels.py>) | 保存数据的诊断图；不替代正式版布局 |
| Code | [complete_topicvelo.py](<../../09_Core/shared/spatial_completion/complete_topicvelo.py>) | 历史计算或绘图实现；保留原路径、环境及外部输入依赖 |
| Code | [evaluate_saved_outputs.py](<../../09_Core/shared/spatial_completion/evaluate_saved_outputs.py>) | 历史计算或绘图实现；保留原路径、环境及外部输入依赖 |
| Code | [evaluate_topicvelo.py](<../../09_Core/shared/spatial_topicvelo_retry/evaluate_topicvelo.py>) | 历史计算或绘图实现；保留原路径、环境及外部输入依赖 |
| Code | [retry_topicvelo.py](<../../09_Core/shared/spatial_topicvelo_retry/retry_topicvelo.py>) | 历史计算或绘图实现；保留原路径、环境及外部输入依赖 |
| Code | [browser_check.py](<../../09_Core/shared/stability_completion/browser_check.py>) | 历史计算或绘图实现；保留原路径、环境及外部输入依赖 |
| Code | [check_sparse_assembly.py](<../../09_Core/shared/stability_completion/check_sparse_assembly.py>) | 历史计算或绘图实现；保留原路径、环境及外部输入依赖 |
| Code | [complete_topicvelo.py](<../../09_Core/shared/stability_completion/complete_topicvelo.py>) | 历史计算或绘图实现；保留原路径、环境及外部输入依赖 |
| Code | [evaluate_stability.py](<../../09_Core/shared/stability_completion/evaluate_stability.py>) | 历史计算或绘图实现；保留原路径、环境及外部输入依赖 |
| Code | [prepare_publication.py](<../../09_Core/shared/stability_completion/prepare_publication.py>) | 历史计算或绘图实现；保留原路径、环境及外部输入依赖 |
| Code | [render_original_style.py](<../../09_Core/shared/stability_completion/render_original_style.py>) | 历史计算或绘图实现；保留原路径、环境及外部输入依赖 |
| Code | [update_html.py](<../../09_Core/shared/stability_completion/update_html.py>) | 历史计算或绘图实现；保留原路径、环境及外部输入依赖 |
| Code | [validate_delivery.py](<../../09_Core/shared/stability_completion/validate_delivery.py>) | 历史计算或绘图实现；保留原路径、环境及外部输入依赖 |
| Code | [redraw_ed15_reference.py](<code/historical/plots/redraw_ed15_reference.py>) | 历史计算或绘图实现；保留原路径、环境及外部输入依赖 |
| Code | [build_figure_order_assets.py](<../../09_Core/layout/build_figure_order_assets.py>) | 当前图件排版、编号或图注处理；不生成新数值 |
| Code | [figure_order_1019.py](<../../09_Core/layout/figure_order_1019.py>) | 当前图件排版、编号或图注处理；不生成新数值 |
| Code | [figure2_legend_style.py](<../../09_Core/layout/figure2_legend_style.py>) | 当前图件排版、编号或图注处理；不生成新数值 |

## Related figures

[Figure05](<../../05_Figure05/README.md>)

[Repository overview](<../../README.md>) · [Current / historical numbering](<../../12_Provenance/figure_index.csv>)
