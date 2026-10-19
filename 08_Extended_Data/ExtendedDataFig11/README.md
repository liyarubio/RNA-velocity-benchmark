# ED11 — KPTracer

Current figure: **ED11**. Historical source ID: **ED17**.

![ED11](<figures/Extended_Data_Fig11.png>)

## Evidence and reproduction

754 细胞、四条 corrected eval1 关系、分母 706；核对 TFvelo 和 TopicVelo 补全输出。

Data status: `available`. Analysis code: `historical_code_available`.

Portable saved-data checks: [commands and dependencies](<../../09_Core/README.md>). Historical scripts retain their original environment and external-input requirements. Layout scripts are provenance, not standalone numerical generators.

## Files

| Role | File | Purpose |
|---|---|---|
| Figures | [Extended_Data_Fig11.pdf](<figures/Extended_Data_Fig11.pdf>) | 当前正式图件或预览 |
| Figures | [Extended_Data_Fig11.png](<figures/Extended_Data_Fig11.png>) | 当前正式图件或预览 |
| Figures | [Extended_Data_Fig11_legend.pdf](<figures/Extended_Data_Fig11_legend.pdf>) | 当前正式图件或预览 |
| Data | [Extended_Data_Fig17.xlsx](<data/source_data/Extended_Data_Fig17.xlsx>) | 源表、保存结果或输入元数据 |
| Data | [GDC_display_order.csv](<data/kptracer_completion/GDC_display_order.csv>) | 源表、保存结果或输入元数据 |
| Data | [code_manifest.json](<data/kptracer_completion/code_manifest.json>) | 源表、保存结果或输入元数据 |
| Data | [completion_audit.json](<data/kptracer_completion/completion_audit.json>) | 源表、保存结果或输入元数据 |
| Data | [control.json](<data/kptracer_completion/control.json>) | 源表、保存结果或输入元数据 |
| Data | [control_directions.csv](<data/kptracer_completion/control_directions.csv>) | 源表、保存结果或输入元数据 |
| Data | [control_unit_scores.csv](<data/kptracer_completion/control_unit_scores.csv>) | 源表、保存结果或输入元数据 |
| Data | [execution.json](<data/kptracer_completion/execution.json>) | 源表、保存结果或输入元数据 |
| Data | [input_manifest.json](<data/kptracer_completion/input_manifest.json>) | 源表、保存结果或输入元数据 |
| Data | [native_topicvelo_transition.npz](<data/kptracer_completion/native_topicvelo_transition.npz>) | 源表、保存结果或输入元数据 |
| Data | [tfvelo.json](<data/kptracer_completion/tfvelo.json>) | 源表、保存结果或输入元数据 |
| Data | [tfvelo_directions.csv](<data/kptracer_completion/tfvelo_directions.csv>) | 源表、保存结果或输入元数据 |
| Data | [tfvelo_gene_mapping.csv](<data/kptracer_completion/tfvelo_gene_mapping.csv>) | 源表、保存结果或输入元数据 |
| Data | [tfvelo_resources.json](<data/kptracer_completion/tfvelo_resources.json>) | 源表、保存结果或输入元数据 |
| Data | [tfvelo_revision.json](<data/kptracer_completion/tfvelo_revision.json>) | 源表、保存结果或输入元数据 |
| Data | [tfvelo_unit_scores.csv](<data/kptracer_completion/tfvelo_unit_scores.csv>) | 源表、保存结果或输入元数据 |
| Data | [tfvelo_validation.json](<data/kptracer_completion/tfvelo_validation.json>) | 源表、保存结果或输入元数据 |
| Data | [topicvelo.json](<data/kptracer_completion/topicvelo.json>) | 源表、保存结果或输入元数据 |
| Data | [topicvelo_directions.csv](<data/kptracer_completion/topicvelo_directions.csv>) | 源表、保存结果或输入元数据 |
| Data | [topicvelo_inference_status.json](<data/kptracer_completion/topicvelo_inference_status.json>) | 源表、保存结果或输入元数据 |
| Data | [topicvelo_unit_scores.csv](<data/kptracer_completion/topicvelo_unit_scores.csv>) | 源表、保存结果或输入元数据 |
| Data | [topicvelo_validation.json](<data/kptracer_completion/topicvelo_validation.json>) | 源表、保存结果或输入元数据 |
| Code | [reproduce.py](<../../09_Core/validation/reproduce.py>) | 保存数据的数值一致性检查；不运行模型 |
| Code | [plot_saved_panels.py](<../../09_Core/visualization/plot_saved_panels.py>) | 保存数据的诊断图；不替代正式版布局 |
| Code | [complete_topicvelo.py](<code/historical/kptracer_completion/complete_topicvelo.py>) | 历史计算或绘图实现；保留原路径、环境及外部输入依赖 |
| Code | [config.py](<code/historical/kptracer_completion/config.py>) | 历史计算或绘图实现；保留原路径、环境及外部输入依赖 |
| Code | [evaluate.py](<code/historical/kptracer_completion/evaluate.py>) | 历史计算或绘图实现；保留原路径、环境及外部输入依赖 |
| Code | [evaluation_metrics.py](<code/historical/kptracer_completion/evaluation_metrics.py>) | 历史计算或绘图实现；保留原路径、环境及外部输入依赖 |
| Code | [import_to_manuscript.py](<code/historical/kptracer_completion/import_to_manuscript.py>) | 历史计算或绘图实现；保留原路径、环境及外部输入依赖 |
| Code | [run_tfvelo.py](<code/historical/kptracer_completion/run_tfvelo.py>) | 历史计算或绘图实现；保留原路径、环境及外部输入依赖 |
| Code | [validate_tfvelo.py](<code/historical/kptracer_completion/validate_tfvelo.py>) | 历史计算或绘图实现；保留原路径、环境及外部输入依赖 |
| Code | [compose_ed16_single_page.py](<code/historical/plots/compose_ed16_single_page.py>) | 历史计算或绘图实现；保留原路径、环境及外部输入依赖 |
| Code | [redraw_ed16_reference.py](<code/historical/plots/redraw_ed16_reference.py>) | 历史计算或绘图实现；保留原路径、环境及外部输入依赖 |
| Code | [ed16_literature_schematic.py](<code/historical/plots/ed16_literature_schematic.py>) | 历史计算或绘图实现；保留原路径、环境及外部输入依赖 |
| Code | [ed16_aligned_fields.py](<code/historical/plots/ed16_aligned_fields.py>) | 历史计算或绘图实现；保留原路径、环境及外部输入依赖 |
| Code | [complete_ed16_fields.py](<code/historical/plots/complete_ed16_fields.py>) | 历史计算或绘图实现；保留原路径、环境及外部输入依赖 |
| Code | [spatial_display_names.py](<code/historical/plots/spatial_display_names.py>) | 历史计算或绘图实现；保留原路径、环境及外部输入依赖 |
| Code | [build_figure_order_assets.py](<../../09_Core/layout/build_figure_order_assets.py>) | 当前图件排版、编号或图注处理；不生成新数值 |
| Code | [figure_order_1019.py](<../../09_Core/layout/figure_order_1019.py>) | 当前图件排版、编号或图注处理；不生成新数值 |
| Code | [figure2_legend_style.py](<../../09_Core/layout/figure2_legend_style.py>) | 当前图件排版、编号或图注处理；不生成新数值 |

## Related figures

[Figure05](<../../05_Figure05/README.md>)

[Repository overview](<../../README.md>) · [Current / historical numbering](<../../12_Provenance/figure_index.csv>)
