# ED06 — 真实数据 GDC

Current figure: **ED06**. Historical source ID: **ED07**.

![ED06](<figures/Extended_Data_Fig06.png>)

## Evidence and reproduction

与主图源值对齐；HuBMAP 18 配置，TopicVelo 无该任务有效值。

Data status: `available`. Analysis code: `historical_code_available`.

Portable saved-data checks: [commands and dependencies](<../../09_Core/README.md>). Historical scripts retain their original environment and external-input requirements. Layout scripts are provenance, not standalone numerical generators.

## Files

| Role | File | Purpose |
|---|---|---|
| Figures | [Extended_Data_Fig06.pdf](<figures/Extended_Data_Fig06.pdf>) | 当前正式图件或预览 |
| Figures | [Extended_Data_Fig06.png](<figures/Extended_Data_Fig06.png>) | 当前正式图件或预览 |
| Data | [Extended_Data_Fig07.xlsx](<data/source_data/Extended_Data_Fig07.xlsx>) | 源表、保存结果或输入元数据 |
| Code | [reproduce.py](<../../09_Core/validation/reproduce.py>) | 保存数据的数值一致性检查；不运行模型 |
| Code | [plot_saved_panels.py](<../../09_Core/visualization/plot_saved_panels.py>) | 保存数据的诊断图；不替代正式版布局 |
| Code | [redraw_ed5.py](<../../09_Core/shared/plots/redraw_ed5.py>) | 历史计算或绘图实现；保留原路径、环境及外部输入依赖 |
| Code | [redraw_extended.py](<../../09_Core/shared/plots/redraw_extended.py>) | 历史计算或绘图实现；保留原路径、环境及外部输入依赖 |
| Code | [build_figure_order_assets.py](<../../09_Core/layout/build_figure_order_assets.py>) | 当前图件排版、编号或图注处理；不生成新数值 |
| Code | [figure_order_1019.py](<../../09_Core/layout/figure_order_1019.py>) | 当前图件排版、编号或图注处理；不生成新数值 |
| Code | [figure2_legend_style.py](<../../09_Core/layout/figure2_legend_style.py>) | 当前图件排版、编号或图注处理；不生成新数值 |

## Related figures

[Figure03](<../../03_Figure03/README.md>)

[Repository overview](<../../README.md>) · [Current / historical numbering](<../../12_Provenance/figure_index.csv>)
