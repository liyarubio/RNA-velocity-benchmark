# Figure06 — 冻结综合分数和工具

Current figure: **Figure06**. Historical source ID: **Figure06**.

![Figure06](<figures/Figure06.png>)

## Evidence and reproduction

19 配置；0.6 准确性排名+0.4 可扩展性排名；稳定性不参与。与后续四领域有 8 项差异。

Data status: `available`. Analysis code: `historical_code_available`.

Portable saved-data checks: [commands and dependencies](<../09_Core/README.md>). Historical scripts retain their original environment and external-input requirements. Layout scripts are provenance, not standalone numerical generators.

## Files

| Role | File | Purpose |
|---|---|---|
| Figures | [Figure06.pdf](<figures/Figure06.pdf>) | 当前正式图件或预览 |
| Figures | [Figure06.png](<figures/Figure06.png>) | 当前正式图件或预览 |
| Data | [Figure06.xlsx](<data/source_data/Figure06.xlsx>) | 源表、保存结果或输入元数据 |
| Code | [reproduce.py](<../09_Core/validation/reproduce.py>) | 保存数据的数值一致性检查；不运行模型 |
| Code | [plot_saved_panels.py](<../09_Core/visualization/plot_saved_panels.py>) | 保存数据的诊断图；不替代正式版布局 |
| Code | [figure6a.py](<code/historical/Figure6_sources/code/panel_a/figure6a.py>) | 历史计算或绘图实现；保留原路径、环境及外部输入依赖 |

## Additional retained files

- [code/legacy/Visualization/Vis_summary/9_Overall_score.R](<code/legacy/Visualization/Vis_summary/9_Overall_score.R>)

## Related figures

[ED12](<../08_Extended_Data/ExtendedDataFig12/README.md>)

[Repository overview](<../README.md>) · [Current / historical numbering](<../12_Provenance/figure_index.csv>)
