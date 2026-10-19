# ED10 — 上下文配对一致性

Current figure: **ED10**. Historical source ID: **ED16**.

![ED10](<figures/Extended_Data_Fig10.png>)

## Evidence and reproduction

19×608=11,552 对空间结果；18 方法有基因空间分数；VeloAE 为 NA。

Data status: `available`. Analysis code: `historical_code_available`.

Portable saved-data checks: [commands and dependencies](<../../09_Core/README.md>). Historical scripts retain their original environment and external-input requirements. Layout scripts are provenance, not standalone numerical generators.

## Files

| Role | File | Purpose |
|---|---|---|
| Figures | [Extended_Data_Fig10.pdf](<figures/Extended_Data_Fig10.pdf>) | 当前正式图件或预览 |
| Figures | [Extended_Data_Fig10.png](<figures/Extended_Data_Fig10.png>) | 当前正式图件或预览 |
| Figures | [Extended_Data_Fig10_legend.pdf](<figures/Extended_Data_Fig10_legend.pdf>) | 当前正式图件或预览 |
| Data | [Extended_Data_Fig16.xlsx](<data/source_data/Extended_Data_Fig16.xlsx>) | 源表、保存结果或输入元数据 |
| Data | [DeepKINET_shared_genes.csv](<data/spatial_context_consistency/DeepKINET_shared_genes.csv>) | 源表、保存结果或输入元数据 |
| Data | [DeepVelo(GCN)_shared_genes.csv](<data/spatial_context_consistency/DeepVelo(GCN)_shared_genes.csv>) | 源表、保存结果或输入元数据 |
| Data | [DeepVelo(VAE)_shared_genes.csv](<data/spatial_context_consistency/DeepVelo(VAE)_shared_genes.csv>) | 源表、保存结果或输入元数据 |
| Data | [LatentVelo_shared_genes.csv](<data/spatial_context_consistency/LatentVelo_shared_genes.csv>) | 源表、保存结果或输入元数据 |
| Data | [SDEvelo_shared_genes.csv](<data/spatial_context_consistency/SDEvelo_shared_genes.csv>) | 源表、保存结果或输入元数据 |
| Data | [STT_shared_genes.csv](<data/spatial_context_consistency/STT_shared_genes.csv>) | 源表、保存结果或输入元数据 |
| Data | [TFvelo_shared_genes.csv](<data/spatial_context_consistency/TFvelo_shared_genes.csv>) | 源表、保存结果或输入元数据 |
| Data | [TopicVelo_shared_genes.csv](<data/spatial_context_consistency/TopicVelo_shared_genes.csv>) | 源表、保存结果或输入元数据 |
| Data | [UniTVelo_shared_genes.csv](<data/spatial_context_consistency/UniTVelo_shared_genes.csv>) | 源表、保存结果或输入元数据 |
| Data | [VeloVAE_shared_genes.csv](<data/spatial_context_consistency/VeloVAE_shared_genes.csv>) | 源表、保存结果或输入元数据 |
| Data | [audit.json](<data/spatial_context_consistency/audit.json>) | 源表、保存结果或输入元数据 |
| Data | [biVI_shared_genes.csv](<data/spatial_context_consistency/biVI_shared_genes.csv>) | 源表、保存结果或输入元数据 |
| Data | [by_cell_type.csv](<data/spatial_context_consistency/by_cell_type.csv>) | 源表、保存结果或输入元数据 |
| Data | [cell2fate_shared_genes.csv](<data/spatial_context_consistency/cell2fate_shared_genes.csv>) | 源表、保存结果或输入元数据 |
| Data | [cellDancer_shared_genes.csv](<data/spatial_context_consistency/cellDancer_shared_genes.csv>) | 源表、保存结果或输入元数据 |
| Data | [display_summary.csv](<data/spatial_context_consistency/display_summary.csv>) | 源表、保存结果或输入元数据 |
| Data | [dynamo_shared_genes.csv](<data/spatial_context_consistency/dynamo_shared_genes.csv>) | 源表、保存结果或输入元数据 |
| Data | [paired_spots.csv](<data/spatial_context_consistency/paired_spots.csv>) | 源表、保存结果或输入元数据 |
| Data | [preprocessing_audit.csv](<data/spatial_context_consistency/preprocessing_audit.csv>) | 源表、保存结果或输入元数据 |
| Data | [run_metadata.json](<data/spatial_context_consistency/run_metadata.json>) | 源表、保存结果或输入元数据 |
| Data | [scVelo(dyn)_shared_genes.csv](<data/spatial_context_consistency/scVelo(dyn)_shared_genes.csv>) | 源表、保存结果或输入元数据 |
| Data | [scVelo(sto)_shared_genes.csv](<data/spatial_context_consistency/scVelo(sto)_shared_genes.csv>) | 源表、保存结果或输入元数据 |
| Data | [summary.csv](<data/spatial_context_consistency/summary.csv>) | 源表、保存结果或输入元数据 |
| Data | [veloVI_shared_genes.csv](<data/spatial_context_consistency/veloVI_shared_genes.csv>) | 源表、保存结果或输入元数据 |
| Data | [velocyto_shared_genes.csv](<data/spatial_context_consistency/velocyto_shared_genes.csv>) | 源表、保存结果或输入元数据 |
| Code | [reproduce.py](<../../09_Core/validation/reproduce.py>) | 保存数据的数值一致性检查；不运行模型 |
| Code | [plot_saved_panels.py](<../../09_Core/visualization/plot_saved_panels.py>) | 保存数据的诊断图；不替代正式版布局 |
| Code | [analyze.py](<code/historical/spatial_context_consistency/analyze.py>) | 历史计算或绘图实现；保留原路径、环境及外部输入依赖 |
| Code | [plot_context_ed16.py](<code/historical/plots/plot_context_ed16.py>) | 历史计算或绘图实现；保留原路径、环境及外部输入依赖 |
| Code | [build_figure_order_assets.py](<../../09_Core/layout/build_figure_order_assets.py>) | 当前图件排版、编号或图注处理；不生成新数值 |
| Code | [figure_order_1019.py](<../../09_Core/layout/figure_order_1019.py>) | 当前图件排版、编号或图注处理；不生成新数值 |
| Code | [figure2_legend_style.py](<../../09_Core/layout/figure2_legend_style.py>) | 当前图件排版、编号或图注处理；不生成新数值 |

## Related figures

[Figure05](<../../05_Figure05/README.md>)

[Repository overview](<../../README.md>) · [Current / historical numbering](<../../12_Provenance/figure_index.csv>)
