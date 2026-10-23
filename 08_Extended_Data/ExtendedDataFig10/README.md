# ED10 — Paired consistency across spatial contexts

Current release: **Manuscript_1021**, including updates through **22 October 2026**. Historical numerical source ID: **ED16**.

20 methods x 608 paired spots after adding TopoVelo. TopoVelo has 47 shared genes and 606 valid pairs; two undefined pairs remain missing. VeloAE gene-space scores remain unavailable.

![ED10](<figures/Extended_Data_Fig10.png>)

## Individual panel code

[One plot per file: current entrypoints and limits](code/README.md).

## Current files

| Role | File |
|---|---|
| results | [figures/Extended_Data_Fig10.pdf](<figures/Extended_Data_Fig10.pdf>) |
| results | [figures/Extended_Data_Fig10.png](<figures/Extended_Data_Fig10.png>) |
| results | [figures/Extended_Data_Fig10.svg](<figures/Extended_Data_Fig10.svg>) |
| results | [figures/legends/Extended_Data_Fig10_legend.pdf](<figures/legends/Extended_Data_Fig10_legend.pdf>) |
| data | [data/source_data/Extended_Data_Fig16.xlsx](<data/source_data/Extended_Data_Fig16.xlsx>) |
| data | [data/spatial_context_consistency/DeepKINET_shared_genes.csv](<data/spatial_context_consistency/DeepKINET_shared_genes.csv>) |
| data | [data/spatial_context_consistency/DeepVelo(GCN)_shared_genes.csv](<data/spatial_context_consistency/DeepVelo(GCN)_shared_genes.csv>) |
| data | [data/spatial_context_consistency/DeepVelo(VAE)_shared_genes.csv](<data/spatial_context_consistency/DeepVelo(VAE)_shared_genes.csv>) |
| data | [data/spatial_context_consistency/LatentVelo_shared_genes.csv](<data/spatial_context_consistency/LatentVelo_shared_genes.csv>) |
| data | [data/spatial_context_consistency/SDEvelo_shared_genes.csv](<data/spatial_context_consistency/SDEvelo_shared_genes.csv>) |
| data | [data/spatial_context_consistency/STT_shared_genes.csv](<data/spatial_context_consistency/STT_shared_genes.csv>) |
| data | [data/spatial_context_consistency/TFvelo_shared_genes.csv](<data/spatial_context_consistency/TFvelo_shared_genes.csv>) |
| data | [data/spatial_context_consistency/TopicVelo_shared_genes.csv](<data/spatial_context_consistency/TopicVelo_shared_genes.csv>) |
| data | [../../10_Shared/data/spatial_context_consistency/TopoVelo_shared_genes.csv](<../../10_Shared/data/spatial_context_consistency/TopoVelo_shared_genes.csv>) |
| data | [data/spatial_context_consistency/UniTVelo_shared_genes.csv](<data/spatial_context_consistency/UniTVelo_shared_genes.csv>) |
| data | [data/spatial_context_consistency/VeloVAE_shared_genes.csv](<data/spatial_context_consistency/VeloVAE_shared_genes.csv>) |
| data | [data/spatial_context_consistency/audit.json](<data/spatial_context_consistency/audit.json>) |
| data | [data/spatial_context_consistency/biVI_shared_genes.csv](<data/spatial_context_consistency/biVI_shared_genes.csv>) |
| data | [data/spatial_context_consistency/by_cell_type.csv](<data/spatial_context_consistency/by_cell_type.csv>) |
| data | [data/spatial_context_consistency/cell2fate_shared_genes.csv](<data/spatial_context_consistency/cell2fate_shared_genes.csv>) |
| data | [data/spatial_context_consistency/cellDancer_shared_genes.csv](<data/spatial_context_consistency/cellDancer_shared_genes.csv>) |
| data | [data/spatial_context_consistency/display_summary.csv](<data/spatial_context_consistency/display_summary.csv>) |
| data | [data/spatial_context_consistency/dynamo_shared_genes.csv](<data/spatial_context_consistency/dynamo_shared_genes.csv>) |
| data | [data/spatial_context_consistency/paired_spots.csv](<data/spatial_context_consistency/paired_spots.csv>) |
| data | [data/spatial_context_consistency/preprocessing_audit.csv](<data/spatial_context_consistency/preprocessing_audit.csv>) |
| data | [data/spatial_context_consistency/run_metadata.json](<data/spatial_context_consistency/run_metadata.json>) |
| data | [data/spatial_context_consistency/scVelo(dyn)_shared_genes.csv](<data/spatial_context_consistency/scVelo(dyn)_shared_genes.csv>) |
| data | [data/spatial_context_consistency/scVelo(sto)_shared_genes.csv](<data/spatial_context_consistency/scVelo(sto)_shared_genes.csv>) |
| data | [data/spatial_context_consistency/summary.csv](<data/spatial_context_consistency/summary.csv>) |
| data | [data/spatial_context_consistency/veloVI_shared_genes.csv](<data/spatial_context_consistency/veloVI_shared_genes.csv>) |
| data | [data/spatial_context_consistency/velocyto_shared_genes.csv](<data/spatial_context_consistency/velocyto_shared_genes.csv>) |
| data | [../../10_Shared/data/topovelo_historical_20261021/README.md](<../../10_Shared/data/topovelo_historical_20261021/README.md>) |
| data | [../../10_Shared/data/topovelo_historical_20261021/derived/ST_raw_directions.csv](<../../10_Shared/data/topovelo_historical_20261021/derived/ST_raw_directions.csv>) |
| data | [../../10_Shared/data/topovelo_historical_20261021/derived/ST_small_directions.csv](<../../10_Shared/data/topovelo_historical_20261021/derived/ST_small_directions.csv>) |
| data | [../../10_Shared/data/topovelo_historical_20261021/derived/TopoVelo_field.pdf](<../../10_Shared/data/topovelo_historical_20261021/derived/TopoVelo_field.pdf>) |
| data | [../../10_Shared/data/topovelo_historical_20261021/derived/TopoVelo_shared_genes.csv](<../../10_Shared/data/topovelo_historical_20261021/derived/TopoVelo_shared_genes.csv>) |
| data | [../../10_Shared/data/topovelo_historical_20261021/derived/by_cell_type.csv](<../../10_Shared/data/topovelo_historical_20261021/derived/by_cell_type.csv>) |
| data | [../../10_Shared/data/topovelo_historical_20261021/derived/evaluation.json](<../../10_Shared/data/topovelo_historical_20261021/derived/evaluation.json>) |
| data | [../../10_Shared/data/topovelo_historical_20261021/derived/paired_spots.csv](<../../10_Shared/data/topovelo_historical_20261021/derived/paired_spots.csv>) |
| data | [../../10_Shared/data/topovelo_historical_20261021/derived/summary.csv](<../../10_Shared/data/topovelo_historical_20261021/derived/summary.csv>) |
| data | [../../10_Shared/data/topovelo_historical_20261021/provenance.json](<../../10_Shared/data/topovelo_historical_20261021/provenance.json>) |
| code | [code/historical/spatial_context_consistency/analyze.py](<code/historical/spatial_context_consistency/analyze.py>) |
| code | [code/historical/plots/plot_context_ed16.py](<code/historical/plots/plot_context_ed16.py>) |
| code | [../../09_Core/shared/latest/derive_topovelo_historical_20261021.py](<../../09_Core/shared/latest/derive_topovelo_historical_20261021.py>) |
| code | [../../09_Core/shared/latest/topovelo_historical_evaluation.py](<../../09_Core/shared/latest/topovelo_historical_evaluation.py>) |

## Code and validation

[Portable checks and workflow boundaries](<../../09_Core/README.md>). Saved-data checks do not rerun models. Historical scripts require their original environments and external inputs.

[Retained figure-specific code](<code>)

[Repository overview](<../../README.md>) · [Release notes](<../../12_Provenance/RELEASE_1021.md>)
