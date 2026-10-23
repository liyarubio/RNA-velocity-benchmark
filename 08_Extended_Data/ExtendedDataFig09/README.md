# ED09 — ST_raw spatial evaluation

Current release: **Manuscript_1021**, including updates through **22 October 2026**. Historical numerical source ID: **ED15**.

1,292 spots and five reference relationships. Archived TopoVelo is included with a 608-cell reference denominator; its training configuration is unavailable.

![ED09](<figures/Extended_Data_Fig09.png>)

## Current files

| Role | File |
|---|---|
| results | [figures/Extended_Data_Fig09.pdf](<figures/Extended_Data_Fig09.pdf>) |
| results | [figures/Extended_Data_Fig09.png](<figures/Extended_Data_Fig09.png>) |
| results | [figures/Extended_Data_Fig09.svg](<figures/Extended_Data_Fig09.svg>) |
| results | [figures/legends/Extended_Data_Fig09_legend.pdf](<figures/legends/Extended_Data_Fig09_legend.pdf>) |
| data | [data/source_data/Extended_Data_Fig15.xlsx](<data/source_data/Extended_Data_Fig15.xlsx>) |
| data | [../../10_Shared/data/spatial_completion/GDC_display_order.csv](<../../10_Shared/data/spatial_completion/GDC_display_order.csv>) |
| data | [../../10_Shared/data/spatial_completion/ST_raw_status.json](<../../10_Shared/data/spatial_completion/ST_raw_status.json>) |
| data | [../../10_Shared/data/spatial_completion/ST_small_status.json](<../../10_Shared/data/spatial_completion/ST_small_status.json>) |
| data | [../../10_Shared/data/spatial_completion/TFvelo_evaluation_summary.csv](<../../10_Shared/data/spatial_completion/TFvelo_evaluation_summary.csv>) |
| data | [../../10_Shared/data/spatial_completion/completion_audit.json](<../../10_Shared/data/spatial_completion/completion_audit.json>) |
| data | [../../10_Shared/data/spatial_completion/evaluation_probe.json](<../../10_Shared/data/spatial_completion/evaluation_probe.json>) |
| data | [../../10_Shared/data/spatial_completion/spatial_cells.csv](<../../10_Shared/data/spatial_completion/spatial_cells.csv>) |
| data | [../../10_Shared/data/spatial_completion/spatial_directions.csv](<../../10_Shared/data/spatial_completion/spatial_directions.csv>) |
| data | [../../10_Shared/data/spatial_completion/topicvelo_retry/ST_raw.json](<../../10_Shared/data/spatial_completion/topicvelo_retry/ST_raw.json>) |
| data | [../../10_Shared/data/spatial_completion/topicvelo_retry/ST_raw_directions.csv](<../../10_Shared/data/spatial_completion/topicvelo_retry/ST_raw_directions.csv>) |
| data | [../../10_Shared/data/spatial_completion/topicvelo_retry/ST_raw_effective_transition_function.py](<../../10_Shared/data/spatial_completion/topicvelo_retry/ST_raw_effective_transition_function.py>) |
| data | [../../10_Shared/data/spatial_completion/topicvelo_retry/ST_raw_native_transition.npz](<../../10_Shared/data/spatial_completion/topicvelo_retry/ST_raw_native_transition.npz>) |
| data | [../../10_Shared/data/spatial_completion/topicvelo_retry/ST_raw_status.json](<../../10_Shared/data/spatial_completion/topicvelo_retry/ST_raw_status.json>) |
| data | [../../10_Shared/data/spatial_completion/topicvelo_retry/ST_small.json](<../../10_Shared/data/spatial_completion/topicvelo_retry/ST_small.json>) |
| data | [../../10_Shared/data/spatial_completion/topicvelo_retry/ST_small_directions.csv](<../../10_Shared/data/spatial_completion/topicvelo_retry/ST_small_directions.csv>) |
| data | [../../10_Shared/data/spatial_completion/topicvelo_retry/ST_small_effective_transition_function.py](<../../10_Shared/data/spatial_completion/topicvelo_retry/ST_small_effective_transition_function.py>) |
| data | [../../10_Shared/data/spatial_completion/topicvelo_retry/ST_small_native_transition.npz](<../../10_Shared/data/spatial_completion/topicvelo_retry/ST_small_native_transition.npz>) |
| data | [../../10_Shared/data/spatial_completion/topicvelo_retry/ST_small_status.json](<../../10_Shared/data/spatial_completion/topicvelo_retry/ST_small_status.json>) |
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
| code | [../../09_Core/shared/spatial_completion/complete_topicvelo.py](<../../09_Core/shared/spatial_completion/complete_topicvelo.py>) |
| code | [../../09_Core/shared/spatial_completion/evaluate_saved_outputs.py](<../../09_Core/shared/spatial_completion/evaluate_saved_outputs.py>) |
| code | [../../09_Core/shared/spatial_topicvelo_retry/evaluate_topicvelo.py](<../../09_Core/shared/spatial_topicvelo_retry/evaluate_topicvelo.py>) |
| code | [../../09_Core/shared/spatial_topicvelo_retry/retry_topicvelo.py](<../../09_Core/shared/spatial_topicvelo_retry/retry_topicvelo.py>) |
| code | [../../09_Core/shared/stability_completion/browser_check.py](<../../09_Core/shared/stability_completion/browser_check.py>) |
| code | [../../09_Core/shared/stability_completion/check_sparse_assembly.py](<../../09_Core/shared/stability_completion/check_sparse_assembly.py>) |
| code | [../../09_Core/shared/stability_completion/complete_topicvelo.py](<../../09_Core/shared/stability_completion/complete_topicvelo.py>) |
| code | [../../09_Core/shared/stability_completion/evaluate_stability.py](<../../09_Core/shared/stability_completion/evaluate_stability.py>) |
| code | [../../09_Core/shared/stability_completion/prepare_publication.py](<../../09_Core/shared/stability_completion/prepare_publication.py>) |
| code | [../../09_Core/shared/stability_completion/render_original_style.py](<../../09_Core/shared/stability_completion/render_original_style.py>) |
| code | [../../09_Core/shared/stability_completion/update_html.py](<../../09_Core/shared/stability_completion/update_html.py>) |
| code | [../../09_Core/shared/stability_completion/validate_delivery.py](<../../09_Core/shared/stability_completion/validate_delivery.py>) |
| code | [../../09_Core/shared/latest/derive_topovelo_historical_20261021.py](<../../09_Core/shared/latest/derive_topovelo_historical_20261021.py>) |
| code | [../../09_Core/shared/latest/topovelo_historical_evaluation.py](<../../09_Core/shared/latest/topovelo_historical_evaluation.py>) |
| code | [../../09_Core/shared/latest/spatial_field_style.py](<../../09_Core/shared/latest/spatial_field_style.py>) |
| code | [code/historical/plots/redraw_ed15_reference.py](<code/historical/plots/redraw_ed15_reference.py>) |

## Code and validation

[Portable checks and workflow boundaries](<../../09_Core/README.md>). Saved-data checks do not rerun models. Historical scripts require their original environments and external inputs.

[Retained figure-specific code](<code>)

[Repository overview](<../../README.md>) · [Release notes](<../../12_Provenance/RELEASE_1021.md>)
