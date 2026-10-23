# ED07 — scNT metabolic-labeling reference

Current release: **Manuscript_1021**, including updates through **22 October 2026**. Historical numerical source ID: **ED13**.

1,947 reference cells and 19 configurations; projected cosines were recomputed per cell. Preprocessing differs across pipelines.

![ED07](<figures/Extended_Data_Fig07.png>)

## Individual panel code

[One plot per file: current entrypoints and limits](code/README.md).

## Current files

| Role | File |
|---|---|
| results | [figures/Extended_Data_Fig07.pdf](<figures/Extended_Data_Fig07.pdf>) |
| results | [figures/Extended_Data_Fig07.png](<figures/Extended_Data_Fig07.png>) |
| results | [figures/Extended_Data_Fig07.svg](<figures/Extended_Data_Fig07.svg>) |
| results | [figures/legends/Extended_Data_Fig07_legend.pdf](<figures/legends/Extended_Data_Fig07_legend.pdf>) |
| data | [data/source_data/Extended_Data_Fig13.xlsx](<data/source_data/Extended_Data_Fig13.xlsx>) |
| data | [data/reference_restoration/scNT_projected_directions.csv](<data/reference_restoration/scNT_projected_directions.csv>) |
| code | [../../09_Core/shared/plots/restore_reference_figures.py](<../../09_Core/shared/plots/restore_reference_figures.py>) |
| code | [../../09_Core/reference_evaluation/scLT_scNT_velogif/evaluate_reference.py](<../../09_Core/reference_evaluation/scLT_scNT_velogif/evaluate_reference.py>) |
| code | [../../09_Core/reference_evaluation/scLT_scNT_velogif/finalize.py](<../../09_Core/reference_evaluation/scLT_scNT_velogif/finalize.py>) |
| code | [../../09_Core/reference_evaluation/scLT_scNT_velogif/prepare.py](<../../09_Core/reference_evaluation/scLT_scNT_velogif/prepare.py>) |
| code | [../../09_Core/reference_evaluation/scLT_scNT_velogif/run_all.py](<../../09_Core/reference_evaluation/scLT_scNT_velogif/run_all.py>) |
| code | [../../09_Core/reference_evaluation/scLT_scNT_velogif/test_reference.py](<../../09_Core/reference_evaluation/scLT_scNT_velogif/test_reference.py>) |
| code | [../../09_Core/shared/restore_tfvelo_projection.py](<../../09_Core/shared/restore_tfvelo_projection.py>) |
| data | [data/reference_restoration/scNT_TFvelo_restored_projection.csv](<data/reference_restoration/scNT_TFvelo_restored_projection.csv>) |
| data | [data/reference_restoration/scNT_method_summary.csv](<data/reference_restoration/scNT_method_summary.csv>) |
| data | [data/reference_restoration/scNT_per_cell_cosine.csv](<data/reference_restoration/scNT_per_cell_cosine.csv>) |
| data | [data/reference_restoration/scNT_reference_cells.csv](<data/reference_restoration/scNT_reference_cells.csv>) |

## Code and validation

[Portable checks and workflow boundaries](<../../09_Core/README.md>). Saved-data checks do not rerun models. Historical scripts require their original environments and external inputs.

[Retained figure-specific code](<code>)

[Repository overview](<../../README.md>) · [Release notes](<../../12_Provenance/RELEASE_1021.md>)
