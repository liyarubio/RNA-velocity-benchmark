# ED04 — Simulation cosine comparisons and full noise series

Current release: **Manuscript_1021**, including figure presentation updates through **23 October 2026**. Historical numerical source ID: **ED05**.

Noise panel c contains 840 slots (822 finite, 18 missing), with 84 cross-noise means. STT forward branching has 1/10 coverage; other groups have 10/10. Conditions are not independent replicates.

![ED04](<figures/Extended_Data_Fig04.png>)

## Individual panel code

[One plot per file: current entrypoints and limits](code/README.md).

## Mean-only summary update

The simulation summary columns now display Mean only, with ranks 1–3 aligned inside the bars. Numerical inputs and normalization are unchanged. Companion caption text is retained for author revision. See the panel code guide for independent formal rendering commands.

## Current files

| Role | File |
|---|---|
| results | [figures/Extended_Data_Fig04.pdf](<figures/Extended_Data_Fig04.pdf>) |
| results | [figures/Extended_Data_Fig04.png](<figures/Extended_Data_Fig04.png>) |
| results | [figures/Extended_Data_Fig04.svg](<figures/Extended_Data_Fig04.svg>) |
| results | [figures/legends/Extended_Data_Fig04_legend.pdf](<figures/legends/Extended_Data_Fig04_legend.pdf>) |
| data | [data/source_data/Extended_Data_Fig05.xlsx](<data/source_data/Extended_Data_Fig05.xlsx>) |
| data | [data/cellDancer_current_30/README.txt](<data/cellDancer_current_30/README.txt>) |
| data | [data/cellDancer_current_30/dataset_inventory.csv](<data/cellDancer_current_30/dataset_inventory.csv>) |
| data | [data/cellDancer_current_30/dataset_summary.csv](<data/cellDancer_current_30/dataset_summary.csv>) |
| data | [data/cellDancer_current_30/evaluation_protocol.json](<data/cellDancer_current_30/evaluation_protocol.json>) |
| data | [data/cellDancer_current_30/method_coverage.csv](<data/cellDancer_current_30/method_coverage.csv>) |
| code | [../../09_Core/shared/plots/redraw_ed5.py](<../../09_Core/shared/plots/redraw_ed5.py>) |
| code | [../../09_Core/shared/plots/redraw_extended.py](<../../09_Core/shared/plots/redraw_extended.py>) |
| data | [data/cellDancer_noise_20261020/README.txt](<data/cellDancer_noise_20261020/README.txt>) |
| data | [data/cellDancer_noise_20261020/dataset_inventory.csv](<data/cellDancer_noise_20261020/dataset_inventory.csv>) |
| data | [data/cellDancer_noise_20261020/method_coverage.csv](<data/cellDancer_noise_20261020/method_coverage.csv>) |
| data | [data/cellDancer_noise_20261020/noise_means.csv](<data/cellDancer_noise_20261020/noise_means.csv>) |
| data | [data/cellDancer_noise_20261020/panel_color_scales.csv](<data/cellDancer_noise_20261020/panel_color_scales.csv>) |
| data | [data/cellDancer_noise_20261020/plot_data.csv](<data/cellDancer_noise_20261020/plot_data.csv>) |
| results | [figures/Extended_Data_Fig04_noise.pdf](<figures/Extended_Data_Fig04_noise.pdf>) |
| results | [figures/Extended_Data_Fig04_noise.svg](<figures/Extended_Data_Fig04_noise.svg>) |
| results | [figures/Extended_Data_Fig04_noise.png](<figures/Extended_Data_Fig04_noise.png>) |
| code | [code/historical/cellDancer_evaluation/evaluate.py](<code/historical/cellDancer_evaluation/evaluate.py>) |
| code | [code/historical/cellDancer_evaluation/validate_data.py](<code/historical/cellDancer_evaluation/validate_data.py>) |

## Code and validation

[Portable checks and workflow boundaries](<../../09_Core/README.md>). Saved-data checks do not rerun models. Historical scripts require their original environments and external inputs.

[Retained figure-specific code](<code>)

[Repository overview](<../../README.md>) · [Release notes](<../../12_Provenance/RELEASE_1021.md>)
