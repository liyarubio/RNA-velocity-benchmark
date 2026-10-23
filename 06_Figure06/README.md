# Figure06 — Runtime, read-depth sensitivity, frozen composite and VeloGIF

Current release: **Manuscript_1021**, including updates through **22 October 2026**. Historical numerical source ID: **Figure06**.

Panels a/b now contain runtime/read-depth sensitivity; c is the unchanged frozen composite (0.6 accuracy rank + 0.4 scalability rank; stability excluded); d is VeloGIF. Matched-protocol limitations remain.

![Figure06](<figures/Figure06.png>)

## Individual panel code

[One plot per file: current entrypoints and limits](code/README.md).

## Current files

| Role | File |
|---|---|
| results | [figures/Figure06.pdf](<figures/Figure06.pdf>) |
| results | [figures/Figure06.png](<figures/Figure06.png>) |
| results | [figures/legends/Figure06_legend.pdf](<figures/legends/Figure06_legend.pdf>) |
| data | [data/source_data/Figure06.xlsx](<data/source_data/Figure06.xlsx>) |
| code | [code/b/historical/compute_metrics.py](<code/b/historical/compute_metrics.py>) |
| code | [code/c/historical_composite.py](<code/c/historical_composite.py>) |

## Code and validation

[Portable checks and workflow boundaries](<../09_Core/README.md>). Saved-data checks do not rerun models. Historical scripts require their original environments and external inputs.

[Retained figure-specific code](<code>)

[Repository overview](<../README.md>) · [Release notes](<../12_Provenance/RELEASE_1021.md>)
