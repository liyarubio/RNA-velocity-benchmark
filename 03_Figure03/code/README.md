# Figure 3 code, by panel

Each panel has its own directory. The original combined a/b/c/d/f notebook has been split into a panel-specific `generate_metrics.ipynb` in each corresponding directory. Each includes its required original preparation and metric cells. Numerical panels have separate plotting entrypoints, including b and c.

| Panel | Content | Directory and entrypoint |
|---|---|---|
| [a](a/README.md) | High-dimensional simulation cosines | [a/figure3a.R](a/figure3a.R) |
| [b](b/README.md) | Method-level high- versus low-dimensional comparison | [b/figure3b.R](b/figure3b.R) |
| [c](c/README.md) | Paired high- versus low-dimensional distributions | [c/figure3c.R](c/figure3c.R) |
| [d](d/README.md) | PCA dimension sensitivity (2-30) | [d/figure3d.R](d/figure3d.R) |
| [e](e/README.md) | Metric-definition schematic | [Retained artwork](../figures/Figure03.pdf) |
| [f](f/README.md) | Agreement between evaluation metrics | [f/figure3f.R](f/figure3f.R) |
| [g](g/README.md) | HuBMAP performance | [g/figure3g.R](g/figure3g.R) |
| [h](h/README.md) | Ten conventional datasets | [h/figure3h.R](h/figure3h.R) |

## Saved-data redraws

Install the R packages in [requirements-R.txt](requirements-R.txt), then run from the repository root:

```bash
FIGURE3_OUT=/tmp/figure3-panels Rscript --vanilla 03_Figure03/code/render_all.R
```

The runner calls the seven individual numerical panel scripts. Each can also run alone. Without `FIGURE3_OUT`, outputs go to `demo_output/figure3`. Panel e is manually composed artwork. The retained formal figure is unchanged.

The renderers were recovered from `Manuscript_2607/Figure3_sources`; all nine numerical input CSVs are byte-identical to the current Manuscript_1021 inputs. R syntax, path bindings and preserved computation blocks were checked. Full R plotting execution was not completed because available environments lack required packages. These checks do not certify identical formal rendering or a model rerun.

## Upstream computation

See each panel README for its metric sources and [shared/README.md](shared/README.md) for reusable helpers and the auxiliary PCA30 calculation. The original notebook code cells are preserved; historical outputs and execution counts were cleared during recovery. These generators require external H5ADs, original paths and method environments. Their method coverage may differ from current panels; current CSVs and the workbook define released values.

[Shared dyngen simulation and evaluation](../../09_Core/simulation/dyngen/README.md) · [Earlier visualization implementations](legacy/Visualization/Vis_summary/)

## Independently check saved data

```bash
python 03_Figure03/code/validate_panel_data.py --out /tmp/figure3-panel-validation.json
```

This recalculates 435 panel-d dimension summaries, 798 panel-f aligned units and their pairwise-complete Pearson correlations. Missing values remain missing.

[Recovery hashes and scope](../../12_Provenance/figure3_code_recovery.json) · [Figure overview](../README.md)
