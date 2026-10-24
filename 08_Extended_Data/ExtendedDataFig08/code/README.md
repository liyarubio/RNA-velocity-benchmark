# ED08: individual panel code

Separate a, c and d files cover timepoints, reference/method fields and median-cosine diagnostics. Panel b clone examples remain retained artwork because original membership/selection records are incomplete. The d diagnostic shows medians rather than formal boxplots.

| Panel / diagnostic | One-plot file | Scope |
|---|---|---|
| d | [ed08d.py](d/ed08d.py) | median-only diagnostic; formal boxplots retained as artwork |
| c | [ed08c.py](c/ed08c.py) | saved-data redraw; not the formal layout |
| a | [ed8a.py](a/ed8a.py) | formal panel with ED8c cell-type colors; original point geometry retained |

Run one file from the repository root:

```bash
python 08_Extended_Data/ExtendedDataFig08/code/d/ed08d.py --out /tmp/ed08-panels
```

Panels c and d use `requirements-demo.txt` and write one diagnostic plot in PDF and SVG plus a source record. They do not refit models or replace the formal artwork. Panel a reproduces the formal panel using [the shared vector-rendering dependencies](../../../09_Core/visualization/mean_top3/requirements.txt); add `--full-figure` to export the complete artwork page. Shared modules do not execute panels on import.

Historical `legacy/` and `historical/` workflows keep their original numbering and external input requirements; they are provenance sources, not current portable plotting entrypoints.

[Complete panel audit](../../../12_Provenance/PANEL_CODE_INDEX.md) · [Figure overview](../README.md)

Panel a writes PDF, SVG, PNG and validation JSON. Its seven shared cell types match the formal panel c palette exactly. The six Mast and two Ccr7_DC cells present only at Day 5 are displayed together in gray as Other, with their original annotations retained in the source CSV. The panel-c saved-data entry now uses the formal palette too; its published artwork is unchanged. [Annotation provenance](../../../12_Provenance/ed08a_celltype_annotations.json).
