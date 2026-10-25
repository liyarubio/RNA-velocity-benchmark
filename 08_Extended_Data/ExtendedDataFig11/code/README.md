# ED11: individual panel code

Panel a has an annotation-only renderer; its literature schematic remains artwork. Panel b has a formal vector-layout entrypoint that preserves saved paths, aligns all 20 fields to common cell bounds, and sorts them by descending displayed GDC. The GDC score check is separate from the current field panels; current score display is Figure 5e.

| Panel / diagnostic | One-plot file | Scope |
|---|---|---|
| score_check | [ed11_scores.py](score_check/ed11_scores.py) | archived score diagnostic; separate from current field panels |
| b | [ed11b.py](b/ed11b.py) | formal vector field layout; no model or metric rerun |
| a | [ed11a.py](a/ed11a.py) | annotation-only redraw; literature schematic retained as artwork |

Run one file from the repository root:

```bash
python 08_Extended_Data/ExtendedDataFig11/code/score_check/ed11_scores.py --out /tmp/ed11-panels
```

Panel a and the score diagnostic use `requirements-demo.txt` and write PDF/SVG plus a source record. They remain saved-data diagnostics. Panel b reproduces the formal vector layout and additionally requires Poppler (`pdftocairo`) for SVG export. None of these entries refit models.

Historical `legacy/` and `historical/` workflows keep their original numbering and external input requirements; they are provenance sources, not current portable plotting entrypoints.

[Complete panel audit](../../../12_Provenance/PANEL_CODE_INDEX.md) · [Figure overview](../README.md)

For the formal panel run `python 08_Extended_Data/ExtendedDataFig11/code/b/ed11b.py --out demo_output/ed11b`, using [the vector-rendering dependencies](../../../09_Core/visualization/mean_top3/requirements.txt). It writes PDF, SVG, PNG, displayed-order CSV and validation JSON. PhyloVelo stays backward, with adapted GDC +1 and no comparable rank. [Layout validation](../../../12_Provenance/ed11_field_layout_validation.json).
