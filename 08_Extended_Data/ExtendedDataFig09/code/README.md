# ED09: individual panel code

The current field panel has a separate 20-method saved-vector renderer. The archived GDC score check is separately named; it is not a current lettered panel.

| Panel / diagnostic | One-plot file | Scope |
|---|---|---|
| score_check | [ed09_scores.py](score_check/ed09_scores.py) | archived score diagnostic; separate from current field panels |
| fields | [ed9_fields.py](fields/ed9_fields.py) | saved-vector redraw; not the formal layout |

Run one file from the repository root:

```bash
python 08_Extended_Data/ExtendedDataFig09/code/score_check/ed09_scores.py --out /tmp/ed09-panels
```

Use `requirements-demo.txt`. Each entrypoint writes exactly one plot in PDF and SVG plus a source record. Shared modules contain readers and plotting primitives; they do not run other panels. These are saved-data checks and do not refit models or replace the formal artwork.

Historical `legacy/` and `historical/` workflows keep their original numbering and external input requirements; they are provenance sources, not current portable plotting entrypoints.

[Complete panel audit](../../../12_Provenance/PANEL_CODE_INDEX.md) · [Figure overview](../README.md)
