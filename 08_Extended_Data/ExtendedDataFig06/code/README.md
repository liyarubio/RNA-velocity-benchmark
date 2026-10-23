# ED06: individual panel code

Separate a and b files check HuBMAP and classical-dataset GDC. These diagnostics use raw-order heatmaps; they do not reconstruct the formal clustering, organ strip or bubble encodings.

| Panel / diagnostic | One-plot file | Scope |
|---|---|---|
| a | [ed6a.py](a/ed6a.py) | raw-order heatmap diagnostic; does not reproduce formal clustering or organ strip |
| b | [ed6b.py](b/ed6b.py) | saved-data redraw; not the formal layout |

Run one file from the repository root:

```bash
python 08_Extended_Data/ExtendedDataFig06/code/a/ed6a.py --out /tmp/ed06-panels
```

Use `requirements-demo.txt`. Each entrypoint writes exactly one plot in PDF and SVG plus a source record. Shared modules contain readers and plotting primitives; they do not run other panels. These are saved-data checks and do not refit models or replace the formal artwork.

Historical `legacy/` and `historical/` workflows keep their original numbering and external input requirements; they are provenance sources, not current portable plotting entrypoints.

[Complete panel audit](../../../12_Provenance/PANEL_CODE_INDEX.md) · [Figure overview](../README.md)
