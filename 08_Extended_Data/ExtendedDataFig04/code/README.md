# ED04: individual panel code

Separate a, b and c scripts use current saved data. Panel c preserves all 840 slots, including missing results; each of its six heatmaps has its own color scale.

| Panel / diagnostic | One-plot file | Scope |
|---|---|---|
| a | [ed4a.py](a/ed4a.py) | saved-data diagnostic; not the formal layout |
| b | [ed4b.py](b/ed4b.py) | saved-data diagnostic; not the formal layout |
| c | [ed4c.py](c/ed4c.py) | saved-data redraw; not the formal layout |

Run one file from the repository root:

```bash
python 08_Extended_Data/ExtendedDataFig04/code/a/ed4a.py --out /tmp/ed04-panels
```

Use `requirements-demo.txt`. Each entrypoint writes exactly one plot in PDF and SVG plus a source record. Shared modules contain readers and plotting primitives; they do not run other panels. These are saved-data checks and do not refit models or replace the formal artwork.

Historical `legacy/` and `historical/` workflows keep their original numbering and external input requirements; they are provenance sources, not current portable plotting entrypoints.

[Complete panel audit](../../../12_Provenance/PANEL_CODE_INDEX.md) · [Figure overview](../README.md)
