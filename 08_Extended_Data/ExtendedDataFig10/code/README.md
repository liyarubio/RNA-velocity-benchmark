# ED10: individual panel code

Separate a and b files draw cell-type medians and all paired spot maps. The former across-method median diagnostic remains separately named; the boxplot is current Figure 5c.

| Panel / diagnostic | One-plot file | Scope |
|---|---|---|
| median_check | [median_check.py](median_check/median_check.py) | median-only diagnostic; formal a-b are separate |
| a | [ed10a.py](a/ed10a.py) | saved-data redraw; not the formal layout |
| b | [ed10b.py](b/ed10b.py) | saved-data redraw; not the formal layout |

Run one file from the repository root:

```bash
python 08_Extended_Data/ExtendedDataFig10/code/median_check/median_check.py --out /tmp/ed10-panels
```

Use `requirements-demo.txt`. Each entrypoint writes exactly one plot in PDF and SVG plus a source record. Shared modules contain readers and plotting primitives; they do not run other panels. These are saved-data checks and do not refit models or replace the formal artwork.

Historical `legacy/` and `historical/` workflows keep their original numbering and external input requirements; they are provenance sources, not current portable plotting entrypoints.

[Complete panel audit](../../../12_Provenance/PANEL_CODE_INDEX.md) · [Figure overview](../README.md)
