# ED08: individual panel code

Separate a, c and d files cover timepoints, reference/method fields and median-cosine diagnostics. Panel b clone examples remain retained artwork because original membership/selection records are incomplete. The d diagnostic shows medians rather than formal boxplots.

| Panel / diagnostic | One-plot file | Scope |
|---|---|---|
| d | [ed08d.py](d/ed08d.py) | median-only diagnostic; formal boxplots retained as artwork |
| c | [ed08c.py](c/ed08c.py) | saved-data redraw; not the formal layout |
| a | [ed8a.py](a/ed8a.py) | saved-data redraw; not the formal layout |

Run one file from the repository root:

```bash
python 08_Extended_Data/ExtendedDataFig08/code/d/ed08d.py --out /tmp/ed08-panels
```

Use `requirements-demo.txt`. Each entrypoint writes exactly one plot in PDF and SVG plus a source record. Shared modules contain readers and plotting primitives; they do not run other panels. These are saved-data checks and do not refit models or replace the formal artwork.

Historical `legacy/` and `historical/` workflows keep their original numbering and external input requirements; they are provenance sources, not current portable plotting entrypoints.

[Complete panel audit](../../../12_Provenance/PANEL_CODE_INDEX.md) · [Figure overview](../README.md)
