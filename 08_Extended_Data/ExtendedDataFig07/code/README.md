# ED07: individual panel code

Separate a, b and c files cover annotation, reference/method fields and median-cosine diagnostics. The c diagnostic shows medians, not the complete formal boxplot distribution.

| Panel / diagnostic | One-plot file | Scope |
|---|---|---|
| c | [ed07c.py](c/ed07c.py) | median-only diagnostic; formal boxplots retained as artwork |
| a | [ed7a.py](a/ed7a.py) | saved-data redraw; not the formal layout |
| b | [ed07b.py](b/ed07b.py) | saved-data redraw; not the formal layout |

Run one file from the repository root:

```bash
python 08_Extended_Data/ExtendedDataFig07/code/c/ed07c.py --out /tmp/ed07-panels
```

Use `requirements-demo.txt`. Each entrypoint writes exactly one plot in PDF and SVG plus a source record. Shared modules contain readers and plotting primitives; they do not run other panels. These are saved-data checks and do not refit models or replace the formal artwork.

Historical `legacy/` and `historical/` workflows keep their original numbering and external input requirements; they are provenance sources, not current portable plotting entrypoints.

[Complete panel audit](../../../12_Provenance/PANEL_CODE_INDEX.md) · [Figure overview](../README.md)
